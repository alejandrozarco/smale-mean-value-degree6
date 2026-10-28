#!/usr/bin/env python3
"""verify_chunks.py — fail-closed, input-bound verification of a stream of subdivision-tree records with arbcheck2.
Replaces check_stream.py for certification (findings 1-7 of the review of the d = 7 chain).

usage: zstd -dc ARCHIVE | python3 c/verify_chunks.py LABEL CANON_TASKS CANON_TASKS_SHA256 CFG16 ARBCHECK2 \
           ARBCHECK2_SHA256 NTHREADS CHUNK_BYTES STREAM_SHA256 WORKDIR [D]

1. The canonical tasks file must have the given sha256, and its header must be
   "# smale_bb_v3 tasks cfg=CFG16 <cfg string> ntask=N half_width=0x1p-2" (s = 4, the run's grid; check_done.py
   verifies that this file lists exactly the canonical grid boxes). The arbcheck2 binary must have the given sha256.
2. The stream (stdin) is split into record-aligned chunks of about CHUNK_BYTES. For each chunk:
   - its bytes are written to WORKDIR/LABEL.NNNN.tree and hashed (chunk_sha256);
   - its subset tasks file is built from the canonical file: the canonical header line, plus the canonical line of
     each id in the chunk;
   - arbcheck2 is run with fixed arguments (-t NTHREADS --quiet D subset chunk);
   - success requires exit status 0, exactly one SUMMARY line, status=PASS, fail=0, missing=0, dup=0, unknown=0,
     badheader=0, tail_incomplete=0, failed_records=0, and records = the number of records in the chunk;
   - a receipt WORKDIR/LABEL.NNNN.receipt.json is written, binding chunk_sha256, chunk byte count, the ids file
     ("id nodes record_sha256", its sha256), the subset tasks sha256, the canonical tasks sha256, cfg16, the arbcheck2
     sha256, the argument list, the SUMMARY line and the exit status.
   - a receipt is reused only if it has PASS and the same chunk_sha256, subset tasks sha256, canonical tasks sha256,
     arbcheck2 sha256 and arguments; any other existing receipt is an error.
3. At the end the sha256 of the whole stream must equal STREAM_SHA256; then LABEL.stream.json is written.
Any failure exits non-zero, and no success receipt is written for that chunk or stream."""
import hashlib, json, os, re, struct, subprocess, sys, time

def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''): h.update(b)
    return h.hexdigest()

def die(msg):
    print('VERIFY FAILED:', msg, file=sys.stderr); sys.exit(2)

(label, canon, canon_sha, cfg16, ab, ab_sha, nth, lim, stream_sha, wd) = sys.argv[1:11]
d = sys.argv[11] if len(sys.argv) > 11 else '7'
lim = int(lim); os.makedirs(wd, exist_ok=True)
if not re.fullmatch(r'[A-Za-z0-9_]+', label): die('bad label')
if sha_file(canon) != canon_sha: die('canonical tasks sha256 mismatch')
if sha_file(ab) != ab_sha: die('arbcheck2 sha256 mismatch')
lines = open(canon).read().split('\n')
header = lines[0]
if not re.fullmatch(r'# smale_bb_v3 tasks cfg=' + cfg16 + r' smale_bb_v3;src=[0-9a-f]{64};d=' + d +
                    r';r_excl=[0-9a-f]{16};excl_mode=2;s=4;seed=1 ntask=\d+ half_width=0x1p-2', header):
    die('canonical tasks header does not match the configuration')
tl = {}
for l in lines[1:]:
    if l.strip():
        i = int(l.split(None, 1)[0])
        if i in tl: die(f'canonical tasks: duplicate id {i}')
        tl[i] = l
args = ['-t', str(nth), '--quiet', d]
inp = sys.stdin.buffer; hs = hashlib.sha256(); k = 0

def verify_chunk(name, path, recs):
    csha = sha_file(path); cbytes = os.path.getsize(path)
    idsp = os.path.join(wd, name + '.ids'); idst = idsp + '.tmp'      # committed only with a passing receipt
    with open(idst, 'w') as f:
        for i, n, rs in recs: f.write(f'{i} {n} {rs}\n')
    tp = os.path.join(wd, name + '.tasks')
    with open(tp, 'w') as f:
        f.write(header + '\n')
        for i in sorted(i for i, _, _ in recs):
            if i not in tl: die(f'{name}: id {i} not in canonical tasks')
            f.write(tl[i] + '\n')
    tsha = sha_file(tp)
    rp = os.path.join(wd, name + '.receipt.json')
    if os.path.exists(rp):
        r = json.load(open(rp))
        if (r.get('ok') and r['chunk_sha256'] == csha and r['subset_tasks_sha256'] == tsha and
                r['canonical_tasks_sha256'] == canon_sha and r['arbcheck2_sha256'] == ab_sha and r['args'] == args and
                r['ids_sha256'] == sha_file(idst) and os.path.exists(idsp) and sha_file(idsp) == r['ids_sha256']):
            os.remove(path); os.remove(tp); os.remove(idst); return
        os.remove(idst)
        die(f'{name}: existing receipt does not match this chunk (refusing to reuse or overwrite)')
    t0 = time.time()
    p = subprocess.run(['nice', '-n', '19', ab] + args + [tp, path], capture_output=True, text=True)
    summ = [l for l in p.stdout.split('\n') if l.startswith('SUMMARY')]
    with open(os.path.join(wd, name + '.err'), 'w') as f: f.write(p.stderr[-500000:])
    ok = p.returncode == 0 and len(summ) == 1
    g = dict(re.findall(r'(\w+)=(\S+)', summ[0])) if len(summ) == 1 else {}
    for kk, v in (('status', 'PASS'), ('fail', '0'), ('missing', '0'), ('dup', '0'), ('unknown', '0'),
                  ('badheader', '0'), ('tail_incomplete', '0'), ('failed_records', '0'), ('records', str(len(recs)))):
        ok = ok and g.get(kk) == v
    r = dict(label=label, chunk=name, chunk_sha256=csha, chunk_bytes=cbytes, n_records=len(recs),
             ids_sha256=sha_file(idst), subset_tasks_sha256=tsha, canonical_tasks_sha256=canon_sha, cfg16=cfg16,
             arbcheck2_sha256=ab_sha, args=args, summary=summ[0] if summ else None, returncode=p.returncode,
             seconds=round(time.time() - t0, 1), ok=ok)
    if not ok:
        json.dump(r, open(rp + '.FAILED', 'w'), indent=1); die(f'{name}: arbcheck2 did not pass: {r["summary"]}')
    os.replace(idst, idsp)
    json.dump(r, open(rp, 'w'), indent=1)
    os.remove(path); os.remove(tp)

while True:
    name = f'{label}.{k:04d}'; path = os.path.join(wd, name + '.tree')
    recs = []; raw = 0
    with open(path, 'wb') as out:
        while raw < lim:
            h = inp.read(14)
            if not h: break
            if len(h) != 14 or h[:2] != b'TK': die('bad record header (stream not record-aligned)')
            i, n = struct.unpack('<IQ', h[2:]); body = inp.read(n)
            if len(body) != n: die('truncated record')
            hs.update(h); hs.update(body); out.write(h); out.write(body)
            recs.append((i, n, hashlib.sha256(h + body).hexdigest())); raw += 14 + n
    if not recs:
        os.remove(path); break
    verify_chunk(name, path, recs); k += 1
if hs.hexdigest() != stream_sha: die(f'stream sha256 {hs.hexdigest()} != expected {stream_sha}')
json.dump(dict(label=label, chunks=k, stream_sha256=stream_sha, canonical_tasks_sha256=canon_sha,
               arbcheck2_sha256=ab_sha, ok=True), open(os.path.join(wd, label + '.stream.json'), 'w'), indent=1)
print(f'{label}: {k} chunks verified, stream sha256 {stream_sha}')
