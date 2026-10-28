#!/usr/bin/env python3
"""cover_check2.py — the input-bound verification (c/verify_chunks.py) covers the whole run exactly once.
usage: python3 c/cover_check2.py MERGED.done CANON_TASKS_SHA256 ARBCHECK2_SHA256[,SHA256...] CFG16 DIR...
(several arbcheck2 sha256 values when builds on several machines are used; each needs its build record)
Collects every *.receipt.json and *.stream.json under the DIRs and checks:
- every receipt has ok=true, a SUMMARY with status=PASS, fail=0, missing=0 and records = n_records, the given canonical
  tasks sha256, arbcheck2 sha256 and cfg16, and its ids file exists with the recorded sha256 (recomputed);
- chunk names are unique across all DIRs; every label has a stream receipt with ok=true and the same hashes, and its
  chunk indices are 0..chunks-1 with no gap;
- the ids (from the ids files) are pairwise disjoint and equal the set of task ids of MERGED.done, and for every id the
  tree record's node count equals `processed` in MERGED.done;
- prints a canonical digest: sha256 over the lines "id record_sha256" in increasing id order.
Prints COVER COMPLETE only if everything holds. *.FAILED receipts are reported."""
import glob, hashlib, json, os, re, sys

def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''): h.update(b)
    return h.hexdigest()

done, tsha, asha, cfg16 = sys.argv[1:5]; dirs = sys.argv[5:]; asha = set(asha.split(','))
proc = {}
for l in open(done):
    m = re.match(r'R (\d+) (\d+) ', l)
    if m: proc[int(m.group(1))] = int(m.group(2))
probs = []; chunks = {}; streams = {}
for dd in dirs:
    for f in glob.glob(os.path.join(dd, '**', '*.receipt.json*'), recursive=True):
        if f.endswith('.FAILED'): probs.append(f'failed receipt present: {f}'); continue
        r = json.load(open(f)); c = r['chunk']
        if c in chunks: probs.append(f'duplicate chunk name {c}'); continue
        chunks[c] = (r, os.path.dirname(f))
    for f in glob.glob(os.path.join(dd, '**', '*.stream.json'), recursive=True):
        r = json.load(open(f))
        if r['label'] in streams: probs.append(f'duplicate stream label {r["label"]}')
        streams[r['label']] = r
seen = {}; rec_sha = {}
for c, (r, base) in sorted(chunks.items()):
    g = dict(re.findall(r'(\w+)=(\S+)', r.get('summary') or ''))
    if not (r.get('ok') and r.get('returncode') == 0 and g.get('status') == 'PASS' and g.get('fail') == '0' and
            g.get('missing') == '0' and g.get('records') == str(r['n_records'])):
        probs.append(f'{c}: not a passing receipt')
    if r['canonical_tasks_sha256'] != tsha or r['arbcheck2_sha256'] not in asha or r['cfg16'] != cfg16:
        probs.append(f'{c}: receipt bound to other inputs')
    idsp = os.path.join(base, c + '.ids')
    if not os.path.exists(idsp) or sha_file(idsp) != r['ids_sha256']:
        probs.append(f'{c}: ids file missing or changed'); continue
    n = 0
    for l in open(idsp):
        i, nodes, rs = l.split(); i = int(i); nodes = int(nodes); n += 1
        if i in seen: probs.append(f'id {i} in {seen[i]} and {c}')
        seen[i] = c; rec_sha[i] = rs
        if proc.get(i) != nodes: probs.append(f'id {i}: tree nodes {nodes} != processed {proc.get(i)}')
    if n != r['n_records']: probs.append(f'{c}: ids file has {n} lines, receipt says {r["n_records"]}')
labels = {}
for c in chunks:
    lab, idx = c.rsplit('.', 1); labels.setdefault(lab, []).append(int(idx))
for lab, idxs in labels.items():
    s = streams.get(lab)
    if not s or not s.get('ok') or s['canonical_tasks_sha256'] != tsha or s['arbcheck2_sha256'] not in asha:
        probs.append(f'label {lab}: no valid stream receipt'); continue
    if sorted(idxs) != list(range(s['chunks'])): probs.append(f'label {lab}: chunks {sorted(idxs)} != 0..{s["chunks"]-1}')
for lab in streams:
    if lab not in labels: probs.append(f'stream {lab} has no chunk receipts')
missing = set(proc) - set(seen); extra = set(seen) - set(proc)
if missing: probs.append(f'{len(missing)} run ids without a verified tree record, e.g. {sorted(missing)[:5]}')
if extra: probs.append(f'{len(extra)} verified ids not in the run, e.g. {sorted(extra)[:5]}')
dig = hashlib.sha256(''.join(f'{i} {rec_sha[i]}\n' for i in sorted(rec_sha)).encode()).hexdigest()
for p in probs[:30]: print('PROBLEM:', p)
print(f'{len(chunks)} chunk receipts, {len(streams)} streams, {len(seen)} ids verified, {len(proc)} ids in the run')
print(f'canonical record digest (sha256 of "id record_sha256" lines in id order): {dig}')
print('COVER COMPLETE' if not probs else f'COVER FAILED ({len(probs)} problems)')
sys.exit(1 if probs else 0)
