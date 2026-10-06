#!/usr/bin/env python3
"""merge_done.py — merge the checkpoint files (.done) of several partial runs of ONE configuration (smale_bb_v3 with
SMALE_PART on several machines) into one .done file, for c/check_done.py.
usage: python3 c/merge_done.py OUT.done IN1.done IN2.done ...
Every input line must have a valid checksum; all headers must be identical; a torn trailing line is ignored (as the
run program would truncate it). A task id present in several inputs must have identical counts (the computation is
deterministic); one copy is kept. The output is the header followed by the records sorted by id. The records keep
their own checksums, so check_done.py validates every line of the merged file exactly as for a single run."""
import hashlib, re, sys

def payload(line):
    if len(line) < 18 or line[-17] != ' ':
        return None
    p, chk = line[:-17], line[-16:]
    return p if hashlib.sha256(p.encode()).hexdigest()[:16] == chk else None

out, ins = sys.argv[1], sys.argv[2:]
header, recs = None, {}
for fn in ins:
    raw = open(fn, 'rb').read().decode('ascii')
    lines = raw.split('\n')
    if not raw.endswith('\n'):
        print(f'{fn}: ignoring torn trailing line'); lines = lines[:-1]
    lines = [l for l in lines if l]
    for k, l in enumerate(lines):
        p = payload(l)
        if p is None:
            sys.exit(f'{fn} line {k + 1}: bad checksum')
        if k == 0:
            if header is None: header = l
            elif l != header: sys.exit(f'{fn}: header differs from {ins[0]}')
            continue
        # full record format (as check_done.py): id, 8 counts, maxdepth, seconds, volume, cfg hash of the header
        m = re.fullmatch(r'R (\d+) ((?:\d+ ){8})(\d+) (\d+\.\d{3}) (\S+) cfg=([0-9a-f]{16})', p)
        if not m or not re.fullmatch(r'H \S+ cfg=' + m.group(6) + r' .*', payload(header)):
            sys.exit(f'{fn} line {k + 1}: malformed record or record from another configuration')
        t = int(m.group(1))
        key = (m.group(2), m.group(3), m.group(5))          # counts, maxdepth, volume: deterministic
        if t in recs:
            if recs[t][1] != key:
                sys.exit(f'task {t}: different deterministic fields in two inputs')
            continue
        recs[t] = (l, key)
with open(out, 'w') as f:
    f.write(header + '\n')
    for t in sorted(recs):
        f.write(recs[t][0] + '\n')
print(f'merged {len(ins)} files: {len(recs)} distinct task records -> {out}')
