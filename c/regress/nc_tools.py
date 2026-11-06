#!/usr/bin/env python3
"""Tree-record tampering helpers for the negative controls of c/regress/nc_verify_chain.sh.
usage: nc_tools.py CMD IN OUT [args]
  largest IN            -> prints the id and size of the largest record
  rootleaf IN OUT ID    -> replace record ID's body with a single F leaf (code 0)
  unres IN OUT ID       -> set the first leaf of record ID to code 6 (unresolved)
  drop IN OUT ID        -> remove record ID
  dup IN OUT ID         -> append a copy of record ID at the end
  tear IN OUT           -> cut the stream in the middle of the last record
  swapcols IN_TASKS OUT_TASKS A B -> swap the coordinates of tasks A and B (header kept)
  firstx IN             -> prints the id of the smallest record holding an excluded leaf (a box near an equality point)
  pickx OUT K MAXBYTES  -> stdin stream: write the K smallest records holding an excluded leaf among the first MAXBYTES
  editids IN OUT        -> change one character of the first line of an ids file
"""
import struct, sys, hashlib

def read(path):
    data = open(path, 'rb').read(); off = 0; recs = []
    while off < len(data):
        h = data[off:off + 14]; assert h[:2] == b'TK'
        i, n = struct.unpack('<IQ', h[2:]); recs.append((i, data[off + 14:off + 14 + n])); off += 14 + n
    return recs

def write(path, recs):
    with open(path, 'wb') as f:
        for i, b in recs: f.write(b'TK' + struct.pack('<IQ', i, len(b)) + b)

cmd = sys.argv[1]
if cmd == 'largest':
    recs = read(sys.argv[2]); i, b = max(recs, key=lambda r: len(r[1])); print(i, len(b))
elif cmd == 'firstx':
    print(min((len(b), i) for i, b in read(sys.argv[2]) if b.count(5))[1])
elif cmd == 'pickx':
    out, k, lim = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]); inp = sys.stdin.buffer; got = []; seen = 0
    while seen < lim:
        h = inp.read(14)
        if len(h) < 14: break
        i, n = struct.unpack('<IQ', h[2:]); b = inp.read(n); seen += 14 + n
        if b.count(5): got.append((len(b), i, b))
    write(out, [(i, b) for _, i, b in sorted(got)[:k]])
elif cmd == 'editids':
    d = open(sys.argv[2]).read(); open(sys.argv[3], 'w').write(('1' if d[0] != '1' else '2') + d[1:])
elif cmd == 'rootleaf':
    recs = read(sys.argv[2]); t = int(sys.argv[4])
    write(sys.argv[3], [(i, b'\x00' if i == t else b) for i, b in recs])
elif cmd == 'unres':
    recs = read(sys.argv[2]); t = int(sys.argv[4]); out = []
    for i, b in recs:
        if i == t:
            k = next(j for j, c in enumerate(b) if c < 16); b = b[:k] + b'\x06' + b[k + 1:]
        out.append((i, b))
    write(sys.argv[3], out)
elif cmd == 'drop':
    recs = read(sys.argv[2]); t = int(sys.argv[4]); write(sys.argv[3], [r for r in recs if r[0] != t])
elif cmd == 'dup':
    recs = read(sys.argv[2]); t = int(sys.argv[4]); write(sys.argv[3], recs + [r for r in recs if r[0] == t])
elif cmd == 'tear':
    data = open(sys.argv[2], 'rb').read(); open(sys.argv[3], 'wb').write(data[:-7])
elif cmd == 'swapcols':
    lines = open(sys.argv[2]).read().split('\n'); a, b = int(sys.argv[4]), int(sys.argv[5])
    la, lb = lines[1 + a].split(None, 1)[1], lines[1 + b].split(None, 1)[1]
    lines[1 + a] = f'{a} {lb}'; lines[1 + b] = f'{b} {la}'
    open(sys.argv[3], 'w').write('\n'.join(lines))
else:
    sys.exit(__doc__)
