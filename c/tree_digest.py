"""tree_digest.py — canonical sha256 of a subdivision tree file written by c/export_tree.c.
The exporter writes records in thread-completion order, so the raw file hash is not reproducible; this digest is the
sha256 of the records ("TK" | uint32 id | uint64 nbytes | codes) concatenated in increasing id order.
usage: python3 c/tree_digest.py FILE.tree [--write CANONICAL.tree]
  --write also writes the records in increasing id order (the canonical file; its sha256 equals the digest).
       python3 c/tree_digest.py --records FILE.tree|-  [more files ...]
  the canonical record digest of c/cover_check2.py: sha256 over the lines "id record_sha256" in increasing id order,
  record_sha256 = sha256("TK" | id | nbytes | codes). Streams the input (- = standard input), so a tree split over
  several files or streams in any order gives the same digest; each id must occur once."""
import hashlib, struct, sys

if sys.argv[1] == '--records':
    rec = {}
    for path in sys.argv[2:]:
        f = sys.stdin.buffer if path == '-' else open(path, 'rb')
        while True:
            h = f.read(14)
            if not h:
                break
            assert len(h) == 14 and h[:2] == b'TK', f'bad record header in {path}'
            rid, nb = struct.unpack('<IQ', h[2:])
            r = hashlib.sha256(h); left = nb
            while left:
                b = f.read(min(left, 1 << 24))
                assert b, f'truncated record {rid} in {path}'
                r.update(b); left -= len(b)
            assert rid not in rec, f'duplicate record id {rid}'
            rec[rid] = r.hexdigest()
    dig = hashlib.sha256(''.join(f'{i} {rec[i]}\n' for i in sorted(rec)).encode()).hexdigest()
    print(f'{dig}  {len(rec)} records (canonical record digest)')
    sys.exit(0)

path = sys.argv[1]
out = open(sys.argv[3], 'wb') if len(sys.argv) > 3 and sys.argv[2] == '--write' else None
idx = []
with open(path, 'rb') as f:
    pos = 0
    while True:
        h = f.read(14)
        if not h:
            break
        assert len(h) == 14 and h[:2] == b'TK', f'bad record header at offset {pos}'
        rid, nb = struct.unpack('<IQ', h[2:])
        idx.append((rid, pos, nb))
        pos += 14 + nb
        f.seek(pos)
    ids = [r[0] for r in idx]
    assert len(set(ids)) == len(ids), 'duplicate record ids'
    d = hashlib.sha256()
    for rid, p, nb in sorted(idx):
        f.seek(p)
        left = 14 + nb
        while left:
            chunk = f.read(min(left, 1 << 24))
            d.update(chunk); left -= len(chunk)
            if out: out.write(chunk)
if out: out.close()
print(f'{d.hexdigest()}  {len(idx)} records  {path}')
