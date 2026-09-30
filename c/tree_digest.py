"""tree_digest.py — canonical sha256 of a subdivision tree file written by c/export_tree.c.
The exporter writes records in thread-completion order, so the raw file hash is not reproducible; this digest is the
sha256 of the records ("TK" | uint32 id | uint64 nbytes | codes) concatenated in increasing id order.
usage: python3 c/tree_digest.py FILE.tree [--write CANONICAL.tree]
  --write also writes the records in increasing id order (the canonical file; its sha256 equals the digest)."""
import hashlib, struct, sys

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
