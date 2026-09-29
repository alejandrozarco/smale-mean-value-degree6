"""Progress of a smale_bb (v1) or smale_bb_v2 run: python3 c/progress.py runs/d6_v2/d6 [ntasks=49152]
(Diagnostic only; the authoritative completion check for v2 is c/check_done.py.)"""
import sys, time, os
pre = sys.argv[1]; nt = int(sys.argv[2]) if len(sys.argv) > 2 else 49152
rows = []
for l in open(pre + ".done"):
    f = l.split()
    if f and f[0] == 'R' and len(f) == 15:        # v2: R id pr F E L out sym excl unres md sec vol cfg chk
        rows.append(f[1:13])
    elif len(f) == 12:                             # v1
        rows.append(f)
ids = {int(r[0]) for r in rows}
pr = sum(int(r[1]) for r in rows); cpu = sum(float(r[10]) for r in rows); un = sum(int(r[8]) for r in rows)
age = time.time() - os.path.getmtime(pre + ".done")
print(f"tasks done {len(ids)}/{nt} ({100*len(ids)/nt:.1f}%), boxes {pr:.3e}, task-cpu {cpu:.0f} s, "
      f"rate {pr/max(cpu,1e-9):.0f} boxes/s/thread, unresolved {un}, last write {age:.0f} s ago")
