"""cmp_export_counts.py — compare the per-task counts written by c/export_tree.c (OUT.counts) with a run record (.done).
usage: python3 c/regress/cmp_export_counts.py OUT/d6.counts runs/d6_v2/d6.done   (exit 0 iff identical for every task)"""
import os, sys
cnt, done = sys.argv[1], sys.argv[2]
rec = {}
for l in open(done):
    f = l.split()
    if f[0] == 'R':
        rec[int(f[1])] = list(map(int, f[2:11]))
got = {}
for l in open(cnt):
    f = l.split()
    got[int(f[0])] = list(map(int, f[1:]))
diff = [k for k in rec if rec[k] != got.get(k)]
print('run records', len(rec), 'exported tasks', len(got), 'tasks with differing counts', len(diff), diff[:5])
tree = cnt[:-len('.counts')] + '.tree'
if os.path.exists(tree):
    print('tree bytes', os.path.getsize(tree), 'expected', sum(v[0] for v in got.values()) + 14 * len(got))
sys.exit(0 if not diff and len(got) == len(rec) else 1)
