"""replay_tasks.py — recompute a sample of top-box tasks of a completed run with a locally built binary and compare
the per-task counts with the run record.
usage: python3 c/regress/replay_tasks.py BIN RUN_PREFIX [n=20] [seed=1] [id ...]
  BIN         a smale_bb_v2 binary built with c/build_v2.sh
  RUN_PREFIX  e.g. runs/d6_v2/d6 (reads .tasks and .done; d, r_excl, excl_mode from the header)
Runs `BIN box d r_excl excl_mode <centres> <half-widths>` for n random task ids (plus any ids given) and compares
processed/F/E/L/outside/sym/excl/unresolved with the .done record.  Exit status 0 iff all agree."""
import sys, json, random, subprocess, re
BIN, PRE = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 20
R = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
extra = [int(x) for x in sys.argv[5:]]
lines = open(PRE + '.tasks').read().split('\n')
hdr = lines[0]
d = int(re.search(r';d=(\d+);', hdr).group(1))
rbits = re.search(r'r_excl=([0-9a-f]{16})', hdr).group(1)
import struct
r = struct.unpack('>d', bytes.fromhex(rbits))[0]
em = int(re.search(r'excl_mode=(\d+)', hdr).group(1))
hw = re.search(r'half_width=(\S+)', hdr).group(1)
task = {int(l.split()[0]): l.split()[1:] for l in lines[1:] if l.strip()}
rec = {}
for l in open(PRE + '.done'):
    f = l.split()
    if f[0] == 'R':
        rec[int(f[1])] = list(map(int, f[2:10]))
ids = R.sample(sorted(rec), min(n, len(rec))) + extra
KEYS = ['processed', 'F', 'E', 'L', 'outside', 'symmetry', 'excluded', 'unresolved']
ok = True
for t in ids:
    c = task[t]
    out = subprocess.run([BIN, 'box', str(d), repr(r), str(em)] + c + [hw] * len(c), capture_output=True, text=True)
    res = json.loads(out.stdout)
    got = [res[k] for k in KEYS]
    same = got == rec[t]
    ok &= same
    print(f'task {t}: processed {got[0]}  {"same as record" if same else "DIFFERENT: record " + str(rec[t]) + " replay " + str(got)}', flush=True)
print('ALL SAME' if ok else 'MISMATCH'); sys.exit(0 if ok else 1)
