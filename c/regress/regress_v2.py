"""regress_v2.py — targeted regression of smale_bb_v2: domain refusals, boxes at extreme scales (evaluated with the
domain check lifted in diagnostic modes), equality points, and the checkpoint protocol.

usage: python3 c/regress/regress_v2.py BIN D5_RUN_PREFIX SCRATCH_DIR
  BIN            a smale_bb_v2 binary (built by c/build_v2.sh)
  D5_RUN_PREFIX  a COMPLETE d=5 run of that binary (s=4, seed 1), e.g. runs/d5_v2/d5 (used for the resume tests)
  SCRATCH_DIR    empty directory for temporary run copies
Prints one line per test and a final PASS/FAIL summary."""
import sys, os, subprocess, shutil, json, math
from fractions import Fraction as Q
import mpmath as mp
mp.mp.dps = 60
BIN, D5, SCR = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(SCR, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.join(HERE, '..', 'check_done.py')
results = []
def rep(name, ok, info=''):
    results.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name + (': ' + info if info else ''), flush=True)
def run(args, inp=None, env=None):
    e = dict(os.environ); e.pop('SMALE_V2_DIAG_NO_DOMAIN', None)
    if env: e.update(env)
    return subprocess.run([BIN] + args, input=inp, capture_output=True, text=True, env=e)
DIAG = {'SMALE_V2_DIAG_NO_DOMAIN': '1'}

def S_all(u, n):
    """log|S_i| for u = (u_2..u_n) (mpc), b_1 = 1, b_j = 1/u_j; S_i = int_0^1 prod_{j != i}(1 - t b_i/b_j) dt."""
    b = [mp.mpc(1)] + [1 / x for x in u]
    out = []
    for i in range(n):
        f = lambda t: mp.fprod([1 - t * b[i] / b[j] for j in range(n) if j != i])
        out.append(mp.log(abs(mp.quad(f, [0, 1]))))
    return out

# ---- 1. affine-model enclosure failure at scale 2^-200 (d = 6, models mode) ----
t, h = 2.0 ** -200, 2.0 ** -210
box = ' '.join(repr(x) for x in [t, 0.0] * 4 + [h, 0.0] * 4) + '\n'
r = run(['models', '6'], box)
rep('models_2^-200 refused (domain A5)', r.stdout.strip() == 'REFUSED', r.stdout.strip()[:40])
r = run(['models', '6'], box, DIAG)
rows = [l.split() for l in r.stdout.strip().split('\n')]
valid = [int(x[0]) for x in rows]
# exact check of any model still claimed valid: S_i(s,s,s,s) = 1/5 - 1/(30 s) at s = t + h (the test point)
sv = mp.mpf(t) + mp.mpf(h)
ok = True; info = []
for i, row in enumerate(rows):
    if int(row[0]):
        a, e = float(row[1]), float(row[2]); G = [float(x) for x in row[5:]]
        Si = (mp.mpf(1) / 5 - 1 / (30 * sv)) if i > 0 else None
        if Si is not None:
            pred = a + sum(G[2 * k] * h for k in range(4))
            err = abs(mp.log(abs(Si)) - pred)
            if err > e: ok = False; info.append(f'i={i} err/e={float(err / e):.3g}')
rep('models_2^-200 (diag, no domain): no invalid enclosure', ok and valid[1:] == [0, 0, 0, 0],
    f'valid flags {valid} (v1: 1 for i=2..5 with err/e = 87)' + (' ' + ' '.join(info) if info else ''))

# ---- 2. modrange symmetry discard at s = 2^-537 (d = 6) ----
s = 2.0 ** -537
C = [s * x for x in (11 / 16, 11 / 16, 3 / 4, 0, 1 / 4, 0, 1 / 8, 0)]
H = [2.0 ** -550] * 8
box = ' '.join(x.hex() for x in C + H)
line = ' '.join(repr(x) for x in C + H) + '\n'
r = run(['geom', '6'], line)
rep('modrange_2^-537 refused (domain A5)', r.stdout.strip() == 'REFUSED')
r = run(['geom', '6'], line, DIAG)
rep('modrange_2^-537 (diag): NOT discarded', r.stdout.strip() == '0', f'geom flag {r.stdout.strip()}')
r = run(['box', '6', '0.05', '2'] + [repr(x) for x in C + H], None, DIAG)
import re
j = json.loads(re.sub(r': (-?nan|-?inf)', ': null', r.stdout))  # vol_frac = 0/0 at this scale
rep('modrange_2^-537 (diag, box mode): symmetry = 0', j['symmetry'] == 0, r.stdout.strip()[:120])
r = run(['box', '6', '0.05', '2'] + [repr(x) for x in C + H])
rep('modrange_2^-537 box mode refused', r.returncode == 4 and 'refused' in r.stdout)
# same shape at an in-domain scale must not be discarded either
C2 = [2.0 ** -20 * x for x in (11 / 16 * 16, 11 / 16 * 16, 12, 0, 4, 0, 2, 0)]
H2 = [2.0 ** -30] * 8
r = run(['geom', '6'], ' '.join(repr(x) for x in C2 + H2) + '\n')
rep('same configuration at scale 2^-16 (in domain): not discarded', r.stdout.strip() == '0', r.stdout.strip())

# ---- 3. domain / configuration refusals ----
for args, why in [(['run', '6', '0.05', '3', '1', f'{SCR}/x'], 's=3 (gaps)'),
                  (['run', '6', '0.05', '4', '1', f'{SCR}/x', '1'], 'excl_mode 1'),
                  (['run', '6', '0.05', '4', '1', f'{SCR}/x', '0'], 'excl_mode 0'),
                  (['run', '6', '0.06', '4', '1', f'{SCR}/x'], 'r_excl 0.06')]:
    r = run(args)
    rep(f'run refuses {why}', r.returncode == 5 and 'REFUSED' in r.stderr, r.stderr.strip()[:90])
for c, hh, why in [(0.3, 0.25, 'centre not a multiple of h'), (0.0, 0.3, 'h not a power of two'),
                   (0.0, 2.0 ** -41, 'h < 2^-40'), (1.0, 0.25, '|c|+h > 1')]:
    r = run(['eval', '6'], ' '.join(repr(x) for x in [c] + [0.0] * 7 + [hh] + [0.25] * 7) + '\n')
    rep(f'eval refuses box with {why}', r.stdout.strip() == 'REFUSED')

# ---- 4. equality points (closed forms) vs mpmath ----
worst = 0
for d in (4, 5, 6, 7):
    n = d - 1
    r = run(['roots', str(d)])
    pts = [[float.fromhex(x) for x in l.split()] for l in r.stdout.strip().split('\n')]
    exact = set()
    import itertools
    got = 0
    for p in itertools.permutations(range(1, n)):
        ex = [z for a in p for z in (mp.cos(2 * mp.pi * a / n), -mp.sin(2 * mp.pi * a / n))]
        best = min(max(abs(mp.mpf(q[k]) - ex[k]) for k in range(len(ex))) for q in pts)
        worst = max(worst, float(best))
    rep(f'd={d}: {len(pts)} equality points = all permutations', len(pts) == math.factorial(n - 1))
rep('equality-point coordinate error <= 1e-15 (exclusion margin 1e-13 >= sqrt(8)*1e-15)', worst <= 1e-15, f'max {worst:.3g}')

# ---- 5. checkpoint protocol (on a copy of a complete d=5 run) ----
def fresh(tag):
    dst = os.path.join(SCR, tag); shutil.rmtree(dst, ignore_errors=True); os.makedirs(dst)
    for ext in ('.done', '.tasks', '.unresolved'):
        shutil.copy(D5 + ext, os.path.join(dst, 'd5' + ext))
    return os.path.join(dst, 'd5')
def resume(pre):
    return run(['run', '5', '0.05', '4', '1', pre])
def checker(pre):
    return subprocess.run([sys.executable, CHECK, pre], capture_output=True, text=True)
lines = open(D5 + '.done').read().split('\n')[:-1]
# (a) torn-line scenario: last 5 records lost, torn "1" (no newline) left behind
pre = fresh('torn')
with open(pre + '.done', 'w') as f:
    f.write('\n'.join(lines[:-5]) + '\n' + '1')
lost = [int(l.split()[1]) for l in lines[-5:]]
r = resume(pre)
ck = checker(pre)
new = open(pre + '.done').read().split('\n')[:-1]
redone = sorted(int(l.split()[1]) for l in new[len(lines) - 5:])
rep('torn trailing "1": truncated, the 5 lost tasks redone, checker passes',
    'truncating torn trailing line' in r.stderr and redone == sorted(lost) and ck.returncode == 0 and len(new) == len(lines),
    f'redone {redone}; checker: {ck.stdout.strip().splitlines()[-1]}')
# (b) torn fragment that would have merged into "1234 ..." under v1
pre = fresh('merge')
with open(pre + '.done', 'w') as f:
    f.write('\n'.join(lines[:-1]) + '\n' + 'R 1')
r = resume(pre)
ck = checker(pre)
rep('torn "R 1" + resume: no merged record, checker passes', ck.returncode == 0, ck.stdout.strip().splitlines()[-1])
# (c) corrupted complete line in the middle
pre = fresh('corrupt')
bad = list(lines); f = bad[100].split(); f[2] = str(int(f[2]) + 2); bad[100] = ' '.join(f)
open(pre + '.done', 'w').write('\n'.join(bad) + '\n')
r = resume(pre); ck = checker(pre)
rep('corrupted middle record: run refuses (exit 7), checker fails', r.returncode == 7 and ck.returncode != 0, r.stderr.strip()[:80])
# (d) duplicate record
pre = fresh('dup')
open(pre + '.done', 'w').write('\n'.join(lines + [lines[50]]) + '\n')
r = resume(pre); ck = checker(pre)
rep('duplicate record: run refuses, checker fails', r.returncode == 7 and ck.returncode != 0, r.stderr.strip()[:80])
# (e) other configuration on the same prefix (seed 2 / r_excl 0.04)
pre = fresh('cfg')
r = run(['run', '5', '0.05', '4', '1', pre, '2', '2'])
rep('different seed on same prefix: refused', r.returncode == 6, r.stderr.strip()[:80])
r = run(['run', '5', '0.04', '4', '1', pre])
rep('different r_excl on same prefix: refused', r.returncode == 6, r.stderr.strip()[:80])
# (f) missing task -> checker reports incomplete
pre = fresh('missing')
open(pre + '.done', 'w').write('\n'.join(lines[:-1]) + '\n')
ck = checker(pre)
rep('one record missing: checker FAILS', ck.returncode != 0 and 'missing' in ck.stdout)
# (g) record with unresolved > 0 (valid checksum) -> checker fails
pre = fresh('unres')
import hashlib
f = lines[10].split(); payload = ' '.join(f[:-1]).split(' ')
# fields: R id pr F E L out sym excl unres ...: bump unresolved and processed consistently (tree count stays valid)
payload[9] = str(int(payload[9]) + 1); payload[2] = str(int(payload[2]) + 2)
p = ' '.join(payload); forged = p + ' ' + hashlib.sha256(p.encode()).hexdigest()[:16]
open(pre + '.done', 'w').write('\n'.join(lines[:10] + [forged] + lines[11:]) + '\n')
ck = checker(pre)
rep('record with unresolved = 1: checker FAILS', ck.returncode != 0 and 'unresolved' in ck.stdout)

print('SUMMARY: %d/%d passed' % (sum(ok for _, ok in results), len(results)))
sys.exit(0 if all(ok for _, ok in results) else 1)
