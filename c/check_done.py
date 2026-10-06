#!/usr/bin/env python3
"""check_done.py — independent completion checker for a smale_bb_v2 or smale_bb_v3 certification run.

usage: python3 c/check_done.py PREFIX [--src c/smale_bb_vN.c] [--build DIR/BUILD.txt] [--bin DIR/smale_bb_vN.bin]
                                     [--partial]
PREFIX is the run prefix (e.g. runs/d6_v2/d6): reads PREFIX.done, PREFIX.tasks, PREFIX.unresolved.

Independent of the C code (pure Python, exact rational arithmetic where it matters). Verifies:
  1. header line: checksum; config string parsed; config hash = sha256(config string)[:16];
     src sha256 in the config = sha256 of --src (default c/smale_bb_vN.c, N from the header) = BUILD.txt's;
     binary sha256 = BUILD.txt (if --bin/--build given); excl_mode = 2; 0 < r_excl <= 0.05 and the
     local-certificate hand-over r/(1-r) <= 0.0527 holds exactly for the IEEE value of r_excl.
  2. every record: checksum, exact field format, cfg hash = header's, id in range, each id exactly once,
     all ids 0..ntask-1 present (unless --partial), unresolved = 0, tree count processed = 2*leaves - 1,
     vol_done within 1e-9 (relative) of the top-box volume (diagnostic: float sum).
  3. PREFIX.unresolved contains no UNRESOLVED line.
  4. top-level grid, recomputed with exact rationals: the s^D grid of [-1,1]^D (s a power of two: the closed
     boxes tile [-1,1]^D exactly), which contains the closed chart domain |u_j| <= 1.  Every grid box is either
     listed in PREFIX.tasks exactly once, or is EXACTLY discardable: some |u_v| > 1 on the whole closed box, or
     |u_v| < |u_{v+1}| on the whole closed box, or Im u_2 < 0 on the whole closed box.
     The task list equals an independent re-implementation of the C shuffle (xorshift Fisher-Yates, seed).
Exit status 0 and "CERTIFICATE COMPLETE" only if everything passes.
"""
import sys, hashlib, struct, re, os, itertools
from fractions import Fraction

args = sys.argv[1:]
if not args:
    print(__doc__); sys.exit(2)
pre = args[0]
opt = {}
i = 1
partial = False
while i < len(args):
    if args[i] == '--partial':
        partial = True; i += 1
    else:
        opt[args[i]] = args[i + 1]; i += 2
fails = []
def fail(msg):
    fails.append(msg)
    print('FAIL:', msg)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def check_line(line):
    """line without '\\n' -> payload or None"""
    if len(line) < 18 or line[-17] != ' ':
        return None
    payload, chk = line[:-17], line[-16:]
    return payload if sha(payload.encode())[:16] == chk else None

# ---------- 1. header / config ----------
raw = open(pre + '.done', 'rb').read().decode('ascii')
if raw and not raw.endswith('\n'):
    fail('torn trailing line in .done (the run would truncate it on resume)')
lines = raw.split('\n')[:-1] if raw.endswith('\n') else raw.split('\n')[:-1]
if not lines:
    fail('empty .done'); print('RESULT: FAIL'); sys.exit(1)
hdr = check_line(lines[0])
if hdr is None:
    fail('header checksum'); print('RESULT: FAIL'); sys.exit(1)
m = re.fullmatch(r'H (smale_bb_v[23]) cfg=([0-9a-f]{16}) ((smale_bb_v[23]);src=([0-9a-f]{64});d=(\d+);r_excl=([0-9a-f]{16});'
                 r'excl_mode=(\d+);s=(\d+);seed=(\d+))', hdr)
if not m or m.group(1) != m.group(4):
    fail('header format'); print('RESULT: FAIL'); sys.exit(1)
version, cfg16, cfgs, _, src_sha, d, rbits, em, s, seed = m.groups()
src = opt.get('--src', os.path.join(os.path.dirname(os.path.abspath(__file__)), version + '.c'))
print('version:', version)
d, em, s, seed = int(d), int(em), int(s), int(seed)
r = struct.unpack('>d', bytes.fromhex(rbits))[0]
if sha(cfgs.encode())[:16] != cfg16:
    fail('config hash does not match config string')
if sha(open(src, 'rb').read()) != src_sha:
    fail(f'source sha256 of {src} differs from the config (records were produced by another source)')
if '--bin' in opt and '--build' not in opt:
    fail('--bin needs --build (the binary hash is taken from BUILD.txt)')
if '--build' in opt:
    b = open(opt['--build']).read()
    mm = re.search(r'source sha256: ([0-9a-f]{64})', b)
    if not mm or mm.group(1) != src_sha:
        fail('BUILD.txt source sha256 differs from config')
    mb = re.search(r'binary sha256: ([0-9a-f]{64})', b)
    if '--bin' in opt:
        if not os.path.exists(opt['--bin']):
            fail(f"--bin {opt['--bin']} does not exist")
        elif not mb or sha(open(opt['--bin'], 'rb').read()) != mb.group(1):
            fail('binary sha256 differs from BUILD.txt')
if em != 2:
    fail(f'excl_mode {em} != 2 (only the Euclidean u-ball matches local/LOCAL_CERT.md)')
rq = Fraction(r)
if not (0 < rq <= Fraction(1, 20) + Fraction(1, 10**17)) or not (rq / (1 - rq) <= Fraction(527, 10000)):
    fail(f'r_excl = {r!r} not covered by the local certificate (need r/(1-r) <= 0.0527)')
if s < 1 or s & (s - 1) or s > 2**20:
    fail(f's = {s} is not a power of two')
n, nv, D = d - 1, d - 2, 2 * (d - 2)
print(f'config: d={d} r_excl={r!r} (r/(1-r) = {float(rq/(1-rq)):.9f}) excl_mode={em} s={s} seed={seed} cfg={cfg16}')
print(f'source sha256 {src_sha}')

# ---------- 4. top-level grid (exact) ----------
th = Fraction(1, s)
def box_of(idx):
    q, c = idx, [None] * D
    for k in range(D - 1, -1, -1):
        c[k] = -1 + th * (2 * (q % s) + 1); q //= s
    return c
def exact_discard(c):
    h = th
    mn2, mx2 = [], []
    for v in range(nv):
        x, y = c[2 * v], c[2 * v + 1]
        dx, dy = max(abs(x) - h, 0), max(abs(y) - h, 0)
        ex, ey = abs(x) + h, abs(y) + h
        mn2.append(dx * dx + dy * dy); mx2.append(ex * ex + ey * ey)
    if any(m > 1 for m in mn2):
        return 'outside'
    if any(mx2[v] < mn2[v + 1] for v in range(nv - 1)):
        return 'symmetry'
    if c[1] + h < 0:
        return 'Im u2 < 0'
    return None
# tasks file
tl = open(pre + '.tasks').read().split('\n')
mt = re.fullmatch(r'# ' + version + r' tasks cfg=([0-9a-f]{16}) (\S+) ntask=(\d+) half_width=(\S+)', tl[0])
if not mt or mt.group(1) != cfg16 or mt.group(2) != cfgs:
    fail('tasks file header does not match the config')
ntask = int(mt.group(3)) if mt else 0
if mt and Fraction(float.fromhex(mt.group(4))) != th:
    fail('tasks file half-width')
tasks = []
for t, l in enumerate(tl[1:1 + ntask]):
    f = l.split()
    if int(f[0]) != t or len(f) != D + 1:
        fail(f'tasks file line {t}'); break
    tasks.append(tuple(Fraction(float.fromhex(x)) for x in f[1:]))
if len(tl) != ntask + 2 or tl[-1] != '':
    fail('tasks file length')
taskset = {}
for t, c in enumerate(tasks):
    if c in taskset:
        fail(f'task box listed twice: {taskset[c]} and {t}')
    taskset[c] = t
kept, disc = [], {}
for idx in range(s ** D):
    c = tuple(box_of(idx))
    if c in taskset:
        kept.append(c)
    else:
        why = exact_discard(c)
        if why is None:
            fail(f'grid box {idx} is neither a task nor exactly discardable'); break
        disc[why] = disc.get(why, 0) + 1
if len(kept) != len(tasks):
    fail(f'{len(tasks) - len(kept)} task boxes are not grid boxes')
# independent re-implementation of the C shuffle
M64 = (1 << 64) - 1
order = [c for c in (tuple(box_of(i)) for i in range(s ** D)) if c in taskset]
rs = (seed * 2654435761 + 12345) & M64
for i in range(len(order) - 1, 0, -1):
    rs ^= (rs << 13) & M64; rs ^= rs >> 7; rs ^= (rs << 17) & M64
    j = rs % (i + 1)
    order[i], order[j] = order[j], order[i]
if order != tasks:
    fail('task order differs from the re-implemented shuffle')
print(f'grid: {s}^{D} = {s ** D} boxes tile [-1,1]^{D}; {len(tasks)} tasks; exactly discardable: {disc}')

# ---------- 2. records ----------
seen = {}
tot = dict(processed=0, F=0, E=0, L=0, outside=0, sym=0, excl=0, unres=0)
maxdepth, cpu = 0, 0.0
topvol = (2 * th) ** D
rec = re.compile(r'R (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+\.\d{3}) (\S+) cfg=([0-9a-f]{16})')
for ln, line in enumerate(lines[1:], start=2):
    p = check_line(line)
    if p is None:
        fail(f'line {ln}: checksum'); continue
    mr = rec.fullmatch(p)
    if not mr:
        fail(f'line {ln}: format'); continue
    g = mr.groups()
    t = int(g[0]); pr, F, E, L, o, sy, ex, un, md = map(int, g[1:10])
    if g[12] != cfg16:
        fail(f'line {ln}: cfg hash {g[12]} != {cfg16}')
    if not 0 <= t < ntask:
        fail(f'line {ln}: id {t} out of range')
    if t in seen:
        fail(f'line {ln}: duplicate id {t} (first at line {seen[t]})')
    seen[t] = ln
    leaves = F + E + L + o + sy + ex + un
    if pr != 2 * leaves - 1:
        fail(f'line {ln}: tree count processed={pr} != 2*leaves-1={2 * leaves - 1}')
    if un:
        fail(f'line {ln}: task {t} has {un} unresolved boxes')
    if abs(Fraction(float(g[11])) - topvol) > topvol * Fraction(1, 10**9):
        fail(f'line {ln}: task {t} volume {g[11]} != top-box volume')
    for k, v in zip(['processed', 'F', 'E', 'L', 'outside', 'sym', 'excl', 'unres'], [pr, F, E, L, o, sy, ex, un]):
        tot[k] += v
    maxdepth = max(maxdepth, md); cpu += float(g[10])
missing = ntask - len(seen)
if missing and not partial:
    fail(f'{missing} of {ntask} task ids missing')
# ---------- 3. unresolved file ----------
if os.path.exists(pre + '.unresolved'):
    nu = sum(1 for l in open(pre + '.unresolved') if l.startswith('UNRESOLVED'))
    if nu:
        fail(f'{nu} UNRESOLVED boxes in {pre}.unresolved')
print(f'records: {len(seen)}/{ntask} tasks, boxes {tot}, maxdepth {maxdepth}, task-cpu {cpu:.0f} s')
if fails:
    print(f'RESULT: FAIL ({len(fails)} problems)'); sys.exit(1)
if partial and missing:
    print(f'RESULT: OK SO FAR ({missing} tasks still to do; all present records valid)'); sys.exit(0)
print('RESULT: CERTIFICATE COMPLETE (all tasks exactly once, 0 unresolved, config/source bound, grid exact)')
