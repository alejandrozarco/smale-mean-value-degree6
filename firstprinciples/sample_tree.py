#!/usr/bin/env python3
"""sample_tree.py -- end-to-end falsification test of the subdivision trees
against the true critical-value quotients V_i (computed in fpcore.py from P
built from its critical points).

usage:  python3 sample_tree.py D_DEGREE STEP [options]
steps:  prep stats geom select eval intense report   (or 'all')

Large intermediate files go to $FP_WORK (default: ./work), small results and
logs are printed (redirect to a log file).
"""
import os, sys, struct, subprocess, time, math, itertools, json
from fractions import Fraction
import numpy as np
import mpmath as mp
from fpcore import crit_values, crit_values_mp, x_to_u, EPS

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TASKS = {5: os.environ.get('FP_TASKS_D5', os.path.join(REPO, 'runs/d5_v2/d5.tasks')), 6: os.environ.get('FP_TASKS_D6', os.path.join(REPO, 'runs/d6_v2/d6.tasks'))}
# tree files: the d = 6 tree is the release asset d6.canonical.tree.zst (decompressed); the d = 5 tree is written by
# c/export_tree.c (README, Reproduce). Override with FP_TREE_D5 / FP_TREE_D6.
TREES = {5: os.environ.get('FP_TREE_D5', os.path.join(REPO, 'd5.tree')),
         6: os.environ.get('FP_TREE_D6', os.path.join(REPO, 'd6.canonical.tree'))}
TYPES = ['F', 'E', 'L', 'outside', 'symmetry', 'excluded', 'unresolved']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


class Cfg:
    def __init__(self, d):
        self.d = d; self.n = d - 1; self.D = 2 * (d - 2)
        self.c = Fraction(d - 1, d); self.cf = (d - 1) / d
        self.work = os.path.join(os.environ.get('FP_WORK', os.path.join(HERE, 'work')), 'd%d' % d)
        os.makedirs(self.work, exist_ok=True)
        self.pre = os.path.join(self.work, 'd%d' % d)


# ---------------------------------------------------------------- prep
def read_tasks(cfg):
    with open(TASKS[cfg.d]) as f:
        hdr = f.readline()
        h = float.fromhex(hdr.split('half_width=')[1].split()[0])
        rows = []
        for line in f:
            p = line.split()
            if not p: continue
            assert len(p) == 1 + cfg.D
            rows.append((int(p[0]), [float.fromhex(x) for x in p[1:]]))
    rows.sort()
    assert [r[0] for r in rows] == list(range(len(rows)))
    return h, np.array([r[1] for r in rows]), hdr


def equality_points(cfg, dps=40):
    """u-coordinates of the equality points p = (conj w^{a_2},...,conj w^{a_n})."""
    n = cfg.n
    pts_mp, perms = [], []
    with mp.workdps(dps):
        for a in itertools.permutations(range(1, n)):
            x = []
            for aj in a:
                z = mp.exp(-2j * mp.pi * aj / n)
                x += [mp.re(z), mp.im(z)]
            pts_mp.append(x); perms.append(a)
    return perms, pts_mp


def step_prep(cfg):
    h, C, hdr = read_tasks(cfg)
    log('tasks', len(C), 'half_width', h, '| header:', hdr.strip()[:200])
    with open(cfg.pre + '.taskbin', 'wb') as f:
        f.write(struct.pack('<iid', len(C), cfg.D, h))
        for t in range(len(C)):
            f.write(struct.pack('<i', t)); f.write(C[t].astype('<f8').tobytes())
    perms, pts = equality_points(cfg)
    with open(cfg.pre + '.eqbin', 'wb') as f:
        f.write(struct.pack('<i', len(pts)))
        for x in pts:
            f.write(np.array([float(t) for t in x], '<f8').tobytes())
    log('equality points', len(pts))
    # coverage of the region by the top boxes: grid of side 2h over [-1,1]^D
    k = int(round(1 / h))
    centers = [-1 + h * (2 * i + 1) for i in range(k)]
    have = set(map(tuple, C))
    nmiss, bad = 0, 0
    for c in itertools.product(centers, repeat=cfg.D):
        if c in have: continue
        nmiss += 1
        # a missing cell must contain no point of R: here Im u_2 <= c+h < 0 or |u_v|>1 on all of cell
        if not (c[1] + h < 0):
            bad += 1
    log('top-box grid: %d cells present, %d missing; missing cells not certified disjoint from R: %d'
        % (len(have), nmiss, bad))
    log('R is contained in [-1,1]^D since |u_j|<=1; grid covers [-1,1]^D:', abs(centers[0] - h + 1) < 1e-15 and abs(centers[-1] + h - 1) < 1e-15)


# ---------------------------------------------------------------- stats (C walker)
def step_stats(cfg):
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '19', os.path.join(HERE, 'walk'), 'stats', TREES[cfg.d],
                        cfg.pre + '.taskbin', cfg.pre + '.eqbin', cfg.pre])
    log('walker stats exit code', r.returncode, 'time %.1fs' % (time.time() - t0))
    print(open(cfg.pre + '.stats.summary.txt').read().split('depth_hist')[0])
    tk = np.loadtxt(cfg.pre + '.stats.tasks.txt', dtype=np.int64, ndmin=2)
    log('records', len(tk), 'distinct ids', len(set(tk[:, 0])), 'all complete', bool((tk[:, 2] == 1).all()),
        'extra bytes', int(tk[:, 3].sum()), 'bad codes', int(tk[:, 4].sum()), 'unresolved leaves', int(tk[:, 12].sum()))
    s = os.path.getsize(TREES[cfg.d])
    log('file size', s, '= sum(14 + nbytes)?', s == int((14 + tk[:, 1]).sum()))
    nl = tk[:, 6:13].sum(1)
    nint = tk[:, 1] - nl
    log('every record: #internal = #leaves - 1 ?', bool((nint == nl - 1).all()))


# ---------------------------------------------------------------- exact geometric recheck
def read_geom(cfg, path):
    D = cfg.D
    rs = 4 + 4 + 8 + 16 * D
    out = []
    with open(path, 'rb') as f:
        while True:
            b = f.read(rs)
            if len(b) < rs: break
            tid, code = struct.unpack('<II', b[:8]); li, = struct.unpack('<q', b[8:16])
            cen = np.frombuffer(b[16:16 + 8 * D], '<f8'); half = np.frombuffer(b[16 + 8 * D:], '<f8')
            out.append((tid, code, li, cen.copy(), half.copy()))
    return out


def exact_geom_ok(cfg, code, cen, half, eq_mp):
    lo = [Fraction(float(c)) - Fraction(float(h)) for c, h in zip(cen, half)]
    hi = [Fraction(float(c)) + Fraction(float(h)) for c, h in zip(cen, half)]
    def mn(a1, a2):
        return Fraction(0) if a1 <= 0 <= a2 else min(abs(a1), abs(a2))
    def mx(a1, a2):
        return max(abs(a1), abs(a2))
    m = cfg.n - 1
    if code == 3:
        vals = [mn(lo[2*v], hi[2*v]) ** 2 + mn(lo[2*v+1], hi[2*v+1]) ** 2 - 1 for v in range(m)]
        return max(vals) > 0, float(max(vals))
    if code == 4:
        vals = [-hi[1]]
        for v in range(m - 1):
            vals.append(mn(lo[2*v+2], hi[2*v+2]) ** 2 + mn(lo[2*v+3], hi[2*v+3]) ** 2
                        - mx(lo[2*v], hi[2*v]) ** 2 - mx(lo[2*v+1], hi[2*v+1]) ** 2)
        return max(vals) > 0, float(max(vals))
    if code == 5:
        with mp.workdps(60):
            best = -mp.inf
            for p in eq_mp:
                s = mp.mpf(0)
                for k in range(cfg.D):
                    l = mp.mpf(lo[k].numerator) / lo[k].denominator
                    u = mp.mpf(hi[k].numerator) / hi[k].denominator
                    t = max(abs(l - p[k]), abs(u - p[k])); s += t * t
                best = max(best, mp.mpf(1) / 400 - s)
            return best >= 0, float(best)


def step_geom(cfg):
    print(open(cfg.pre + '.stats.summary.txt').read().split('depth_hist')[0])
    g = read_geom(cfg, cfg.pre + '.stats.geom_flag.bin')
    log('geometric leaves flagged by the double check (borderline or failing):', len(g))
    _, eq = equality_points(cfg, 60)
    nbad = 0
    for tid, code, li, cen, half in g:
        ok, marg = exact_geom_ok(cfg, code, cen, half, eq)
        if not ok:
            nbad += 1
            if nbad <= 20:
                log('  EXACT VIOLATION', TYPES[code], 'task', tid, 'leaf', li, 'margin', marg)
        else:
            log('  exact recheck ok', TYPES[code], 'task', tid, 'leaf', li, 'exact margin %.3e' % marg)
    log('exact geometric violations:', nbad)


# ---------------------------------------------------------------- select
def step_select(cfg, target):
    tk = np.loadtxt(cfg.pre + '.stats.tasks.txt', dtype=np.int64, ndmin=2)
    ids = tk[:, 0]; nFEL = tk[:, 6:9].sum(1)
    total = int(nFEL.sum())
    p = min(1.0, target / max(total, 1))
    prob = np.zeros(int(ids.max()) + 1)
    kmin = 8
    with np.errstate(divide='ignore'):
        pt = np.where(nFEL > 0, np.maximum(p, np.minimum(1.0, kmin / np.maximum(nFEL, 1))), 0.0)
    prob[ids] = pt
    prob.astype('<f8').tofile(cfg.pre + '.prob')
    # deep threshold: smallest depth T such that #F/E/L leaves with depth >= T is <= 20000
    txt = open(cfg.pre + '.stats.summary.txt').read().split('depth_hist')[1].split('\n')[1:]
    hist = {}
    for line in txt:
        if ':' not in line: continue
        dpt, rest = line.split(':'); v = list(map(int, rest.split()))
        hist[int(dpt)] = v[0] + v[1] + v[2]
    T = max(hist) + 1; acc = 0
    for dpt in sorted(hist, reverse=True):
        if acc + hist[dpt] > 20000: break
        acc += hist[dpt]; T = dpt
    if os.environ.get('FP_DEEP'):
        T = int(os.environ['FP_DEEP']); acc = sum(v for k2, v in hist.items() if k2 >= T)
    log('F/E/L leaves total', total, 'base probability %.3g' % p, 'expected sample %.0f' % (pt * nFEL).sum(),
        'deep threshold', T, '(%d leaves at depth >= T)' % acc, 'max depth', max(hist))
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '19', os.path.join(HERE, 'walk'), 'sample', TREES[cfg.d],
                        cfg.pre + '.taskbin', cfg.pre + '.eqbin', cfg.pre, cfg.pre + '.prob', str(T), os.environ.get('FP_SEED', '12345')])
    log('walker sample exit', r.returncode, 'time %.1fs' % (time.time() - t0))
    print(open(cfg.pre + '.sample.summary.txt').read())


# ---------------------------------------------------------------- eval
def leaf_dtype(D):
    return np.dtype([('id', '<u4'), ('code', '<u4'), ('leaf', '<i8'), ('depth', '<u4'), ('deep', '<u4'),
                     ('cen', '<f8', (D,)), ('half', '<f8', (D,))])


def corner_signs(D):
    return np.array(list(itertools.product([-1.0, 1.0], repeat=D)))


def leaf_points(cen, half, rng, nrand, corners, S):
    """Return (B, m, D) sample points for a batch of leaves."""
    B, D = cen.shape
    pts = [cen[:, None, :]]
    if corners:
        pts.append(cen[:, None, :] + S[None, :, :] * half[:, None, :])
    if nrand:
        r = rng.uniform(-1, 1, size=(B, nrand, D))
        pts.append(cen[:, None, :] + r * half[:, None, :])
    return np.concatenate(pts, axis=1)


GAMMA = 64 * EPS   # generous constant for the floating-point error estimate


def margins(cfg, code, absV, errV, valid):
    """absV, errV: (B, m, n); valid: (B, m).  Returns per-leaf
    (margin, pointwise_margin, aux, errest, maxmin) with margin > 0 meaning the
    necessary condition of the claim holds at all valid sample points."""
    n = cfg.n; c = cfg.cf
    B, m, _ = absV.shape
    big = 1e300
    minV = absV.min(2)
    maxmin = np.where(valid, minV, -big).max(1)
    if code == 0:   # F
        per_i = np.where(valid[:, :, None], c - absV, big).min(1)       # (B, n)
        aux = per_i.argmax(1); marg = per_i.max(1)
        pw = np.where(valid, c - minV, big).min(1)
        err = np.where(valid[:, :, None], errV, 0).max((1, 2))
        return marg, pw, aux, err, maxmin
    logV = np.log(absV)
    relerr = errV / absV
    if code == 2:   # L
        s = logV.sum(2)
        marg = np.where(valid, n * math.log(c) - s, big).min(1)
        err = np.where(valid, relerr.sum(2), 0).max(1)
        return marg, marg, np.zeros(B, int), err, maxmin
    # E
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    best = np.full(B, -big); aux = np.zeros(B, int)
    for pi, (i, j) in enumerate(pairs):
        dl = logV[:, :, i] - logV[:, :, j]
        lo = np.where(valid, dl, big).min(1); hi = np.where(valid, dl, -big).max(1)
        mpair = np.maximum(lo, -hi)          # >0 iff constant strict sign at all points
        better = mpair > best
        best = np.where(better, mpair, best); aux = np.where(better, pi, aux)
    err = np.where(valid, relerr.max(2) * 2, 0).max(1)
    return best, best, aux, err, maxmin


def mp_leaf_margin(cfg, code, pts, aux_i=None):
    """Recompute the leaf's sampled margin in 50-digit arithmetic."""
    n = cfg.n; c = mp.mpf(cfg.d - 1) / cfg.d
    rows = []
    with mp.workdps(50):
        for x in pts:
            u = [mp.mpc(mp.mpf(float(x[2*k])), mp.mpf(float(x[2*k+1]))) for k in range(n - 1)]
            if any(z == 0 for z in u): continue
            rows.append([abs(V) for V in crit_values_mp(u)])
        if not rows: return None
        if code == 0:
            per_i = [min(c - r[i] for r in rows) for i in range(n)]
            return max(per_i), min(c - min(r) for r in rows)
        if code == 2:
            return min(n * mp.log(c) - mp.fsum(mp.log(v) for v in r) for r in rows), None
        best = -mp.inf
        for i in range(n):
            for j in range(i + 1, n):
                dl = [mp.log(r[i]) - mp.log(r[j]) for r in rows]
                best = max(best, max(min(dl), -max(dl)))
        return best, None


def step_eval(cfg, nrand, corner_every, batch_pts=600000, seed=7):
    D, n = cfg.D, cfg.n
    dt = leaf_dtype(D)
    path = cfg.pre + '.sample.leaves.bin'
    nrec = os.path.getsize(path) // dt.itemsize
    log('sampled leaves', nrec, 'nrand', nrand, 'corners for every %d-th leaf and all deep leaves' % corner_every)
    S = corner_signs(D)
    rng = np.random.default_rng(seed)
    res = np.zeros(nrec, dtype=[('code', 'u1'), ('marg', 'f8'), ('pw', 'f8'), ('aux', 'i2'), ('err', 'f8'),
                                ('maxmin', 'f8'), ('npts', 'i4'), ('ndeg', 'i4'), ('mp', 'u1')])
    counts = np.zeros(3, np.int64); viol = []
    rechecked = 0
    t0 = time.time()
    mm = np.memmap(path, dtype=dt, mode='r')
    idx_all = np.arange(nrec)
    use_corners = (idx_all % corner_every == 0) | (np.asarray(mm['deep']) == 1)
    for code in (0, 1, 2):
        for withc in (False, True):
            sel = np.nonzero((np.asarray(mm['code']) == code) & (use_corners == withc))[0]
            m = 1 + nrand + (len(S) if withc else 0)
            B = max(1, batch_pts // m)
            for s0 in range(0, len(sel), B):
                ii = sel[s0:s0 + B]
                rec = mm[ii]
                cen = np.asarray(rec['cen']); half = np.asarray(rec['half'])
                P = leaf_points(cen, half, rng, nrand, withc, S)          # (b, m, D)
                U = x_to_u(P).reshape(-1, n - 1)
                valid = (U != 0).all(1)
                Uv = np.where(valid[:, None], U, 1.0)
                V, A = crit_values(Uv)
                absV = np.abs(V).reshape(len(ii), m, n)
                errV = (GAMMA * A).reshape(len(ii), m, n)
                valid = valid.reshape(len(ii), m)
                marg, pw, aux, err, maxmin = margins(cfg, code, absV, errV, valid)
                res['code'][ii] = code; res['marg'][ii] = marg; res['pw'][ii] = pw; res['aux'][ii] = aux
                res['err'][ii] = err; res['maxmin'][ii] = maxmin; res['npts'][ii] = valid.sum(1)
                res['ndeg'][ii] = (~valid).sum(1)
                counts[code] += len(ii)
                # recheck in mpmath when within 1e-9 (+ error estimate) of a threshold, or failing
                near = np.nonzero((marg <= 1e-9 + err) | ((code == 0) & (pw <= 1e-9 + err)))[0]
                for k in near:
                    rechecked += 1
                    r = mp_leaf_margin(cfg, code, P[k])
                    res['mp'][ii[k]] = 1
                    if r is None: continue
                    mmarg = r[0]
                    if not (mmarg > 0) or (code == 0 and not (r[1] > 0)):
                        viol.append((int(rec['id'][k]), int(rec['leaf'][k]), TYPES[code], float(mmarg)))
                        log('  VIOLATION (mp-confirmed)', viol[-1])
                    else:
                        res['marg'][ii[k]] = float(mmarg)
            log('  type %s corners=%s: %d leaves done (%.0fs)' % (TYPES[code], withc, len(sel), time.time() - t0))
    np.save(cfg.pre + '.evalres.npy', res)
    log('evaluated per type: F %d E %d L %d; mp rechecks %d; violations %d' % (counts[0], counts[1], counts[2], rechecked, len(viol)))
    summarize(cfg, mm, res)


def summarize(cfg, mm, res):
    n = cfg.n
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ids = np.asarray(mm['id'])
    log('distinct task ids among sampled leaves:', len(np.unique(ids)))
    for code in (0, 1, 2):
        s = np.nonzero(res['code'] == code)[0]
        s = s[np.asarray(mm['code'])[s] == code]
        if len(s) == 0: continue
        r = res[s]
        nfail = int((r['marg'] <= 0).sum())
        o = np.argsort(r['marg'])[:5]
        log('%s: %d leaves, %d sampled points, %d degenerate points skipped, failing %d, min margin %.3e'
            % (TYPES[code], len(s), int(r['npts'].sum()), int(r['ndeg'].sum()), nfail, r['marg'].min()))
        for k in o:
            q = s[k]
            extra = ''
            if code == 0: extra = 'best i=%d (V_%d), pointwise margin %.3e' % (r['aux'][k], r['aux'][k] + 1, r['pw'][k])
            if code == 1: extra = 'pair V_%d,V_%d' % (pairs[r['aux'][k]][0] + 1, pairs[r['aux'][k]][1] + 1)
            log('   task %d leaf %d depth %d margin %.3e %s' % (mm['id'][q], mm['leaf'][q], mm['depth'][q], r['marg'][k], extra))
        if code == 0:
            log('   F: leaves where no single i works at all samples: %d; pointwise min_i|V_i|<c fails: %d'
                % (int((r['marg'] <= 0).sum()), int((r['pw'] <= 0).sum())))
    deep = np.asarray(mm['deep']) == 1
    if deep.any():
        log('deep leaves (all leaves at depth >= threshold): %d, depths %d..%d, failing %d, min margin %.3e'
            % (int(deep.sum()), int(np.asarray(mm['depth'])[deep].min()), int(np.asarray(mm['depth'])[deep].max()),
               int((res['marg'][deep] <= 0).sum()), res['marg'][deep].min()))
    k = int(np.argmax(res['maxmin']))
    log('largest min_i|V_i| at any sampled point: %.15f (c = %.15f, c - value = %.3e), task %d leaf %d type %s'
        % (res['maxmin'][k], cfg.cf, cfg.cf - res['maxmin'][k], mm['id'][k], mm['leaf'][k], TYPES[res['code'][k]]))


# ---------------------------------------------------------------- intense test of the closest leaves
def make_obj(cfg, lo, hi):
    n = cfg.n
    def V_of(X):
        U = x_to_u(np.atleast_2d(X))
        return crit_values(U)[0]
    return V_of


def maximize(fun_batch, lo, hi, starts, iters=80):
    """Maximise a smooth scalar function over a box with L-BFGS-B (finite-difference
    gradient evaluated as one batch).  fun_batch: (k, D) -> (k,)."""
    from scipy.optimize import minimize
    D = len(lo); best = -np.inf; bx = None
    scale = (hi - lo)
    for x0 in starts:
        def f(y):
            x = lo + y * scale
            # finite-difference step of relative size 1e-7, taken backwards at the upper face
            steps = np.where(y + 1e-7 <= 1, 1e-7, -1e-7) * scale
            Xs = np.vstack([x] + [x + steps[k] * e for k, e in enumerate(np.eye(D))])
            Xs = np.clip(Xs, lo, hi)
            vals = fun_batch(Xs)
            dy = (Xs[1:] - x)[np.arange(D), np.arange(D)] / scale
            g = np.where(dy != 0, (vals[1:] - vals[0]) / np.where(dy != 0, dy, 1), 0.0)
            if not np.all(np.isfinite(g)): g = np.nan_to_num(g, nan=0.0, posinf=1e30, neginf=-1e30)
            return -vals[0], -g
        y0 = (np.asarray(x0) - lo) / scale
        try:
            r = minimize(f, np.clip(y0, 0, 1), jac=True, method='L-BFGS-B', bounds=[(0, 1)] * D,
                         options={'maxiter': iters})
            val = -r.fun; x = lo + r.x * scale
        except Exception:
            continue
        if val > best: best, bx = val, x
    return best, bx


def step_intense(cfg, K, seed=11):
    """For the K leaves of each type with the smallest sampled margin (plus the 50
    deepest), search adversarially inside the leaf box for a point breaking the
    necessary condition of the claim."""
    n, D, c = cfg.n, cfg.D, cfg.cf
    dt = leaf_dtype(D)
    mm = np.memmap(cfg.pre + '.sample.leaves.bin', dtype=dt, mode='r')
    res = np.load(cfg.pre + '.evalres.npy')
    rng = np.random.default_rng(seed)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    S = corner_signs(D)
    chosen = []
    for code in (0, 1, 2):
        s = np.nonzero(res['code'] == code)[0]
        s = s[np.argsort(res['marg'][s])[:K]]
        chosen += [(int(q), 'closest') for q in s]
    dep = np.asarray(mm['depth'])
    for q in np.argsort(-dep)[:50]:
        chosen.append((int(q), 'deepest'))
    log('intense search on %d leaves' % len(chosen))
    worst = {0: (np.inf, None), 1: (np.inf, None), 2: (np.inf, None)}
    viol = 0
    t0 = time.time()
    for q, why in chosen:
        rec = mm[q]; code = int(rec['code'])
        cen = np.asarray(rec['cen']); half = np.asarray(rec['half'])
        lo, hi = cen - half, cen + half
        # dense random sample + corners
        Xs = np.vstack([cen, cen + S * half, cen + rng.uniform(-1, 1, (4000, D)) * half])
        U = x_to_u(Xs); ok = (U != 0).all(1); Xs = Xs[ok]; U = U[ok]
        absV = np.abs(crit_values(U)[0])
        def fb(X, f):
            Uq = x_to_u(X); okq = (Uq != 0).all(1)
            out = np.full(len(X), -1e300)
            if okq.any():
                out[okq] = f(np.abs(crit_values(Uq[okq])[0]))
            return out
        if code == 0:
            # claim: exists i with |V_i| < c on B.  For each candidate i, maximise |V_i|.
            res_i = []
            for i in range(n):
                v = absV[:, i]
                if v.max() >= c: res_i.append(v.max()); continue
                starts = np.vstack([Xs[np.argsort(-v)[:4]], Xs[-2:]])
                best, _ = maximize(lambda X: fb(X, lambda a: a[:, i]), lo, hi, starts)
                res_i.append(max(best, v.max()))
            margin = c - min(res_i)
        elif code == 2:
            v = np.log(absV).sum(1)
            starts = np.vstack([Xs[np.argsort(-v)[:4]], Xs[-2:]])
            best, _ = maximize(lambda X: fb(X, lambda a: np.log(a).sum(1)), lo, hi, starts)
            margin = n * math.log(c) - max(best, v.max())
        else:
            lv = np.log(absV); margin = -np.inf
            for (i, j) in pairs:
                dl = lv[:, i] - lv[:, j]
                if dl.min() > 0: sgn = 1
                elif dl.max() < 0: sgn = -1
                else: continue
                # try to push sgn*dl below zero
                starts = np.vstack([Xs[np.argsort(sgn * dl)[:3]], Xs[-2:]])
                best, _ = maximize(lambda X: fb(X, lambda a: -sgn * (np.log(a[:, i]) - np.log(a[:, j]))), lo, hi, starts)
                mpair = min(-best, (sgn * dl).min())
                margin = max(margin, mpair)
        if margin <= 1e-9:
            # recheck the dense set in mp only if float says (almost) failing
            log('  candidate failure', TYPES[code], 'task', int(rec['id']), 'leaf', int(rec['leaf']), 'float margin %.3e' % margin)
            if margin <= 0: viol += 1
        if margin < worst[code][0]:
            worst[code] = (margin, (int(rec['id']), int(rec['leaf']), int(rec['depth']), why))
    for code in (0, 1, 2):
        log('intense %s: smallest margin after adversarial search %.3e at (task, leaf, depth, why) %s'
            % (TYPES[code], worst[code][0], worst[code][1]))
    log('intense candidate violations (float):', viol, '(%.0fs)' % (time.time() - t0))


def step_maxmin(cfg, K=40, seed=3):
    """Unconstrained (over R^D, kept inside |u_j|<=1 by clipping) local maximisation of
    min_i |V_i| started from the sampled points with the largest min_i |V_i|."""
    from scipy.optimize import minimize
    n, D, c = cfg.n, cfg.D, cfg.cf
    dt = leaf_dtype(D)
    mm = np.memmap(cfg.pre + '.sample.leaves.bin', dtype=dt, mode='r')
    res = np.load(cfg.pre + '.evalres.npy')
    top = np.argsort(-res['maxmin'])[:K]
    _, eq = equality_points(cfg, 20)
    eqf = np.array([[float(t) for t in p] for p in eq])
    rng = np.random.default_rng(seed)
    best_all = -1
    for q in top:
        cen = np.asarray(mm[q]['cen']); half = np.asarray(mm[q]['half'])
        X = cen + rng.uniform(-1, 1, (2000, D)) * half
        U = x_to_u(X); ok = (U != 0).all(1); X = X[ok]
        v = np.abs(crit_values(x_to_u(X))[0]).min(1)
        x0 = X[np.argmax(v)]
        f = lambda x: -np.abs(crit_values(x_to_u(np.atleast_2d(x)))[0]).min()
        r = minimize(f, x0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 20000, 'maxfev': 40000})
        dist = np.sqrt(((eqf - r.x) ** 2).sum(1)).min()
        best_all = max(best_all, -r.fun)
        log('  start task %d leaf %d: local max of min_i|V_i| = %.15f (c - value %.2e), dist to nearest equality point %.2e, max|u_j| %.4f'
            % (mm[q]['id'], mm[q]['leaf'], -r.fun, c + r.fun, dist, np.abs(x_to_u(r.x)).max()))
    log('largest local max of min_i|V_i| found: %.15f, c = %.15f' % (best_all, c))


if __name__ == '__main__':
    d = int(sys.argv[1]); step = sys.argv[2]
    cfg = Cfg(d)
    steps = ['prep', 'stats', 'geom', 'select', 'eval', 'intense', 'maxmin'] if step == 'all' else [step]
    target = float(os.environ.get('FP_TARGET', '1.2e6'))
    nrand = int(os.environ.get('FP_NRAND', '12'))
    cev = int(os.environ.get('FP_CORNER_EVERY', '1'))
    K = int(os.environ.get('FP_K', '200'))
    for s in steps:
        log('=== step', s, 'd =', d)
        if s == 'prep': step_prep(cfg)
        elif s == 'stats': step_stats(cfg)
        elif s == 'geom': step_geom(cfg)
        elif s == 'select': step_select(cfg, target)
        elif s == 'eval': step_eval(cfg, nrand, cev)
        elif s == 'intense': step_intense(cfg, K)
        elif s == 'maxmin': step_maxmin(cfg)
