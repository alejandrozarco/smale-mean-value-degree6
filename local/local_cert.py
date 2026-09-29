"""Rigorous local certificate at the equality configuration, degree d (default 6).  See local/LOCAL_CERT.md.

Chart: b_1 = 1, b_j = omega^{j-1} exp(eps_j) (j = 2..n, n = d-1, omega = e^{2 pi i/n}), eps in C^{n-1},
real coordinates x = (Re eps_2, Im eps_2, Re eps_3, ...), t = ||x||_2.
l_i = log|S_i| = Re h_i, h_i = log S_i - log c (holomorphic near 0), L = (1/n) sum_i l_i, G_i = l_i - L,
E = {|S_1| = ... = |S_n|} = {G = 0};  on E, log F = L.

Exact facts (machine-checked in Q(omega) by local/exact_facts.py, asserted below): S_i(0) = c = (d-1)/d, so
h_i(0) = 0 and L(0) = log c; sum_i S_i'(0) = 0, so grad L(0) = 0 and the linear part of G_i is J_i := D l_i(0).

Taylor models (python-flint acb, 128-bit) on the polydisk ||eps||_inf <= T: h_i = P_i + R_i with P_i the
exact degree-<=N Taylor polynomial (ball coefficients) and |R_i| <= r_i on the polydisk; R_i vanishes to order
N+1, so |R_i(eps)| <= r_i (||eps||_inf/T)^{N+1} <= r_i (t/T)^2 for t <= T (Schwarz; the code uses the weaker (t/T)^2).
T is the exact decimal (arb ball containing it), not the nearest double.

Certificate (tube + augmented Hessian), for x in E, 0 < t <= T:
  tube:   |Jx| <= kappa t^2, from Jx = -G^{>=2}(x), x = Pi x + w, |w| <= beta |Jx| + gamma t, and
          y <= q_Pi t^2 + 2 b2 |Pi| t (beta y + gamma t) + b2 (beta y + gamma t)^2 + g3(t)  (fixed point, y = |Jx|);
          Pi = V V^T computed in arb from the stored float basis V (exact product), beta = |J^+|, gamma = |I - Pi - J^+ J|.
  quad:   1/2 x^T H_L x = 1/2 x^T (H_L - sigma J^T J) x + sigma/2 |Jx|^2 <= -mu/2 t^2 + sigma/2 kappa^2 t^4,
          mu certified by arb LDL^T of -(H_L - sigma J^T J) - mu I.
  cubic:  L^3(Pi x) bounded by a rigorous sphere B&B on K (exact rational sphere-intersection tests),
          |L^3(a + w) - L^3(a)| <= C3 ((|a| + |w|)^3 - |a|^3).
  higher: sum_{k=4}^N (sum_a |c_a| m_a) t^k with m_a = max_{|rho|_2 = 1} rho^a = prod_j (a_j/k)^{a_j/2} (arb),
          remainder (1/n) sum_i r_i (t/T)^2.
  => L(x) - log c <= t^2 g(t), g increasing in t;  g(T) < 0 certifies F < c on E within |x|_2 <= T.

usage: python3 local_cert.py [d=6] [T=0.0527 (decimal string)] [N=6] [sigma=5]
"""
import sys, itertools, math, json
from fractions import Fraction
import numpy as np
from flint import acb, arb, ctx
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact_facts


def ipow(x, k):
    """x**k by repeated multiplication (k >= 0 integer).  python-flint 0.6 returns nan for arb ** k when the ball
    contains 0; products are exact ball operations and never produce nan."""
    r = x * 0 + 1 if not isinstance(x, (int, float, Fraction)) else 1
    for _ in range(int(k)):
        r = r * x
    return r


def sq(x):
    return x * x

ctx.prec = 128


# ---------------- Taylor model arithmetic (polydisk radius T, degree <= N) ----------------
class TM:
    __slots__ = ("p", "r")

    def __init__(self, p, r):
        self.p = p          # dict: exponent tuple -> acb
        self.r = r          # arb >= 0: remainder bound on the polydisk


def tm_const(c, nv):
    return TM({(0,) * nv: acb(c)}, arb(0))


def tm_bound_poly(p, T):
    s = arb(0)
    for k, c in p.items():
        s += abs(c) * ipow(T, sum(k))
    return s


def tm_bound(f, T):
    return tm_bound_poly(f.p, T) + f.r


def tm_add(f, g, a=1, b=1):
    p = dict()
    for k, c in f.p.items():
        p[k] = c * a
    for k, c in g.p.items():
        p[k] = p[k] + c * b if k in p else c * b
    return TM(p, f.r * abs(arb(a)) + g.r * abs(arb(b)))


def tm_scale(f, a):
    return TM({k: c * a for k, c in f.p.items()}, f.r * abs(acb(a)))


def tm_mul(f, g, N, T):
    p = {}
    trunc = arb(0)
    for ka, va in f.p.items():
        da = sum(ka)
        for kb, vb in g.p.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            v = va * vb
            if da + sum(kb) <= N:
                p[k] = p[k] + v if k in p else v
            else:
                trunc += abs(v) * ipow(T, da + sum(kb))
    Bf, Bg = tm_bound_poly(f.p, T), tm_bound_poly(g.p, T)
    r = f.r * Bg + g.r * Bf + f.r * g.r + trunc
    return TM(p, r)


def tm_exp_linear(w, N, T):
    """exp(w) for a TM w with zero constant term and zero remainder."""
    nv = len(next(iter(w.p)))
    out = tm_const(1, nv)
    term = tm_const(1, nv)
    for k in range(1, N + 1):
        term = tm_scale(tm_mul(term, w, N, T), arb(1) / k)
        out = tm_add(out, term)
    Bw = tm_bound(w, T)
    # the powers w^k with k > N are dropped: |sum_{k>N} w^k/k!| <= Bw^{N+1}/(N+1)! e^{Bw}; also the
    # truncation remainders accumulated inside term (already in out.r)
    out.r = out.r + ipow(Bw, N + 1) / arb(math.factorial(N + 1)) * Bw.exp()
    return out


def tm_log1p(y, N, T):
    By = tm_bound(y, T)
    assert By < 1, "log series: |y| bound >= 1"
    nv = len(next(iter(y.p)))
    out = TM({(0,) * nv: acb(0)}, arb(0))
    pw = tm_const(1, nv)
    for k in range(1, N + 1):
        pw = tm_mul(pw, y, N, T)
        out = tm_add(out, pw, 1, arb((-1) ** (k + 1)) / k)
    out.r = out.r + ipow(By, N + 1) / ((N + 1) * (1 - By))
    return out


# ---------------- the problem ----------------
def build_h(d, N, Ta):
    n = d - 1
    nv = n - 1
    om = [acb(arb(2 * k) / n).exp_pi_i() for k in range(n)]      # omega^k exactly enclosed
    eps = []
    for j in range(n):
        if j == 0:
            eps.append(TM({(0,) * nv: acb(0)}, arb(0)))
        else:
            e = [0] * nv; e[j - 1] = 1
            eps.append(TM({tuple(e): acb(1)}, arb(0)))
    E1 = [tm_exp_linear(eps[j], N, Ta) if j else tm_const(1, nv) for j in range(n)]
    Em1 = [tm_exp_linear(tm_scale(eps[j], -1), N, Ta) if j else tm_const(1, nv) for j in range(n)]
    c = arb(d - 1) / d
    hs = []
    for i in range(n):
        # z_j = omega^{i-j} exp(eps_i) exp(-eps_j), j != i
        zs = []
        for j in range(n):
            if j == i:
                continue
            z = tm_mul(E1[i], Em1[j], N, Ta)
            zs.append(tm_scale(z, om[(i - j) % n]))
        # elementary symmetric functions
        Ek = [tm_const(1, nv)] + [TM({(0,) * nv: acb(0)}, arb(0)) for _ in range(len(zs))]
        for z in zs:
            for k in range(len(zs), 0, -1):
                Ek[k] = tm_add(Ek[k], tm_mul(Ek[k - 1], z, N, Ta))
        S = TM({(0,) * nv: acb(0)}, arb(0))
        for k in range(len(zs) + 1):
            S = tm_add(S, Ek[k], 1, arb((-1) ** k) / ((k + 1) * (k + 2)))
        # y = S/c - 1
        y = tm_add(S, tm_const(1, nv), 1 / c, -1)
        h = tm_log1p(y, N, Ta)          # h_i - log c
        hs.append(h)
    return hs


def real_parts(h, nv):
    """gradient J_i (length 2nv) and Hessian matrix Q_i (Re quadratic part = 1/2 x^T Q x) as arb."""
    D = 2 * nv
    g = [arb(0)] * D
    Q = [[arb(0)] * D for _ in range(D)]
    for k, c in h.p.items():
        if sum(k) == 1:
            j = k.index(1)
            g[2 * j] += c.real; g[2 * j + 1] += -c.imag
        elif sum(k) == 2:
            idx = [j for j in range(nv) for _ in range(k[j])]
            j, l = idx
            # Re(c (x_j + i y_j)(x_l + i y_l)) = Re c (x_j x_l - y_j y_l) - Im c (x_j y_l + y_j x_l)
            a, b = c.real, c.imag
            terms = [((2 * j, 2 * l), a), ((2 * j + 1, 2 * l + 1), -a), ((2 * j, 2 * l + 1), -b), ((2 * j + 1, 2 * l), -b)]
            for (p, q), v in terms:
                # contributes v * x_p x_q to the quadratic form = 1/2 x^T Q x
                Q[p][q] += v; Q[q][p] += v
    return g, Q


def ldl_posdef(M):
    """True if the symmetric arb matrix M is certainly positive definite (LDL^T with all pivots > 0)."""
    D = len(M)
    A = [row[:] for row in M]
    for k in range(D):
        if not (A[k][k] > 0):
            return False
        for i in range(k + 1, D):
            f = A[i][k] / A[k][k]
            for j in range(k + 1, D):
                A[i][j] = A[i][j] - f * A[k][j]
    return True


def sphere_max(a):
    """max_{rho >= 0, |rho|_2 = 1} rho^a = prod_j (a_j/k)^{a_j/2}  (Lagrange), as an arb enclosure."""
    k = sum(a)
    v = arb(1)
    for aj in a:
        if aj:
            v *= ipow((arb(aj) / k).sqrt(), aj)
    return v


def opnorm_bound(M):
    """rigorous upper bound for the spectral norm of a (possibly non-square) arb matrix M (list of rows):
    lam >= ||M||^2 certified by LDL of lam I - M^T M."""
    Mm = np.array([[float(x.mid()) for x in r] for r in M])
    lam = float(np.linalg.norm(Mm, 2)) ** 2 * (1 + 1e-6) + 1e-20
    cols = len(M[0])
    MtM = [[sum((M[k][p] * M[k][q] for k in range(len(M))), arb(0)) for q in range(cols)] for p in range(cols)]
    A = [[(arb(lam) if p == q else arb(0)) - MtM[p][q] for q in range(cols)] for p in range(cols)]
    assert ldl_posdef(A), "opnorm LDL failed"
    return arb(lam).sqrt() * (1 + arb(1e-15))


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), arb(0)) for j in range(len(B[0]))] for i in range(len(A))]


def poly_real_quadratic(p, nv):
    """Q with Re(sum_{|a|=2} c_a eps^a) = 1/2 x^T Q x."""
    return real_parts(TM({k: c for k, c in p.items() if sum(k) == 2}, arb(0)), nv)[1]


def weighted_abs(p, k):
    s = arb(0)
    for a, c in p.items():
        if sum(a) == k:
            s += abs(c) * sphere_max(a)
    return s


def eval_re_poly(p, eps):
    """Re(sum_a c_a eps^a) for acb (interval) eps."""
    s = acb(0)
    for a, c in p.items():
        term = c
        for j, e in enumerate(a):
            if e:
                term = term * ipow(eps[j], e)
        s += term
    return s.real


def rp_mul(P, Q):
    R = {}
    for a, u in P.items():
        for b, v in Q.items():
            k = tuple(x + y for x, y in zip(a, b))
            R[k] = R[k] + u * v if k in R else u * v
    return R


def rp_add(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R[k] + s * v if k in R else s * v
    return R


def to_real_z_poly(p, V):
    """Re(sum_a c_a eps^a) with eps_j = x_{2j} + i x_{2j+1}, x = V z, as a real polynomial in z (arb coeffs)."""
    m = V.shape[1]
    nvar = V.shape[0] // 2
    def e(q):
        z = [0] * m; z[q] = 1; return tuple(z)
    X = [{e(q): arb(V[r, q]) for q in range(m)} for r in range(2 * nvar)]
    one = {(0,) * m: arb(1)}
    out = {}
    for a, c in p.items():
        re, im = one, {}
        for j, ex in enumerate(a):
            for _ in range(ex):
                xr, xi = X[2 * j], X[2 * j + 1]
                re, im = rp_add(rp_mul(re, xr), rp_mul(im, xi), -1), rp_add(rp_mul(re, xi), rp_mul(im, xr))
        out = rp_add(out, rp_add({k: v * c.real for k, v in re.items()}, {k: v * c.imag for k, v in im.items()}, -1))
    return out


def rp_eval(P, z):
    s = arb(0)
    for a, c in P.items():
        t = c
        for q, ex in enumerate(a):
            if ex:
                t = t * ipow(z[q], ex)
        s += t
    return s


def rp_diff(P, q):
    R = {}
    for a, c in P.items():
        if a[q]:
            b = list(a); b[q] -= 1
            R[tuple(b)] = R.get(tuple(b), arb(0)) + c * a[q]
    return R


def sphere_sup(polys, V, k, rel=1.05, h0=0.25, maxboxes=400000, seed=0):
    """Rigorous upper bound for sup_{z in R^m, |z| = 1} || (Re P_i(V z))_i ||_2 (P_i holomorphic polynomials
    in eps, eps_j = x_{2j} + i x_{2j+1}, x = V z).  Branch and bound over z-boxes meeting the unit sphere,
    mean-value form: |q(z)| <= |q(c)| + sum_q |d_q q(box)| h  (arb)."""
    m = V.shape[1]
    Rp = [to_real_z_poly(p, V) for p in polys]
    dR = [[rp_diff(P, q) for q in range(m)] for P in Rp]
    Rf = [{a: float(c.mid()) for a, c in P.items()} for P in Rp]
    def fval(z):
        return np.sqrt(sum(sum(c * np.prod(z ** np.array(a)) for a, c in P.items()) ** 2 for P in Rf))
    rng = np.random.default_rng(seed)
    zs = rng.normal(size=(20000, m)); zs /= np.linalg.norm(zs, axis=1)[:, None]
    vals = [fval(z) for z in zs]
    est = max(vals)
    z = zs[int(np.argmax(vals))]
    for step in [0.05, 0.02, 0.01, 0.005, 0.002, 0.001]:      # crude hill climb
        for _ in range(200):
            w = z + step * rng.normal(size=m); w /= np.linalg.norm(w)
            fw = fval(w)
            if fw > est:
                est, z = fw, w
    target = arb(est * rel + 1e-12)
    grid = np.arange(-1 + h0, 1, 2 * h0)
    stack = [(np.array(c), h0) for c in itertools.product(grid, repeat=m)]
    nb = 0
    while stack:
        c, h = stack.pop()
        # exact rational sphere-intersection test (c, h are dyadic floats: Fraction(float) is exact)
        cq = [Fraction(float(ci)) for ci in c]; hq = Fraction(float(h))
        lo2 = sum(max(abs(ci) - hq, 0) ** 2 for ci in cq)
        hi2 = sum((abs(ci) + hq) ** 2 for ci in cq)
        if lo2 > 1 or hi2 < 1:
            continue
        nb += 1
        if nb > maxboxes:
            return None, est
        zc = [arb(ci) for ci in c]
        zb = [arb(ci, h) for ci in c]
        tot = arb(0)
        for P, dP in zip(Rp, dR):
            v = abs(rp_eval(P, zc)) + sum((abs(rp_eval(dq, zb)) for dq in dP), arb(0)) * arb(h)
            tot += v * v
        if tot < target * target:
            continue
        for sgn in itertools.product((-0.5, 0.5), repeat=m):
            stack.append((c + np.array(sgn) * h, h / 2))
    return float(target.upper()), est


def certify(d=6, T="0.0527", N=6, sigma=5.0, verbose=True):
    """Certificate on E: for x in E with 0 < |x|_2 = t <= T:  L(x) < log c  (L = mean log|S_i| = log F on E).
    Steps (all constants rigorous, arb):
      Jx = -g(x) on E, g = G^{>=2}, G_i = l_i - L;  |Jx| <= kappa t^2  (fixed-point bound, see LOCAL_CERT.md)
      1/2 x^T H_L x = 1/2 x^T (H_L - sigma J^T J) x + sigma/2 |Jx|^2 <= -mu/2 t^2 + sigma/2 kappa^2 t^4
      L - log c <= that + sum_{k=3}^N B_k t^k + rem_L (t/T)^{N+1}."""
    n = d - 1; nv = n - 1; D = 2 * nv
    Ta = arb(str(T))                     # ball containing the exact decimal T (upper end >= T)
    assert Ta.lower() <= arb(Fraction(str(T)).numerator) / Fraction(str(T)).denominator <= Ta.upper()
    ex = exact_facts.check(d)            # S_i(0) = c and grad L(0) = 0, exactly (Q(omega))
    assert ex["ok"], ex
    hs = build_h(d, N, Ta)
    keys = set().union(*[h.p.keys() for h in hs])
    Lp = {k: sum((h.p.get(k, acb(0)) for h in hs), acb(0)) / n for k in keys}
    remL = sum((h.r for h in hs), arb(0)) / n
    Gp = [{k: h.p.get(k, acb(0)) - Lp[k] for k in keys} for h in hs]
    remG = [h.r + remL for h in hs]
    # exact-in-arb linear/quadratic parts
    J = [real_parts(h, nv)[0] for h in hs]                       # n x D
    Q = [real_parts(h, nv)[1] for h in hs]
    HL = [[sum((Q[i][p][q] for i in range(n)), arb(0)) / n for q in range(D)] for p in range(D)]
    A = [[[(Q[i][p][q] - HL[p][q]) / 2 for q in range(D)] for p in range(D)] for i in range(n)]  # G^2_i = x^T A_i x
    # gradient of L at 0 must vanish (it does exactly, by symmetry); report the enclosure
    gradL = [sum((J[i][k] for i in range(n)), arb(0)) / n for k in range(D)]
    # consistency of the enclosures with the exact facts (the exact values are what the proof uses)
    assert all(x.contains(0) for x in gradL)
    assert all(h.p[(0,) * nv].contains(0) for h in hs)
    # mu: LDL of -(H_L - sigma J^T J) - mu I
    sa = arb(sigma)
    JtJ = [[sum((J[i][p] * J[i][q] for i in range(n)), arb(0)) for q in range(D)] for p in range(D)]
    Hs = [[HL[p][q] - sa * JtJ[p][q] for q in range(D)] for p in range(D)]
    Hm = np.array([[float(x.mid()) for x in r] for r in Hs])
    mu = -float(np.linalg.eigvalsh(Hm).max()) * (1 - 1e-6)
    assert mu > 0
    assert ldl_posdef([[-Hs[p][q] - (arb(mu) if p == q else arb(0)) for q in range(D)] for p in range(D)])
    # projector data (float choices, rigorous norms)
    Jm = np.array([[float(x.mid()) for x in r] for r in J])
    Pp = np.linalg.pinv(Jm, rcond=1e-10)                          # D x n
    _, sv, Vt_ = np.linalg.svd(Jm)
    V = Vt_[n - 1:].T                                              # D x (D-n+1), approx. orthonormal basis of K
    V_a = [[arb(float(v)) for v in r] for r in V]                  # exact dyadic entries of the stored basis
    Pi_a = matmul(V_a, [list(col) for col in zip(*V_a)])          # Pi = V V^T in arb (the exact product)
    Pp_a = [[arb(v) for v in r] for r in Pp]
    beta = opnorm_bound(Pp_a)
    piN = opnorm_bound(Pi_a)
    PpJ = matmul(Pp_a, J)
    gam2 = sum((sq(arb(1 if p == q else 0) - Pi_a[p][q] - PpJ[p][q]) for p in range(D) for q in range(D)), arb(0))
    assert gam2.is_finite()
    gam = arb(gam2.upper()).sqrt()                 # sqrt of an exact upper bound of the (>= 0) sum of squares
    gam = arb(gam.upper()) + arb(1e-30)
    Jn = opnorm_bound(J)
    Vn = opnorm_bound([[arb(v) for v in r] for r in V])
    G2polys = [{k: c for k, c in Gp[i].items() if sum(k) == 2} for i in range(n)]
    qs, qest = sphere_sup(G2polys, V, 2)
    assert qs is not None
    qPi = arb(qs) * sq(Vn)                     # >= sup_{|v|<=1} |G^2(Pi v)|
    L3 = {k: c for k, c in Lp.items() if sum(k) == 3}
    cs, cest = sphere_sup([L3], V, 3)
    assert cs is not None
    cstar = arb(cs)
    C3raw = sum((abs(c) for c in L3.values()), arb(0))
    b2 = sum((sq(opnorm_bound(A[i])) for i in range(n)), arb(0)).sqrt()
    # higher-order bounds
    BL = {k: weighted_abs(Lp, k) for k in range(3, N + 1)}
    AG = [{k: weighted_abs(Gp[i], k) for k in range(3, N + 1)} for i in range(n)]
    g3 = sum((sq(sum((AG[i][k] * ipow(Ta, k) for k in range(3, N + 1)), arb(0)) + remG[i]) for i in range(n)), arb(0)).sqrt()
    # kappa fixed point:  y <= q t^2 + 2 b2 pi t (beta y + gam t) + b2 (beta y + gam t)^2 + g3(t),  a priori y <= y0
    num = qPi + 2 * b2 * piN * gam + b2 * sq(gam) + g3 / sq(Ta)
    kap = None
    y0_over = Jn / Ta          # y0 = |J| t  =>  y0 <= (|J|/T) t^2 * (T/t)...: handle first step with y0 = |J| t
    den = 1 - Ta * (2 * b2 * piN * beta + b2 * sq(beta) * Jn + 2 * b2 * beta * gam)
    hist = []
    if den > 0:
        kap = num / den
        hist.append(float(kap.upper()))
        for _ in range(30):
            den = 1 - 2 * b2 * piN * beta * Ta - b2 * sq(beta) * kap * sq(Ta) - 2 * b2 * beta * gam * Ta
            if not den > 0:
                break
            kap2 = num / den
            if kap2 < kap:
                kap = kap2
            hist.append(float(kap.upper()))
    res = dict(d=d, T=str(T), exact_facts=ex, N=N, sigma=sigma, mu=mu, gradL_max_abs=max(float(abs(x).upper()) for x in gradL),
               beta=float(beta.upper()), pi=float(piN.upper()), gamma=float(gam.upper()), normJ=float(Jn.upper()),
               q_Pi=float(qPi.upper()), b2=float(b2.upper()), g3_over_T2=float((g3 / sq(Ta)).upper()),
               B_L={k: float(v.upper()) for k, v in BL.items()}, remL=float(remL.upper()),
               q_sphere=qs, q_est=qest, cubic_sphere=cs, cubic_est=cest, C3raw=float(C3raw.upper()), normV=float(Vn.upper()),
               remG=[float(r.upper()) for r in remG], kappa_iter=hist[:3] + hist[-1:])
    if kap is None:
        res["certified"] = False; res["reason"] = "a-priori denominator <= 0"
    else:
        cub = cstar * ipow(Vn, 3) * Ta
        dcub = C3raw * Ta * (ipow(piN + beta * kap * Ta + gam, 3) - ipow(piN, 3))
        hig = sum((BL[k] * ipow(Ta, k - 2) for k in range(4, N + 1)), arb(0))
        g = -arb(mu) / 2 + sa / 2 * sq(kap) * sq(Ta) + cub + dcub + hig + remL / sq(Ta)
        res.update(kappa=float(kap.upper()), g_upper=float(g.upper()), certified=bool(g < 0),
                   terms=dict(quad=-mu / 2, penalty=float((sa / 2 * sq(kap) * sq(Ta)).upper()),
                              cubic_on_K=float(cub.upper()), cubic_offK=float(dcub.upper()),
                              cubic_crude_for_comparison=float((BL[3] * Ta).upper()),
                              higher=float(sum((BL[k] * ipow(Ta, k - 2) for k in range(4, N + 1)), arb(0)).upper()),
                              rem=float((remL / sq(Ta)).upper())))
    if verbose:
        print(json.dumps(res, indent=1), flush=True)
    return res


if __name__ == "__main__":
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    T = sys.argv[2] if len(sys.argv) > 2 else "0.0527"
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    sig = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    certify(d, T, N, sig)
