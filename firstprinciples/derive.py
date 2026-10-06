#!/usr/bin/env python3
"""derive.py -- first-principles checks of claims (2)-(5).

Everything is derived from:  P'(z) = prod_j (1 - z/b_j),  P(0) = 0,
V_i = P(b_i)/b_i.  P is built by multiplying the linear factors and
integrating the coefficients; no integral formula over t is assumed except
where the claimed formula itself is being tested.
"""
import itertools, math, random, sys, time
import sympy as sp
import numpy as np
import mpmath as mp
from fpcore import crit_values, crit_values_mp, crit_values_from_roots_np

NFAIL = 0


def report(ok, msg):
    global NFAIL
    if not ok: NFAIL += 1
    print(('PASS  ' if ok else 'FAIL  ') + msg, flush=True)


z, t, lam = sp.symbols('z t lambda')


def P_from_inverse_roots(v):
    """v_j = 1/b_j.  Returns P(z) as a sympy polynomial expression in z with
    P'(z) = prod (1 - v_j z), P(0) = 0, built by expansion + term-wise integration."""
    Pd = sp.Poly(sp.Integer(1), z)
    for vj in v:
        Pd = Pd * sp.Poly(1 - vj * z, z)
    coeffs = Pd.all_coeffs()[::-1]          # ascending
    return sum(sp.expand(c) * z ** (k + 1) / (k + 1) for k, c in enumerate(coeffs))


# ------------------------------------------------------------------ (2)
def check_formulas(d):
    n = d - 1
    u = sp.symbols('u2:%d' % (n + 1))           # u_2..u_n
    v = (sp.Integer(1),) + u                     # 1/b_j, with b_1 = 1
    P = P_from_inverse_roots(v)
    # V_1 = P(b_1)/b_1 with b_1 = 1
    V = [sp.expand(P.subs(z, 1))]
    for i in range(2, n + 1):
        ui = u[i - 2]
        V.append(sp.expand(P.subs(z, 1 / ui) * ui))          # P(b_i)/b_i, b_i = 1/u_i
    # claimed formulas
    S1 = sp.integrate(sp.expand((1 - t) * sp.prod([1 - t * uj for uj in u])), (t, 0, 1))
    ok = sp.expand(S1 - V[0]) == 0
    report(ok, 'd=%d (2) S_1 == V_1 exactly (P built from critical points)' % d)
    for i in range(2, n + 1):
        ui = u[i - 2]
        others = [u[j - 2] for j in range(2, n + 1) if j != i]
        Ti = sp.integrate(sp.expand((1 - t) * (ui - t) * sp.prod([ui - t * uj for uj in others])), (t, 0, 1))
        diff = sp.cancel(sp.together(Ti / ui ** (n - 1) - V[i - 1]))
        report(diff == 0, 'd=%d (2) S_%d == V_%d exactly' % (d, i, i))
    # numerical cross-check by quadrature against two independent evaluators
    rng = random.Random(d)
    worst = 0; worst_np = 0
    with mp.workdps(30):
        for _ in range(5):
            uu = [mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(n - 1)]
            Vmp = crit_values_mp(uu, 30)
            Vnp = crit_values_from_roots_np(np.array([1] + [1 / complex(x) for x in uu]))
            s1 = mp.quad(lambda tt: (1 - tt) * mp.fprod([1 - tt * x for x in uu]), [0, 1])
            worst = max(worst, abs(s1 - Vmp[0]))
            for i in range(2, n + 1):
                ui = uu[i - 2]; oth = [uu[j - 2] for j in range(2, n + 1) if j != i]
                Ti = mp.quad(lambda tt: (1 - tt) * (ui - tt) * mp.fprod([ui - tt * x for x in oth]), [0, 1])
                worst = max(worst, abs(Ti / ui ** (n - 1) - Vmp[i - 1]))
            worst_np = max(worst_np, max(abs(complex(a) - b) / max(1, abs(b)) for a, b in zip(Vmp, Vnp)))
    report(worst < 1e-20, 'd=%d (2) quadrature of S_i vs 30-digit critical values, max diff %.1e' % (d, worst))
    report(worst_np < 1e-10, 'd=%d (2) numpy.poly/polyint evaluator agrees with the mp evaluator, max rel diff %.1e' % (d, worst_np))


# ------------------------------------------------------------------ (3)
def V_general(b):
    """V_i as rational functions of general critical points b (sympy)."""
    P = P_from_inverse_roots([1 / bj for bj in b])
    return [sp.cancel(sp.together(P.subs(z, bi) / bi)) for bi in b]


def check_symmetries(d):
    n = d - 1
    b = sp.symbols('b1:%d' % (n + 1))
    V = V_general(b)
    # scale invariance: P_lam(z) = P(lam z)/lam has critical points b/lam, P_lam(0)=0, P_lam'(0)=1
    P = P_from_inverse_roots([1 / bj for bj in b])
    Pl = sp.expand(P.subs(z, lam * z) / lam)
    dPl = sp.diff(Pl, z)
    ok = all(sp.simplify(dPl.subs(z, bj / lam)) == 0 for bj in b) and Pl.subs(z, 0) == 0 \
        and sp.simplify(dPl.subs(z, 0) - 1) == 0
    report(ok, 'd=%d (3) P(lam z)/lam is normalised with critical points b_j/lam' % d)
    Vs = [sp.cancel(Vi.subs({bj: bj / lam for bj in b}, simultaneous=True) - Vi) for Vi in V]
    report(all(x == 0 for x in Vs), 'd=%d (3) V_i(b/lam) == V_i(b) for every i (complex lam != 0)' % d)
    # permutation equivariance: generators transposition (1 2) and n-cycle
    for name, sig in [('transposition(1 2)', [1, 0] + list(range(2, n))), ('n-cycle', list(range(1, n)) + [0])]:
        bp = [b[sig[j]] for j in range(n)]                       # b'_j = b_{sig(j)}
        ok = all(sp.cancel(V[i].subs(dict(zip(b, bp)), simultaneous=True) - V[sig[i]]) == 0 for i in range(n))
        report(ok, 'd=%d (3) V_i(b o sigma) == V_sigma(i)(b) for %s (so permutations of all n critical points permute the V_i)' % (d, name))
    # conjugation: V_i is a rational function with rational coefficients
    ok = True
    for Vi in V:
        num, den = sp.fraction(Vi)
        for pol in (sp.Poly(num, *b), sp.Poly(den, *b)):
            ok &= all(c.is_rational for c in pol.coeffs())
    report(ok, 'd=%d (3) V_i in Q(b_1..b_n), hence V_i(conj b) == conj V_i(b)' % d)
    # numerical: full reduction map applied to random configurations preserves the multiset {|V_i|}
    rng = np.random.default_rng(100 + d)
    worst = 0; inR = True
    for _ in range(2000):
        bb = (rng.normal(size=n) + 1j * rng.normal(size=n)) * np.exp(rng.normal(size=n))
        Vb = np.abs(crit_values_from_roots_np(bb))
        m = np.argmin(np.abs(bb))
        perm = [m] + [j for j in range(n) if j != m]
        bb2 = bb[perm] / bb[m]                  # b_1 = 1, |b_j| >= 1
        u = 1 / bb2[1:]
        u = u[np.argsort(-np.abs(u), kind='stable')]
        if u[0].imag < 0: u = np.conj(u)
        inR &= bool(np.all(np.abs(u) <= 1 + 1e-12) and np.all(np.diff(np.abs(u)) <= 1e-12) and u[0].imag >= 0)
        Vu = np.abs(crit_values(u[None, :])[0][0])
        worst = max(worst, np.max(np.abs(np.sort(Vu) - np.sort(Vb)) / np.sort(Vb)))
    report(inR and worst < 1e-8, 'd=%d (3) reduction (min-modulus b_1, scale, sort, conjugate) lands in R and preserves sorted |V| (max rel diff %.1e)' % (d, worst))


# ------------------------------------------------------------------ (4)
def check_equality_points(d, maxperms=None):
    n = d - 1
    c = sp.Rational(d - 1, d)
    w = sp.symbols('w')
    Phi = sp.Poly(sp.cyclotomic_poly(n, w), w)
    def red(e):
        return sp.Poly(sp.expand(e), w).rem(Phi).as_expr()
    perms = list(itertools.permutations(range(1, n)))
    if maxperms: perms = perms[:maxperms]
    allok = True
    for a in perms:
        A = (0,) + a                                  # b_1 = w^0 = 1, b_j = w^{a_j}
        # 1/b_j = w^{n - a_j}; u_j = conj(w^{a_j}) = w^{-a_j} = 1/b_j, consistent with b_j = 1/u_j
        v = [w ** ((n - aj) % n) for aj in A]
        Pd = sp.Poly(1, z)
        for vj in v:
            Pd = sp.Poly(sp.expand((Pd.as_expr()) * (1 - vj * z)), z)
        q = [red(cf) for cf in Pd.all_coeffs()[::-1]]
        for i in range(n):
            Vi = sum(q[k] * w ** ((A[i] * k) % n) / (k + 1) for k in range(n + 1))
            if red(Vi - c) != 0:
                allok = False
        # P' == 1 - z^n
        if not (q[0] == 1 and q[n] == -1 and all(qq == 0 for qq in q[1:n])):
            allok = False
    report(allok, "d=%d (4) all %d equality points: P' = 1 - z^%d and V_i = %s exactly in Q(w)/Phi_%d" % (d, len(perms), n, c, n))


# ------------------------------------------------------------------ (5)
def check_handover():
    # proof (see REPORT.md): eps_j = -Log(1 + w_j), w_j = u_j/p_j - 1, |w_j| = |u_j - p_j|;
    # |Log(1+w)| <= -log(1-|w|) = |w| g(|w|) with g(s) = -log(1-s)/s increasing on [0,1).
    # Numerical check:
    rng = np.random.default_rng(5)
    worst_ratio = 0
    for n in (3, 4, 5, 6):
        om = np.exp(2j * np.pi / n)
        p = np.array([np.conj(om ** (j - 1)) for j in range(2, n + 1)])
        for r in (1e-3, 0.01, 0.05, 1 / 20, 0.2, 0.5, 0.9):
            for trial in range(4000):
                wv = rng.normal(size=n - 1) + 1j * rng.normal(size=n - 1)
                if trial % 4 == 0:                    # concentrate on one coordinate
                    wv = np.zeros(n - 1, complex); wv[rng.integers(n - 1)] = np.exp(1j * rng.uniform(0, 2 * np.pi))
                if trial % 4 == 1:
                    wv = -np.abs(wv)                  # real negative directions (worst case)
                wv = wv / np.linalg.norm(wv) * r * rng.uniform(0.0, 1.0) ** (1 / (2 * n))
                u = p * (1 + wv)                      # |u - p| = |w| since |p_j| = 1
                assert abs(np.linalg.norm(u - p) - np.linalg.norm(wv)) < 1e-12
                b = 1 / u
                eps = np.log(b / np.array([om ** (j - 1) for j in range(2, n + 1)]))   # principal branch
                worst_ratio = max(worst_ratio, np.linalg.norm(eps) / (-math.log(1 - r)))
    report(worst_ratio <= 1 + 1e-12, '(5) ||eps||_2 <= -log(1-r) whenever ||u - p||_2 <= r (sampled, max ratio %.12f)' % worst_ratio)
    # sharpness: one coordinate with w = -r, i.e. u_2 = (1 - r) p_2, all other u_j = p_j
    r = 1 / 20
    n = 5; om = np.exp(2j * np.pi / n)
    p = np.array([np.conj(om ** (j - 1)) for j in range(2, n + 1)])
    u = p.copy(); u[0] *= (1 - r)
    eps = np.log((1 / u) / np.array([om ** (j - 1) for j in range(2, n + 1)]))
    report(abs(np.linalg.norm(eps) - (-math.log(1 - r))) < 1e-15,
           '(5) bound is attained at u_2 = (1-r) p_2: ||eps|| = %.15f = -log(1-1/20)' % np.linalg.norm(eps))
    # relabelling: permutation sigma of indices 2..n maps p_a to p_id, is an isometry of R^D,
    # and permutes (V_2..V_n) (V_1 fixed).
    ok = True; worst = 0
    for n in (3, 4, 5, 6):
        om = np.exp(2j * np.pi / n)
        for a in itertools.permutations(range(1, n)):
            pa = np.array([np.conj(om ** aj) for aj in a])
            # sigma: new index j (value a = j-1) <- old position where a == j-1
            sig = [a.index(j - 1) for j in range(2, n + 1)]
            pid = np.array([np.conj(om ** (j - 1)) for j in range(2, n + 1)])
            ok &= np.allclose(pa[sig], pid)
            for _ in range(3):
                x = pa + 0.05 * (rng.normal(size=n - 1) + 1j * rng.normal(size=n - 1)) / 3
                y = x[sig]
                ok &= abs(np.linalg.norm(y - pid) - np.linalg.norm(x - pa)) < 1e-14
                Vx = crit_values(x[None, :])[0][0]; Vy = crit_values(y[None, :])[0][0]
                worst = max(worst, abs(Vx[0] - Vy[0]), np.max(np.abs(Vx[1:][sig] - Vy[1:])))
    report(ok and worst < 1e-12, '(5) relabelling b_2..b_n maps the ball around each equality point isometrically onto the ball around the identity point and permutes V_2..V_n (max diff %.1e)' % worst)
    # r_excl as stored in the task headers
    r_hdr = float.fromhex('0x1.999999999999ap-5')
    from fractions import Fraction
    print('INFO  task-header r_excl = 0x3fa999999999999a = %.20f; exceeds 1/20 by %.3e (exact)' % (r_hdr, float(Fraction(r_hdr) - Fraction(1, 20))))


if __name__ == '__main__':
    t0 = time.time()
    for d in (4, 5, 6, 7):
        check_formulas(d)
    for d in (4, 5, 6, 7):
        check_symmetries(d)
    for d in (4, 5, 6, 7):
        check_equality_points(d)
    check_handover()
    print('total failures: %d  (%.0fs)' % (NFAIL, time.time() - t0))
