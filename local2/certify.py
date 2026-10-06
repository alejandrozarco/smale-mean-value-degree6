"""Rigorous local certificate near the regular configuration, d in {5, 6, 7} (see README.md for the mathematics).

Usage:  python3 certify.py [--d=D] [--redesign] [--joint] [--K=K] [T] [N] [gamma]
    --d=D       polynomial degree d (default 6); n = d - 1 critical points, c = (d-1)/d, omega = exp(2 pi i/n)
    --redesign  ignore a saved design file and recompute (and overwrite) it
    --joint     also run the sharper final step joint_certificate() (K intervals in s, default K = 64)
    T           radius in the x-chart, a rational like 1/19 (default 1/19)
    N           Taylor truncation degree (default 6)
    gamma       design parameter (float) for the multipliers (default: automatic)

For d = 6 the saved designs are design_T<p>_<q>.json (as before); for d != 6 they are design_d<D>_T<p>_<q>_N<N>.json.

Everything that enters a final inequality is either exact (Fractions / exact Q(omega)) or an arb ball;
floating point (numpy) is used only to CHOOSE the multipliers Gamma, Psi and the number lambda, which are then
converted exactly to arb and used as given.
"""
import os
import sys
import time
from fractions import Fraction
from math import factorial, comb

import numpy as np
from flint import arb, acb, fmpq, ctx

import exact_qomega as EQ

ctx.prec = 200
D = 6           # polynomial degree (set by setup())
NPTS = 5        # n = d - 1 critical points
NM = 4          # number of constraints D_a / multipliers mu_a = n - 1
NV = 8          # real variables x_1..x_{2(n-1)}
C = Fraction(5, 6)
BASE = 16       # key = sum e_k * 16^k ; all degrees stay < 16
t0 = time.time()


def log(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------------------------------------
# exponent keys for real monomials in x (8 variables)
# ----------------------------------------------------------------------------------------------------------
def key_of(e):
    k = 0
    for i in reversed(range(NV)):
        k = k * BASE + e[i]
    return k


def exps_of(k):
    e = []
    for _ in range(NV):
        e.append(k % BASE)
        k //= BASE
    return tuple(e)


_deg_cache = {}


def deg(k):
    d = _deg_cache.get(k)
    if d is None:
        d = sum(exps_of(k))
        _deg_cache[k] = d
    return d


def unit_key(i):
    return BASE ** i


# ----------------------------------------------------------------------------------------------------------
# real polynomials: dict key -> arb
# ----------------------------------------------------------------------------------------------------------
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r[k] + s * v if k in r else s * v
    return r


def pscale(p, a):
    return {k: a * v for k, v in p.items()}


def by_degree(p):
    out = {}
    for k, v in p.items():
        out.setdefault(deg(k), []).append((k, v))
    return out


def pmul(p, q, Nt):
    """Product truncated to total degree <= Nt."""
    qd = by_degree(q)
    r = {}
    for ka, va in p.items():
        da = deg(ka)
        for db in range(0, Nt - da + 1):
            for kb, vb in qd.get(db, ()):
                k = ka + kb
                if k in r:
                    r[k] += va * vb
                else:
                    r[k] = va * vb
    return r


def homog(p, m):
    return {k: v for k, v in p.items() if deg(k) == m}


def frob(p, m):
    """Upper bound (arb, exact endpoint) for sup_{|u|=1} |p_m(u)|, p_m the degree-m part of p:
    |p_m(u)| <= sqrt( sum_alpha p_alpha^2 alpha!/m! )   (Cauchy-Schwarz with the multinomial identity)."""
    s = arb(0)
    fm = factorial(m)
    for k, v in p.items():
        if deg(k) != m:
            continue
        af = 1
        for a in exps_of(k):
            af *= factorial(a)
        s += v * v * arb(fmpq(af, fm))
    return _sqrt_up(s)


def _sqrt_up(s):
    """Upper bound for sqrt of a (true) nonnegative quantity enclosed by the ball s."""
    u = s.upper()
    if u.is_zero():
        return arb(0)
    assert u > 0
    return u.sqrt().upper()


def maxdeg(p):
    return max((deg(k) for k in p), default=0)


# ----------------------------------------------------------------------------------------------------------
# Step 1: exact Taylor coefficients of S_i (in Q(omega)) and their arb enclosures, split into real polys in x
# ----------------------------------------------------------------------------------------------------------
def omega_pows():
    return [acb(arb.cos_pi_fmpq(fmpq(2 * n, NPTS)), arb.sin_pi_fmpq(fmpq(2 * n, NPTS))) for n in range(NPTS)]


W = omega_pows()


def setup(d):
    """Select the degree d (5, 6 or 7).  Must be called before certify() when d != 6."""
    global D, NPTS, NM, NV, C, W
    EQ.setup(d)
    D, NPTS = d, d - 1
    NM, NV = NPTS - 1, 2 * (NPTS - 1)
    C = Fraction(NPTS, d)
    W = omega_pows()
    _deg_cache.clear()


def qomega_to_acb(v):
    s = acb(0)
    for n in range(NPTS):
        if v[n] != 0:
            s += W[n] * arb(fmpq(v[n].numerator, v[n].denominator))
    return s


def complex_series_to_real(cs):
    """cs: dict beta((n-1)-tuple) -> acb coefficient of eps^beta, eps_k = x_{2k-1} + i x_{2k}.
    Returns (U, V): real polys with U + iV = sum c_beta eps^beta."""
    U, V = {}, {}
    I_POW = [acb(1), acb(0, 1), acb(-1), acb(0, -1)]
    for beta, c in cs.items():
        # expand prod_k (xr_k + i xi_k)^{beta_k}
        terms = [((), 1, 0)]  # (partial exps list, integer multiplier, power of i)
        for kk, bk in enumerate(beta):
            new = []
            for ex, mul, ip in terms:
                for j in range(bk + 1):   # j = power of the imaginary coordinate
                    new.append((ex + (bk - j, j), mul * comb(bk, j), ip + j))
            terms = new
        for ex, mul, ip in terms:
            val = c * I_POW[ip % 4] * mul
            k = key_of(ex)
            re, im = val.real, val.imag
            U[k] = U[k] + re if k in U else re
            V[k] = V[k] + im if k in V else im
    return U, V


def cfrob(cs, m):
    """Upper bound for sup_{|eps|=1} |p_m(eps)| for complex homogeneous part of degree m (n-1 complex vars)."""
    s = arb(0)
    fm = factorial(m)
    for beta, c in cs.items():
        if sum(beta) != m:
            continue
        bf = 1
        for b in beta:
            bf *= factorial(b)
        s += (c.real * c.real + c.imag * c.imag) * arb(fmpq(bf, fm))
    return _sqrt_up(s)


def tail_exp(r, N):
    """Upper bound for sum_{m>N} r^m/m!  (0 <= r < N+2):  r^{N+1}/(N+1)! / (1 - r/(N+2))."""
    assert r < N + 2
    p = arb(1)
    for _ in range(N + 1):
        p = p * r
    return (p / factorial(N + 1) / (1 - r / (N + 2))).upper()


def certify(Tq, N=6, gamma_arg=None, qweight=None, verbose=True, design_out=None, design_in=None):
    """Run the certificate for radius Tq (a Fraction).  Returns a dict with the certified constants."""
    global log
    if not verbose:
        log = lambda *a: None
    Nt = N
    T = arb(fmpq(Tq.numerator, Tq.denominator))
    c = C
    n = NPTS
    log(f"=== Local certificate, d={D}: T = {Tq} (= {float(Tq):.6f}), Taylor degree N = {N}, prec = {ctx.prec} bits")

    # ---------------- exact facts
    log("[1] exact facts at x = 0 (exact arithmetic in Q(omega))")
    ok, _ = EQ.exact_facts(verbose=verbose)
    assert ok, "exact facts failed"

    # ---------------- Taylor polynomials of sigma_i = S_i / c
    log("[2] exact Taylor coefficients of S_i up to degree", N, "-> arb; real polynomials U_i + i V_i = S_i/c")
    sig_c = []      # complex coefficient dicts of sigma_i
    UV = []
    tau = []        # tail bounds tau_i(T) for |S_i/c - sigma_i| on |x| <= T
    for i in range(1, n + 1):
        tc = EQ.taylor_S(i, N)
        cs = {}
        for beta, v in tc.items():
            vv = [x / c for x in v]
            if sum(beta) == 0:
                assert EQ.qeq(vv, EQ.unit(1))
                cs[beta] = acb(1)            # exact: S_i(0)/c = 1 (fact F1)
            else:
                cs[beta] = qomega_to_acb(vv)
        sig_c.append(cs)
        UV.append(complex_series_to_real(cs))
        # tail: |sum_{m>N} P_{i,m}(eps)| <= sum_J |coef_J| tail_exp(|lambda_J| T, N)
        tsum = arb(0)
        for coef, e, lam in EQ.subsets_data(i):
            nl = arb(sum(l * l for l in lam)).sqrt().upper()
            if sum(l * l for l in lam) == 0:
                continue            # constant exponential: no Taylor tail
            tsum += arb(fmpq(abs(coef.numerator), coef.denominator)) * tail_exp(nl * T, N)
        tau.append((tsum / arb(fmpq(c.numerator, c.denominator))).upper())
    log("    tail bounds tau_i(T):", [t.str(5) for t in tau])

    # complex Frobenius norms n_{i,m}
    nrm = [[cfrob(sig_c[i], m) if m > 0 else arb(1) for m in range(N + 1)] for i in range(n)]

    # ---------------- A_i = |S_i/c|^2 : explicit part alpha_i (deg <= Nt) + error E_i(T)
    log("[3] alpha_i = [U_i^2 + V_i^2]_{<=Nt} and error bounds E_i(T)")
    alpha, E = [], []
    for i in range(n):
        U, V = UV[i]
        a = padd(pmul(U, U, Nt), pmul(V, V, Nt))
        # exact constant term 1 (fact F1)
        a[0] = arb(1)
        alpha.append(a)
        hi = arb(0)
        for m in range(N + 1):
            for mp in range(N + 1):
                if m + mp > Nt:
                    p = nrm[i][m] * nrm[i][mp]
                    for _ in range(m + mp):
                        p = p * T
                    hi += p
        sabs = arb(0)
        for m in range(N + 1):
            p = nrm[i][m]
            for _ in range(m):
                p = p * T
            sabs += p
        Ei = (hi + 2 * sabs * tau[i] + tau[i] * tau[i]).upper()
        E.append(Ei)
    log("    E_i(T):", [e.str(5) for e in E])

    # Phi = (1/n) sum alpha_i - 1 ; D_i = alpha_i - alpha_1
    phi = {}
    for i in range(n):
        phi = padd(phi, pscale(alpha[i], arb(fmpq(1, n))))
    phi[0] = phi.get(0, arb(0)) - 1
    # constant and linear terms of Phi are exactly 0 (facts F1, F2); check the enclosures contain 0 and set them
    for k in list(phi):
        if deg(k) <= 1:
            assert phi[k].contains(0), ("enclosure of an exactly-zero coefficient does not contain 0", k, phi[k])
            del phi[k]
    d = []
    for i in range(1, n):
        di = padd(alpha[i], alpha[0], -1)
        assert di[0].contains(0)
        del di[0]           # D_i(0) = 0 exactly (F1)
        d.append(di)

    # ---------------- design (floating point, non-rigorous choices)
    log("[4] design of multipliers (floating point; choices only)")

    def fl(v):
        return float(v.mid())

    P = np.zeros((NV, NV))
    for k, v in homog(phi, 2).items():
        ex = exps_of(k)
        idx = [j for j in range(NV) for _ in range(ex[j])]
        if idx[0] == idx[1]:
            P[idx[0], idx[0]] = fl(v)
        else:
            P[idx[0], idx[1]] = P[idx[1], idx[0]] = fl(v) / 2
    Lam = np.zeros((NM, NV))
    for a in range(NM):
        for k, v in homog(d[a], 1).items():
            Lam[a, exps_of(k).index(1)] = fl(v)
    _, sv, Vt = np.linalg.svd(Lam)
    Kb, Kp = Vt[NM:].T, Vt[:NM].T
    M = Lam @ Kp
    Pzz, Pzw, Pww = Kb.T @ P @ Kb, Kb.T @ P @ Kp, Kp.T @ P @ Kp
    lamK = -np.linalg.eigvalsh(Pzz).max()
    log(f"    singular values of Lambda: {sv}")
    log(f"    Phi_2 on ker Lambda: max eigenvalue {-lamK:.6f}")
    gamma = gamma_arg if gamma_arg is not None else 1.1 * lamK
    Gz = -2 * np.linalg.inv(M).T @ Pzw.T
    Gw = -np.linalg.inv(M).T @ (gamma * np.eye(NM) + Pww)
    Gam = Gz @ Kb.T + Gw @ Kp.T          # NM x NV, mu_lin(x) = Gam x
    Gam = np.round(Gam * 2 ** 40) / 2 ** 40

    # Psi by weighted least squares on the cubic (and quartic) parts of F
    def fpoly(p):
        return {k: fl(v) for k, v in p.items()}

    def fmul(p, q, Nt_):
        r = {}
        for ka, va in p.items():
            for kb, vb in q.items():
                if deg(ka) + deg(kb) <= Nt_:
                    r[ka + kb] = r.get(ka + kb, 0.0) + va * vb
        return r

    lin_keys = [unit_key(j) for j in range(NV)]
    quad_keys = sorted({unit_key(a) + unit_key(b) for a in range(NV) for b in range(a, NV)})
    cub_keys = sorted({q + unit_key(j) for q in quad_keys for j in range(NV)})
    quart_keys = sorted({q + r for q in quad_keys for r in quad_keys})

    def weight(k):
        m = deg(k)
        af = 1
        for a in exps_of(k):
            af *= factorial(a)
        return np.sqrt(af / factorial(m))

    fphi = fpoly(phi)
    fd = [fpoly(x) for x in d]
    gam_polys = [{lin_keys[j]: Gam[a, j] for j in range(NV)} for a in range(NM)]
    base = dict(fphi)
    for a in range(NM):
        for k, v in fmul(gam_polys[a], fd[a], Nt).items():
            base[k] = base.get(k, 0.0) + v
    ci = {k: n for n, k in enumerate(cub_keys)}
    qi = {k: n for n, k in enumerate(quart_keys)}
    nunk = NM * len(quad_keys)
    A3 = np.zeros((len(cub_keys), nunk))
    A4 = np.zeros((len(quart_keys), nunk))
    for a in range(NM):
        d1 = homog(fd[a], 1)
        d2 = homog(fd[a], 2)
        for qn, qk in enumerate(quad_keys):
            col = a * len(quad_keys) + qn
            for k, v in d1.items():
                A3[ci[qk + k], col] += v
            for k, v in d2.items():
                A4[qi[qk + k], col] += v
    b3 = np.array([base.get(k, 0.0) for k in cub_keys])
    b4 = np.array([base.get(k, 0.0) for k in quart_keys])
    w3 = np.array([weight(k) for k in cub_keys])
    w4 = np.array([weight(k) for k in quart_keys])
    tf = float(Tq) if qweight is None else qweight
    Abig = np.vstack([w3[:, None] * A3, tf * w4[:, None] * A4])
    bbig = np.concatenate([w3 * b3, tf * w4 * b4])
    psi, *_ = np.linalg.lstsq(Abig, -bbig, rcond=None)
    psi = np.round(psi * 2 ** 40) / 2 ** 40
    log(f"    gamma = {gamma:.5f};  |Gamma|_F = {np.linalg.norm(Gam):.4f};  |psi| = {np.linalg.norm(psi):.4f}")
    log(f"    float |G3|_F before/after Psi: {np.linalg.norm(w3 * b3):.5f} / {np.linalg.norm(w3 * (b3 + A3 @ psi)):.5f}")

    if design_in:
        # reuse a saved design (exact doubles), so the certified constants do not depend on the numpy/BLAS build
        import json
        with open(design_in) as fh:
            dz = json.load(fh)
        assert dz["T"] == str(Tq) and dz["N"] == N, "saved design is for another T or N"
        assert dz.get("d", 6) == D, "saved design is for another d"
        Gam = np.array([[float.fromhex(v) for v in row] for row in dz["Gamma"]])
        psi = np.array([float.fromhex(dz["Psi"][a][str(exps_of(qk))]) for a in range(NM) for qk in quad_keys])
        gamma = dz["gamma"]
        log(f"    design loaded from {design_in} (Gamma, Psi as saved; the float design above is not used)")

    # ---------------- rigorous assembly of F = Phi + sum_a mu_a D_a
    log("[5] rigorous assembly of the explicit part G of F = Phi + sum_a mu_a(x) D_a(x)")
    Gam_a = [{lin_keys[j]: arb(float(Gam[a, j])) for j in range(NV) if Gam[a, j] != 0} for a in range(NM)]
    Psi_a = [{qk: arb(float(psi[a * len(quad_keys) + qn])) for qn, qk in enumerate(quad_keys)
              if psi[a * len(quad_keys) + qn] != 0} for a in range(NM)]
    G = dict(phi)
    for a in range(NM):
        mu = padd(Gam_a[a], Psi_a[a])
        G = padd(G, pmul(mu, d[a], Nt))
    # G has no terms of degree 0 or 1: phi has none (exact), mu*d starts at degree 2.
    assert all(deg(k) >= 2 for k in G)

    # ---------------- G_2 <= -lambda |x|^2 via an arb Cholesky of -H2 - lambda I
    H2 = [[arb(0)] * NV for _ in range(NV)]
    for k, v in homog(G, 2).items():
        ex = exps_of(k)
        idx = [j for j in range(NV) for _ in range(ex[j])]
        if idx[0] == idx[1]:
            H2[idx[0]][idx[0]] = v
        else:
            H2[idx[0]][idx[1]] = v / 2
            H2[idx[1]][idx[0]] = v / 2
    H2f = np.array([[fl(v) for v in row] for row in H2])
    ev = np.linalg.eigvalsh(H2f)
    log(f"    float eigenvalues of H2 (quadratic part of F): {np.round(ev, 5)}")
    lam_f = -ev.max()
    lam = Fraction(int(lam_f * 0.999 * 2 ** 30), 2 ** 30)
    lam_a = arb(fmpq(lam.numerator, lam.denominator))
    Mx = [[(-H2[i][j] - (lam_a if i == j else 0)) for j in range(NV)] for i in range(NV)]
    Lc = [[arb(0)] * NV for _ in range(NV)]
    pd = True
    for j in range(NV):
        s = Mx[j][j]
        for k in range(j):
            s -= Lc[j][k] * Lc[j][k]
        if not (s > 0):
            pd = False
            break
        Lc[j][j] = s.sqrt()
        for i in range(j + 1, NV):
            s2 = Mx[i][j]
            for k in range(j):
                s2 -= Lc[i][k] * Lc[j][k]
            Lc[i][j] = s2 / Lc[j][j]
    log(f"    certified: G_2(x) <= -lambda |x|^2 with lambda = {lam} = {float(lam):.8f}: {pd}")
    assert pd

    # ---------------- higher homogeneous parts and error terms
    log("[6] bounds on higher-order terms at s = T")
    Gn = {m: frob(G, m) for m in range(3, Nt + 1)}
    budget = arb(0)
    for m in range(3, Nt + 1):
        p = Gn[m]
        for _ in range(m - 2):
            p = p * T
        budget += p
        log(f"    |G_{m}|_F = {Gn[m].str(6)}   contribution |G_{m}| T^{m-2} = {p.str(6)}")
    # errors
    errF = arb(0)
    for i in range(n):
        errF += E[i] / n
    for a in range(NM):
        ng = _sqrt_up(sum((arb(float(Gam[a, j])) * arb(float(Gam[a, j])) for j in range(NV)), arb(0)))
        npsi = frob(Psi_a[a], 2)
        mu_T = ng * T + npsi * T * T
        errF += mu_T * (E[a + 1] + E[0])
        dNt = frob(d[a], Nt)
        dNt1 = frob(d[a], Nt - 1)
        TN1 = arb(1)
        for _ in range(Nt + 1):
            TN1 = TN1 * T
        errF += ng * dNt * TN1 + npsi * (dNt1 * TN1 + dNt * TN1 * T)
    errF = errF.upper()
    err_term = (errF / (T * T)).upper()
    log(f"    total error bound errF(T) = {errF.str(6)};  errF(T)/T^2 = {err_term.str(6)}")
    total = (lam_a - budget - err_term)
    log(f"    lambda - sum_m |G_m| T^(m-2) - errF/T^2 = {total.str(8)}")
    ok = total > 0
    log(f"[7] RESULT: F(x) <= -({(total/1).lower().str(8)}) |x|^2 for all |x| <= T : {ok}")
    kap = (total / 2).lower()
    log(f"    => on E with |x| <= T:  Phi <= -{(total).lower().str(8)} |x|^2,  "
        f"L - log({c}) <= Phi/2 <= -kappa |x|^2 with kappa = {kap.str(8)}")
    res = dict(T=Tq, N=N, ok=bool(ok), lam=lam, total=total, kappa=kap, budget=budget, err_term=err_term,
               Gn=Gn, G=G, phi=phi, d=d, Gam=Gam, psi=psi, quad_keys=quad_keys, E=E, tau=tau, errF=errF,
               gamma=gamma)
    if design_out:
        import json
        with open(design_out, "w") as fh:
            json.dump({"d": D, "T": str(Tq), "N": N, "gamma": gamma, "lambda": str(lam),
                       "Gamma": [[float(v).hex() for v in row] for row in Gam],
                       "Psi": [{str(exps_of(qk)): float(psi[a * len(quad_keys) + qn]).hex()
                                for qn, qk in enumerate(quad_keys)} for a in range(NM)]}, fh, indent=1)
    return res


# ----------------------------------------------------------------------------------------------------------
# Optional sharper final step ("joint" bound, flag --joint; not used by the default d = 6 run)
# ----------------------------------------------------------------------------------------------------------
def sym_basis(k):
    """Exponent tuples alpha (|alpha| = k) in NV variables: an orthonormal basis of Sym^k(R^NV).
    The coordinates of u^{(x)k} in this basis are y_alpha = sqrt(k!/alpha!) u^alpha, and |y| = |u|^k."""
    return [exps_of(kk) for kk in sorted({sum(unit_key(j) for j in idx)
                                          for idx in _combs_with_rep(range(NV), k)})]


def _combs_with_rep(it, k):
    from itertools import combinations_with_replacement
    return combinations_with_replacement(list(it), k)


def _fact_t(a):
    f = 1
    for x in a:
        f *= factorial(x)
    return f


def sym_tensor_entry(G, alpha, m):
    """T_alpha = p_alpha alpha!/m!  (entry of the symmetric tensor of the degree-m part of G)."""
    v = G.get(key_of(alpha))
    if v is None:
        return None
    return v * arb(fmpq(_fact_t(alpha), factorial(m)))


def _sqrt_fq(p, q):
    return arb(fmpq(p, q)).sqrt()


def flattening(G, m, p):
    """Matrix A (dim Sym^p x dim Sym^(m-p), arb) with G_m(u) = y^T A w, y = u^{(x)p}, w = u^{(x)(m-p)} in the
    orthonormal symmetric bases:  A_{beta,gamma} = sqrt(p!/beta!) sqrt(q!/gamma!) T_{beta+gamma}."""
    q = m - p
    Bp, Bq = sym_basis(p), sym_basis(q)
    A = [[arb(0)] * len(Bq) for _ in Bp]
    for i, be in enumerate(Bp):
        sb = _sqrt_fq(factorial(p), _fact_t(be))
        for j, ga in enumerate(Bq):
            al = tuple(a + b for a, b in zip(be, ga))
            t = sym_tensor_entry(G, al, m)
            if t is not None:
                A[i][j] = sb * _sqrt_fq(factorial(q), _fact_t(ga)) * t
    return A


def cert_lmax(Mx, lam_guess_up):
    """Certify Mx <= mu I (Mx symmetric arb matrix) by a ball Cholesky of mu I - Mx; mu = lam_guess_up (float)
    is converted exactly.  Returns mu as an arb (exact) or None."""
    n = len(Mx)
    mu = arb(lam_guess_up)
    Lc = [[arb(0)] * n for _ in range(n)]
    for j in range(n):
        s = mu - Mx[j][j]
        for k in range(j):
            s -= Lc[j][k] * Lc[j][k]
        if not (s > 0):
            return None
        Lc[j][j] = s.sqrt()
        for i in range(j + 1, n):
            s2 = -Mx[i][j]
            for k in range(j):
                s2 -= Lc[i][k] * Lc[j][k]
            Lc[i][j] = s2 / Lc[j][j]
    return mu


def _up_guess(x):
    """A float slightly above the float x (for certification targets)."""
    return float(x + 1e-9 * max(1.0, abs(x)) + 1e-12)


def flat_norm(G, m, p=2):
    """Certified upper bound for sup_{|u|=1} |G_m(u)| <= ||A||_2 (flattening p x (m-p)), via AA^T <= c I."""
    A = flattening(G, m, p)
    nr = len(A)
    AAt = [[arb(0)] * nr for _ in range(nr)]
    nz = [[(j, v) for j, v in enumerate(row) if not v.is_zero()] for row in A]
    cols = {}
    for i, row in enumerate(nz):
        for j, v in row:
            cols.setdefault(j, []).append((i, v))
    for j, lst in cols.items():
        for a in range(len(lst)):
            i1, v1 = lst[a]
            for b in range(a, len(lst)):
                i2, v2 = lst[b]
                AAt[i1][i2] += v1 * v2
    for i1 in range(nr):
        for i2 in range(i1):
            AAt[i1][i2] = AAt[i2][i1]
    f = np.array([[float(v.mid()) for v in row] for row in AAt])
    c = _up_guess(np.linalg.eigvalsh(f).max())
    mu = cert_lmax(AAt, c)
    assert mu is not None
    return mu.sqrt().upper()


def joint_certificate(res, K=64):
    """Sharper bound for sup_{|u|=1} G_2(u) + s G_3(u) + s^2 G_4(u), s in [0, T].

    With y = u (x) u in the orthonormal basis of Sym^2 (|y| = 1) and z = (u, y) (|z|^2 = 2):
        G_2(u) = G_2(u)|u|^2 = y^T Hb y,  Hb = matrix of X -> (H X + X H)/2 on Sym^2 (H = Hessian/2 of G_2),
        G_3(u) = u^T A3 y,   G_4(u) = y^T B4 y   (flattenings of the symmetric tensors),
    so for every real t:  G_2 + s G_3 + s^2 G_4 = z^T M(s,t) z  with
        M(s,t) = [[t I, s A3/2], [s A3^T/2, Hb + s^2 B4 - t I]],
    hence <= 2 lambda_max(M(s,t)).  [0, T] is cut into K intervals [s_m - rho, s_m + rho]; M(s_m, t_k) <= mu_k I is
    certified by a ball Cholesky, and M(s) = M(s_m) + delta M1 + (2 s_m delta + delta^2) M2 with ||M1|| <= |G_3|_F/2,
    ||M2|| <= |G_4|_F covers the rest of the interval.  Degrees m >= 5 use certified flattening norms."""
    G, Tq, N = res["G"], res["T"], res["N"]
    T = arb(fmpq(Tq.numerator, Tq.denominator))
    log(f"[6b] joint bound on G_2 + s G_3 + s^2 G_4 over the sphere (K = {K} intervals in s), flattening norms for m >= 5")
    tj = time.time()
    B2 = sym_basis(2)
    ns = len(B2)
    # H: symmetric matrix of G_2
    H = [[arb(0)] * NV for _ in range(NV)]
    for k, v in homog(G, 2).items():
        ex = exps_of(k)
        idx = [j for j in range(NV) for _ in range(ex[j])]
        if idx[0] == idx[1]:
            H[idx[0]][idx[0]] = v
        else:
            H[idx[0]][idx[1]] = v / 2
            H[idx[1]][idx[0]] = v / 2
    # orthonormal basis of Sym^2 as sparse matrices: beta = 2e_i -> E_ii ; e_i + e_j -> (E_ij + E_ji)/sqrt 2
    r2 = arb(2).sqrt()
    E = []
    for be in B2:
        idx = [j for j in range(NV) for _ in range(be[j])]
        if idx[0] == idx[1]:
            E.append([(idx[0], idx[0], arb(1))])
        else:
            E.append([(idx[0], idx[1], 1 / r2), (idx[1], idx[0], 1 / r2)])
    # Hb_{b,b'} = <E_b, (H E_b' + E_b' H)/2> = (tr(E_b H E_b') + tr(E_b E_b' H))/2
    Hb = [[arb(0)] * ns for _ in range(ns)]
    for b in range(ns):
        for bp in range(ns):
            acc = arb(0)
            for (a1, a2, v1) in E[b]:
                for (c1, c2, v2) in E[bp]:
                    # tr(E_b H E_b'): sum E_b[a1,a2] H[a2,c1] E_b'[c1,c2] [c2 == a1]
                    if c2 == a1:
                        acc += v1 * H[a2][c1] * v2
                    # tr(E_b E_b' H): sum E_b[a1,a2] E_b'[a2,c2] H[c2,a1]  (requires c1 == a2)
                    if c1 == a2:
                        acc += v1 * v2 * H[c2][a1]
            Hb[b][bp] = acc / 2
    A3 = flattening(G, 3, 1)       # NV x ns
    B4 = flattening(G, 4, 2)       # ns x ns
    n3, n4 = res["Gn"][3], res["Gn"][4]
    # higher degrees: flattening norms (m >= 5), never worse than Frobenius (take the smaller certified bound)
    hn = {}
    for m in range(5, N + 1):
        fnm = flat_norm(G, m, 2)
        hn[m] = fnm if fnm < res["Gn"][m] else res["Gn"][m]
        log(f"    sup|G_{m}| <= {hn[m].str(6)}  (flattening (2,{m - 2}); Frobenius {res['Gn'][m].str(6)})")
    fl = lambda v: float(v.mid())
    A3f = np.array([[fl(v) for v in row] for row in A3])
    B4f = np.array([[fl(v) for v in row] for row in B4])
    Hbf = np.array([[fl(v) for v in row] for row in Hb])
    log(f"    float: lambda_max(Hb) = {np.linalg.eigvalsh(Hbf).max():.6f} (= max eigenvalue of H)")
    # non-rigorous self-test of the three identities at random unit vectors (guards against indexing mistakes)
    rng = np.random.default_rng(3)
    Gf = {k: fl(v) for k, v in G.items()}

    def pev(u, m):
        return sum(v * np.prod([u[j] ** e for j, e in enumerate(exps_of(k))]) for k, v in Gf.items() if deg(k) == m)
    errs = []
    for _ in range(5):
        u = rng.normal(size=NV)
        u /= np.linalg.norm(u)
        y = np.array([np.sqrt(2 / _fact_t(be)) * np.prod([u[j] ** e for j, e in enumerate(be)]) for be in B2])
        errs += [abs(y @ y - 1), abs(y @ Hbf @ y - pev(u, 2)), abs(u @ A3f @ y - pev(u, 3)), abs(y @ B4f @ y - pev(u, 4))]
    log(f"    float self-test of |y| = 1, G_2 = y'Hb y, G_3 = u'A3 y, G_4 = y'B4 y: max error {max(errs):.1e}")
    assert max(errs) < 1e-10

    def Mfloat(s, t):
        return np.block([[t * np.eye(NV), s * A3f / 2], [s * A3f.T / 2, Hbf + s * s * B4f - t * np.eye(ns)]])

    worst = None
    Tf = float(Tq)
    for kk in range(K):
        sm_q = Fraction(2 * kk + 1, 2 * K) * Tq          # interval midpoint (exact rational)
        rho = arb(fmpq(Tq.numerator, Tq.denominator * 2 * K))
        sm = arb(fmpq(sm_q.numerator, sm_q.denominator))
        s_hi = sm + rho
        smf = float(sm_q)
        # choose t (float) minimising lambda_max
        ts = np.linspace(-0.2, 0.2, 81)
        vals = [np.linalg.eigvalsh(Mfloat(smf, t)).max() for t in ts]
        t0_ = ts[int(np.argmin(vals))]
        lo, hi = t0_ - 0.005, t0_ + 0.005
        for _ in range(40):                              # golden-section refinement
            m1, m2 = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
            if np.linalg.eigvalsh(Mfloat(smf, m1)).max() < np.linalg.eigvalsh(Mfloat(smf, m2)).max():
                hi = m2
            else:
                lo = m1
        t_k = float(np.round(((lo + hi) / 2) * 2 ** 30) / 2 ** 30)
        lf = np.linalg.eigvalsh(Mfloat(smf, t_k)).max()
        ta = arb(t_k)
        Mx = [[arb(0)] * (NV + ns) for _ in range(NV + ns)]
        for i in range(NV):
            Mx[i][i] = ta
            for j in range(ns):
                v = sm * A3[i][j] / 2
                Mx[i][NV + j] = v
                Mx[NV + j][i] = v
        s2 = sm * sm
        for i in range(ns):
            for j in range(ns):
                Mx[NV + i][NV + j] = Hb[i][j] + s2 * B4[i][j] - (ta if i == j else 0)
        mu = cert_lmax(Mx, _up_guess(lf))
        assert mu is not None, "Cholesky certification of lambda_max(M) failed"
        mu_int = mu + rho * (n3 / 2 + 2 * sm * n4) + rho * rho * n4   # covers the whole interval
        hi_terms = arb(0)
        for m in range(5, N + 1):
            p = hn[m]
            for _ in range(m - 2):
                p = p * s_hi
            hi_terms += p
        bound_k = (2 * mu_int + hi_terms).upper()
        if worst is None or bound_k > worst:
            worst, kw = bound_k, kk
        if kk in (0, K // 2, K - 1):
            log(f"    interval {kk}: s in [{float(sm_q) - Tf / (2 * K):.5f}, {float(sm_q) + Tf / (2 * K):.5f}], t = {t_k:.5f}, "
                f"2 mu = {(2 * mu).str(6)}, + width/higher terms -> {bound_k.str(6)}")
    total = -(worst + res["err_term"])
    log(f"    worst interval {kw}: sup_(|u|=1) sum_m G_m(u) s^(m-2) <= {worst.str(8)}; + errF/T^2 = {res['err_term'].str(6)}")
    ok = total > 0
    log(f"[7b] RESULT (joint): F(x) <= -({total.lower().str(8)}) |x|^2 for all |x| <= T : {ok}   ({time.time() - tj:.1f} s)")
    kap = (total / 2).lower()
    log(f"    => on E with |x| <= T:  L - log({C}) <= Phi/2 <= -kappa |x|^2 with kappa = {kap.str(8)}")
    return dict(ok=bool(ok), total=total, kappa=kap)


def design_file(Tq, N):
    if D == 6:
        return f"design_T{Tq.numerator}_{Tq.denominator}.json"
    return f"design_d{D}_T{Tq.numerator}_{Tq.denominator}_N{N}.json"


OPTS = {"joint": False, "K": 64}


def parse_args(argv):
    """Returns (d, redesign, positional args); the flags --d=, --redesign, --joint, --K= may appear anywhere
    (--joint and --K= are stored in OPTS)."""
    d, redesign, pos = 6, False, []
    for a in argv:
        if a.startswith("--d="):
            d = int(a[4:])
        elif a == "--joint":
            OPTS["joint"] = True
        elif a.startswith("--K="):
            OPTS["K"] = int(a[4:])
        elif a == "--redesign":
            redesign = True
        else:
            pos.append(a)
    return d, redesign, pos


def main():
    d, redesign, argv = parse_args(sys.argv[1:])
    setup(d)
    Tq = Fraction(argv[0]) if len(argv) > 0 else Fraction(1, 19)
    N = int(argv[1]) if len(argv) > 1 else 6
    gamma_arg = float(argv[2]) if len(argv) > 2 else None
    dfile = design_file(Tq, N)
    if os.path.exists(dfile) and not redesign:
        res = certify(Tq, N, gamma_arg, design_in=dfile)
    else:
        res = certify(Tq, N, gamma_arg, design_out=dfile)
    if OPTS["joint"]:
        jr = joint_certificate(res, OPTS["K"])
        if jr["ok"] and (not res["ok"] or jr["kappa"] > res["kappa"]):
            res = dict(res, ok=True, total=jr["total"], kappa=jr["kappa"])
            log("    (the joint bound is used for the final statement)")
    # a clean rational kappa below the certified enclosure
    kq = Fraction(int(float(res["kappa"].mid()) * 10 ** 4), 10 ** 4)
    kq_ok = res["kappa"] > arb(fmpq(kq.numerator, kq.denominator))
    log(f"[8] FINAL: for all x in E with |x|_2 <= {Tq}:  L(x) <= log({C}) - {kq} |x|_2^2   (certified: {kq_ok and res['ok']})")
    log(f"    total time {time.time() - t0:.1f} s")
    return res["ok"] and kq_ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
