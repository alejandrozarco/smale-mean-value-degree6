"""Rigorous local certificate near the regular configuration, d = 6 (see README.md for the mathematics).

Usage:  python3 certify.py [T] [N] [gamma]
    T      radius in the x-chart, a rational like 1/19 (default 1/19)
    N      Taylor truncation degree (default 6)
    gamma  design parameter (float) for the multipliers (default: automatic)

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
NV = 8          # real variables x_1..x_8
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
    return [acb(arb.cos_pi_fmpq(fmpq(2 * n, 5)), arb.sin_pi_fmpq(fmpq(2 * n, 5))) for n in range(5)]


W = omega_pows()


def qomega_to_acb(v):
    s = acb(0)
    for n in range(5):
        if v[n] != 0:
            s += W[n] * arb(fmpq(v[n].numerator, v[n].denominator))
    return s


def complex_series_to_real(cs):
    """cs: dict beta(4-tuple) -> acb coefficient of eps^beta, eps_k = x_{2k-1} + i x_{2k}.
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
    """Upper bound for sup_{|eps|=1} |p_m(eps)| for complex homogeneous part of degree m (4 complex vars)."""
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
    c = Fraction(5, 6)
    log(f"=== Local certificate, d=6: T = {Tq} (= {float(Tq):.6f}), Taylor degree N = {N}, prec = {ctx.prec} bits")

    # ---------------- exact facts
    log("[1] exact facts at x = 0 (exact arithmetic in Q(omega))")
    ok, _ = EQ.exact_facts(verbose=verbose)
    assert ok, "exact facts failed"

    # ---------------- Taylor polynomials of sigma_i = S_i / c
    log("[2] exact Taylor coefficients of S_i up to degree", N, "-> arb; real polynomials U_i + i V_i = S_i/c")
    sig_c = []      # complex coefficient dicts of sigma_i
    UV = []
    tau = []        # tail bounds tau_i(T) for |S_i/c - sigma_i| on |x| <= T
    for i in range(1, 6):
        tc = EQ.taylor_S(i, N)
        cs = {}
        for beta, v in tc.items():
            vv = [x / c for x in v]
            if sum(beta) == 0:
                assert EQ.qeq(vv, [1, 0, 0, 0, 0])
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
        tau.append((tsum / arb(fmpq(5, 6))).upper())
    log("    tail bounds tau_i(T):", [t.str(5) for t in tau])

    # complex Frobenius norms n_{i,m}
    nrm = [[cfrob(sig_c[i], m) if m > 0 else arb(1) for m in range(N + 1)] for i in range(5)]

    # ---------------- A_i = |S_i/c|^2 : explicit part alpha_i (deg <= Nt) + error E_i(T)
    log("[3] alpha_i = [U_i^2 + V_i^2]_{<=Nt} and error bounds E_i(T)")
    alpha, E = [], []
    for i in range(5):
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

    # Phi = (1/5) sum alpha_i - 1 ; D_i = alpha_i - alpha_1
    phi = {}
    for i in range(5):
        phi = padd(phi, pscale(alpha[i], arb(fmpq(1, 5))))
    phi[0] = phi.get(0, arb(0)) - 1
    # constant and linear terms of Phi are exactly 0 (facts F1, F2); check the enclosures contain 0 and set them
    for k in list(phi):
        if deg(k) <= 1:
            assert phi[k].contains(0), ("enclosure of an exactly-zero coefficient does not contain 0", k, phi[k])
            del phi[k]
    d = []
    for i in range(1, 5):
        di = padd(alpha[i], alpha[0], -1)
        assert di[0].contains(0)
        del di[0]           # D_i(0) = 0 exactly (F1)
        d.append(di)

    # ---------------- design (floating point, non-rigorous choices)
    log("[4] design of multipliers (floating point; choices only)")

    def fl(v):
        return float(v.mid())

    P = np.zeros((8, 8))
    for k, v in homog(phi, 2).items():
        ex = exps_of(k)
        idx = [j for j in range(8) for _ in range(ex[j])]
        if idx[0] == idx[1]:
            P[idx[0], idx[0]] = fl(v)
        else:
            P[idx[0], idx[1]] = P[idx[1], idx[0]] = fl(v) / 2
    Lam = np.zeros((4, 8))
    for a in range(4):
        for k, v in homog(d[a], 1).items():
            Lam[a, exps_of(k).index(1)] = fl(v)
    _, sv, Vt = np.linalg.svd(Lam)
    Kb, Kp = Vt[4:].T, Vt[:4].T
    M = Lam @ Kp
    Pzz, Pzw, Pww = Kb.T @ P @ Kb, Kb.T @ P @ Kp, Kp.T @ P @ Kp
    lamK = -np.linalg.eigvalsh(Pzz).max()
    log(f"    singular values of Lambda: {sv}")
    log(f"    Phi_2 on ker Lambda: max eigenvalue {-lamK:.6f}")
    gamma = gamma_arg if gamma_arg is not None else 1.1 * lamK
    Gz = -2 * np.linalg.inv(M).T @ Pzw.T
    Gw = -np.linalg.inv(M).T @ (gamma * np.eye(4) + Pww)
    Gam = Gz @ Kb.T + Gw @ Kp.T          # 4 x 8, mu_lin(x) = Gam x
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

    lin_keys = [unit_key(j) for j in range(8)]
    quad_keys = sorted({unit_key(a) + unit_key(b) for a in range(8) for b in range(a, 8)})
    cub_keys = sorted({q + unit_key(j) for q in quad_keys for j in range(8)})
    quart_keys = sorted({q + r for q in quad_keys for r in quad_keys})

    def weight(k):
        m = deg(k)
        af = 1
        for a in exps_of(k):
            af *= factorial(a)
        return np.sqrt(af / factorial(m))

    fphi = fpoly(phi)
    fd = [fpoly(x) for x in d]
    gam_polys = [{lin_keys[j]: Gam[a, j] for j in range(8)} for a in range(4)]
    base = dict(fphi)
    for a in range(4):
        for k, v in fmul(gam_polys[a], fd[a], Nt).items():
            base[k] = base.get(k, 0.0) + v
    ci = {k: n for n, k in enumerate(cub_keys)}
    qi = {k: n for n, k in enumerate(quart_keys)}
    nunk = 4 * len(quad_keys)
    A3 = np.zeros((len(cub_keys), nunk))
    A4 = np.zeros((len(quart_keys), nunk))
    for a in range(4):
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
        Gam = np.array([[float.fromhex(v) for v in row] for row in dz["Gamma"]])
        psi = np.array([float.fromhex(dz["Psi"][a][str(exps_of(qk))]) for a in range(4) for qk in quad_keys])
        gamma = dz["gamma"]
        log(f"    design loaded from {design_in} (Gamma, Psi as saved; the float design above is not used)")

    # ---------------- rigorous assembly of F = Phi + sum_a mu_a D_a
    log("[5] rigorous assembly of the explicit part G of F = Phi + sum_a mu_a(x) D_a(x)")
    Gam_a = [{lin_keys[j]: arb(float(Gam[a, j])) for j in range(8) if Gam[a, j] != 0} for a in range(4)]
    Psi_a = [{qk: arb(float(psi[a * len(quad_keys) + qn])) for qn, qk in enumerate(quad_keys)
              if psi[a * len(quad_keys) + qn] != 0} for a in range(4)]
    G = dict(phi)
    for a in range(4):
        mu = padd(Gam_a[a], Psi_a[a])
        G = padd(G, pmul(mu, d[a], Nt))
    # G has no terms of degree 0 or 1: phi has none (exact), mu*d starts at degree 2.
    assert all(deg(k) >= 2 for k in G)

    # ---------------- G_2 <= -lambda |x|^2 via an arb Cholesky of -H2 - lambda I
    H2 = [[arb(0)] * 8 for _ in range(8)]
    for k, v in homog(G, 2).items():
        ex = exps_of(k)
        idx = [j for j in range(8) for _ in range(ex[j])]
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
    Mx = [[(-H2[i][j] - (lam_a if i == j else 0)) for j in range(8)] for i in range(8)]
    Lc = [[arb(0)] * 8 for _ in range(8)]
    pd = True
    for j in range(8):
        s = Mx[j][j]
        for k in range(j):
            s -= Lc[j][k] * Lc[j][k]
        if not (s > 0):
            pd = False
            break
        Lc[j][j] = s.sqrt()
        for i in range(j + 1, 8):
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
    for i in range(5):
        errF += E[i] / 5
    for a in range(4):
        ng = _sqrt_up(sum((arb(float(Gam[a, j])) * arb(float(Gam[a, j])) for j in range(8)), arb(0)))
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
        f"L - log(5/6) <= Phi/2 <= -kappa |x|^2 with kappa = {kap.str(8)}")
    res = dict(T=Tq, N=N, ok=bool(ok), lam=lam, total=total, kappa=kap, budget=budget, err_term=err_term,
               Gn=Gn, G=G, phi=phi, d=d, Gam=Gam, psi=psi, quad_keys=quad_keys, E=E, tau=tau, errF=errF,
               gamma=gamma)
    if design_out:
        import json
        with open(design_out, "w") as fh:
            json.dump({"T": str(Tq), "N": N, "gamma": gamma, "lambda": str(lam),
                       "Gamma": [[float(v).hex() for v in row] for row in Gam],
                       "Psi": [{str(exps_of(qk)): float(psi[a * len(quad_keys) + qn]).hex()
                                for qn, qk in enumerate(quad_keys)} for a in range(4)]}, fh, indent=1)
    return res


def main():
    Tq = Fraction(sys.argv[1]) if len(sys.argv) > 1 else Fraction(1, 19)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    gamma_arg = float(sys.argv[3]) if len(sys.argv) > 3 else None
    dfile = f"design_T{Tq.numerator}_{Tq.denominator}.json"
    if os.path.exists(dfile):
        res = certify(Tq, N, gamma_arg, design_in=dfile)
    else:
        res = certify(Tq, N, gamma_arg, design_out=dfile)
    # a clean rational kappa below the certified enclosure
    kq = Fraction(int(float(res["kappa"].mid()) * 10 ** 4), 10 ** 4)
    kq_ok = res["kappa"] > arb(fmpq(kq.numerator, kq.denominator))
    log(f"[8] FINAL: for all x in E with |x|_2 <= {Tq}:  L(x) <= log(5/6) - {kq} |x|_2^2   (certified: {kq_ok and res['ok']})")
    log(f"    total time {time.time() - t0:.1f} s")
    return res["ok"] and kq_ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
