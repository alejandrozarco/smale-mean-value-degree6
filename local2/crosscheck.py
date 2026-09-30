"""Non-rigorous consistency checks (floating point) of the certificate pipeline.

An independent float evaluator of S_i (expand prod_j (1 - t b_i/b_j) as a polynomial in t and integrate) is
compared with the explicit Taylor polynomials produced by certify.py:
  (a) |F_direct(x) - G(x)| against the rigorous remainder bound errF(|x|), at random x with |x| = s;
  (b) the ratio F_direct(x)/|x|^2 against the certified bound -(lambda - budget - err);
  (c) the relabelling symmetry used in the hand-over, S_k(b') = S_{tau(k)}(b), at random points.
None of this is part of the proof; it guards against coding mistakes.
"""
import sys
from fractions import Fraction
import numpy as np
import certify as C

w = np.exp(2j * np.pi / 5)


def S_direct(b):
    out = []
    for i in range(5):
        p = np.poly1d([1.0 + 0j])
        for j in range(5):
            p = p * np.poly1d([-b[i] / b[j], 1.0])
        co = p.c[::-1]
        out.append(sum(co[k] / (k + 1) for k in range(len(co))))
    return np.array(out)


def b_of_x(x):
    eps = np.concatenate([[0], x[0::2] + 1j * x[1::2]])
    return np.array([w ** j for j in range(5)]) * np.exp(eps)


def peval(p, x):
    s = 0.0
    for k, v in p.items():
        e = C.exps_of(k)
        s += float(v.mid()) * np.prod([x[j] ** e[j] for j in range(8)])
    return s


def main():
    Tq = Fraction(sys.argv[1]) if len(sys.argv) > 1 else Fraction(1, 19)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    res = C.certify(Tq, N, verbose=False)
    G, Gam, psi, qk = res["G"], res["Gam"], res["psi"], res["quad_keys"]
    T = float(Tq)
    errT = float(res["errF"].mid())
    rng = np.random.default_rng(12345)
    bound_ratio = -float(res["total"].mid())
    print(f"crosscheck at T = {Tq}, N = {N}: certified F/|x|^2 <= {bound_ratio:.6f}, errF(T) = {errT:.3e}")
    worst_rel, worst_ratio = 0.0, -1e9
    for s_frac in [1.0, 0.5, 0.25]:
        s = T * s_frac
        maxdiff = 0.0
        for trial in range(150):
            u = rng.normal(size=8)
            u /= np.linalg.norm(u)
            x = s * u
            A = np.abs(S_direct(b_of_x(x)) / (5 / 6)) ** 2
            Phi = A.mean() - 1
            D = A[1:] - A[0]
            mu = Gam @ x + np.array([sum(psi[a * len(qk) + n] * np.prod([x[j] ** C.exps_of(q)[j] for j in range(8)])
                                         for n, q in enumerate(qk)) for a in range(4)])
            F = Phi + mu @ D
            diff = abs(F - peval(G, x))
            maxdiff = max(maxdiff, diff)
            worst_ratio = max(worst_ratio, F / s ** 2)
        # errF(s) <= errF(T) (s/T)^(N+1)  (majorant series of order >= N+1)
        allowed = errT * s_frac ** (N + 1)
        print(f"  s = {s:.5f}: max |F_direct - G| = {maxdiff:.3e}   (rigorous bound errF(s) <= {allowed:.3e})")
        worst_rel = max(worst_rel, maxdiff / allowed)
    print(f"  worst sampled F/|x|^2 = {worst_ratio:.6f}  (certified <= {bound_ratio:.6f})")
    # relabelling check
    rng2 = np.random.default_rng(7)
    import itertools
    maxerr = 0.0
    for perm in itertools.permutations(range(1, 5)):
        b = np.concatenate([[1.0], rng2.normal(size=4) + 1j * rng2.normal(size=4)])
        bp = np.concatenate([[b[0]], b[1:][list(p - 1 for p in perm)]])
        S, Sp = S_direct(b), S_direct(bp)
        idx = [0] + list(perm)          # b'_k = b_{perm(k)}  =>  S_k(b') = S_{perm(k)}(b)
        maxerr = max(maxerr, np.max(np.abs(Sp - S[idx])))
    print(f"  relabelling: max |S_k(b') - S_tau(k)(b)| over 24 permutations = {maxerr:.2e}")
    ok = worst_rel < 1 and worst_ratio < bound_ratio and maxerr < 1e-10
    print("CROSSCHECK", "PASSED" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
