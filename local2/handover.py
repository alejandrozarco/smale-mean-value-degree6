"""Hand-over from the global computation (SPEC section 3), verified exactly / in arb.  Works for d in {5, 6, 7}.

Usage:  python3 handover.py [--d=D] [T] [r]        (defaults: d = 6, T = 1/19, r = 1/20)
        python3 handover.py --d=D --rmax T1 T2 ...  (largest exclusion radius supported by each certified T)

Notation: n = d - 1, omega = exp(2 pi i/n), b_1 = 1, u_j = 1/b_j (j = 2..n).

(H1) Relabelling.  For a permutation (a_2..a_n) of (1..n-1) let tau(1) = 1 and tau(k) = the index j with
     a_j = k - 1 (k = 2..n).  Put b'_k = b_{tau(k)}.  Then u'_k = u_{tau(k)} and
         sum_k |u'_k - conj(omega)^(k-1)|^2 = sum_j |u_j - conj(omega)^(a_j)|^2,
     so the ball around p = (conj(w)^a_2, ..., conj(w)^a_n) is mapped onto the ball around
     p0 = (conj(w), conj(w)^2, ..., conj(w)^(n-1)) (the inverse relabelling maps it back).
     S_k(b') = int prod_j (1 - t b_tau(k)/b_tau(j)) dt = S_tau(k)(b), since j -> tau(j) is a bijection of {1..n}:
     the S_i are permuted, so E and min_i |S_i| are invariant.  Checked below for all (n-1)! permutations
     (the index identities exactly; the S-identity is the displayed one-line argument, which does not use n).

(H2) Chart.  On the ball around p0 write w_j = u_j / p0_j - 1, |w_j| = |u_j - p0_j| (|p0_j| = 1), sum_j |w_j|^2 <= r^2.
     eps_j := -Log(1 + w_j) (principal branch) gives b_j = 1/u_j = omega^(j-1) exp(eps_j).
     |eps_j| <= sum_m |w_j|^m / m = -log(1 - |w_j|) <= C |w_j| with C = -log(1 - r)/r, because
     g(q) = -log(1-q)/q = sum_m q^(m-1)/m is increasing on [0, 1) and |w_j| <= r.  Hence
         |x|_2^2 = sum_j |eps_j|^2 <= C^2 sum_j |w_j|^2 <= C^2 r^2,   |x|_2 <= C r = -log(1 - r).
     This bound does not depend on n.  Certified below: -log(1 - r) < T.
     Conversely the largest r with -log(1 - r) <= T is r* = 1 - exp(-T); --rmax prints r* and certifies a
     rational r_q <= r* (4 significant digits, rounded down) with -log(1 - r_q) <= T.
"""
import itertools
import sys
from fractions import Fraction
from math import factorial
from flint import arb, fmpq, ctx

ctx.prec = 200


def check_relabelling(n):
    ok = True
    count = 0
    for a in itertools.permutations(range(1, n)):
        amap = {j: a[j - 2] for j in range(2, n + 1)}            # a_j
        tau = {1: 1}
        for k in range(2, n + 1):
            tau[k] = [j for j in range(2, n + 1) if amap[j] == k - 1][0]
        bij = sorted(tau.values()) == list(range(1, n + 1))
        centre = all(amap[tau[k]] == k - 1 for k in range(2, n + 1))   # u'_k is compared with conj(w)^(k-1)
        ok &= bij and centre
        count += bij and centre
    return ok, count, factorial(n - 1)


def main(T=Fraction(1, 19), d=6, r=Fraction(1, 20)):
    n = d - 1
    ok = True
    total = factorial(n - 1)
    print(f"[H1] relabelling of the {total} balls" + ("" if d == 6 else f"  (d = {d}, n = {n})"))
    okp, count, total = check_relabelling(n)
    ok &= okp
    print(f"     permutations checked: {count}/{total}; each relabelling is a bijection fixing index 1 and sends the")
    print("     centre conj(w)^(a_j) of coordinate j to the centre conj(w)^(k-1) of coordinate k = tau^-1(j).")

    print("[H2] ball around p0 in the eps-chart")
    ra = arb(fmpq(r.numerator, r.denominator))
    C = -(1 - ra).log() / ra
    bound = C * ra                     # = -log(1 - r)
    Ta = arb(fmpq(T.numerator, T.denominator))
    h2 = bound < Ta
    ok &= h2
    q = 1 - r
    print(f"     C = -log(1-r)/r = {C.str(12)} ;  |x|_2 <= C r = log({q.denominator}/{q.numerator}) = {bound.str(12)}")
    print(f"     log({q.denominator}/{q.numerator}) < T = {T} = {Ta.str(12)} : {h2}")
    print("HANDOVER", "PASSED" if ok else "FAILED")
    return ok


def rmax(Ts, d):
    n = d - 1
    okp, count, total = check_relabelling(n)
    print(f"[H1] d = {d}: relabelling checked for {count}/{total} permutations: {okp}")
    ok = okp
    print("[H2] largest exclusion radius r with -log(1 - r) <= T  (r* = 1 - exp(-T))")
    for T in Ts:
        Ta = arb(fmpq(T.numerator, T.denominator))
        rstar = 1 - (-Ta).exp()
        # 4 significant digits, rounded down, as an exact rational
        f = float(rstar.mid())
        e = 0
        while f * 10 ** e < 1000:
            e += 1
        rq = Fraction(int(f * 10 ** e), 10 ** e)
        rqa = arb(fmpq(rq.numerator, rq.denominator))
        good = (-(1 - rqa).log()) <= Ta
        ok &= good
        print(f"     T = {T}:  r* = {rstar.str(10)};  certified r = {rq} = {float(rq)}:  -log(1-r) = "
              f"{(-(1 - rqa).log()).str(10)} <= T : {good}")
    print("HANDOVER-RMAX", "PASSED" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    d, pos, want_rmax = 6, [], False
    for a in sys.argv[1:]:
        if a.startswith("--d="):
            d = int(a[4:])
        elif a == "--rmax":
            want_rmax = True
        else:
            pos.append(a)
    if want_rmax:
        sys.exit(0 if rmax([Fraction(p) for p in pos], d) else 1)
    T = Fraction(pos[0]) if pos else Fraction(1, 19)
    r = Fraction(pos[1]) if len(pos) > 1 else Fraction(1, 20)
    sys.exit(0 if main(T, d, r) else 1)
