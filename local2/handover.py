"""Hand-over from the global computation (SPEC section 3), verified exactly / in arb.

(H1) Relabelling.  For a permutation (a_2..a_5) of (1,2,3,4) let tau(1) = 1 and tau(k) = the index j with
     a_j = k - 1 (k = 2..5).  Put b'_k = b_{tau(k)}.  Then u'_k = u_{tau(k)} and
         sum_k |u'_k - conj(omega)^(k-1)|^2 = sum_j |u_j - conj(omega)^(a_j)|^2,
     so the ball around p = (conj(w)^a_2, ..., conj(w)^a_5) is mapped onto the ball around
     p0 = (conj(w), conj(w)^2, conj(w)^3, conj(w)^4) (the inverse relabelling maps it back).
     S_k(b') = int prod_j (1 - t b_tau(k)/b_tau(j)) dt = S_tau(k)(b), since j -> tau(j) is a bijection of {1..5}:
     the S_i are permuted, so E and min_i |S_i| are invariant.  Checked below for all 24 permutations
     (the index identities exactly; the S-identity is the displayed one-line argument).

(H2) Chart.  On the ball around p0 write w_j = u_j / p0_j - 1, |w_j| = |u_j - p0_j| <= r = 1/20.
     eps_j := -Log(1 + w_j) (principal branch) gives b_j = 1/u_j = omega^(j-1) exp(eps_j).
     |eps_j| <= sum_m |w_j|^m / m = -log(1 - |w_j|) <= C |w_j| with C = -log(1 - r)/r, because
     g(q) = -log(1-q)/q = sum_m q^(m-1)/m is increasing on [0, 1).  Hence
         |x|_2^2 = sum_j |eps_j|^2 <= C^2 sum_j |w_j|^2 <= C^2 r^2,   |x|_2 <= C r = log(20/19).
     Certified below: log(20/19) < 1/19.
"""
import itertools
from fractions import Fraction
from flint import arb, fmpq, ctx

ctx.prec = 200


def main(T=Fraction(1, 19)):
    ok = True
    print("[H1] relabelling of the 24 balls")
    count = 0
    for a in itertools.permutations([1, 2, 3, 4]):
        amap = {j: a[j - 2] for j in range(2, 6)}            # a_j
        tau = {1: 1}
        for k in range(2, 6):
            tau[k] = [j for j in range(2, 6) if amap[j] == k - 1][0]
        bij = sorted(tau.values()) == [1, 2, 3, 4, 5]
        centre = all(amap[tau[k]] == k - 1 for k in range(2, 6))   # u'_k is compared with conj(w)^(k-1)
        ok &= bij and centre
        count += bij and centre
    print(f"     permutations checked: {count}/24; each relabelling is a bijection fixing index 1 and sends the")
    print("     centre conj(w)^(a_j) of coordinate j to the centre conj(w)^(k-1) of coordinate k = tau^-1(j).")

    print("[H2] ball around p0 in the eps-chart")
    r = arb(fmpq(1, 20))
    C = -(1 - r).log() / r
    bound = C * r                     # = log(20/19)
    Ta = arb(fmpq(T.numerator, T.denominator))
    h2 = bound < Ta
    ok &= h2
    print(f"     C = -log(1-r)/r = {C.str(12)} ;  |x|_2 <= C r = log(20/19) = {bound.str(12)}")
    print(f"     log(20/19) < T = {T} = {Ta.str(12)} : {h2}")
    print("HANDOVER", "PASSED" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
