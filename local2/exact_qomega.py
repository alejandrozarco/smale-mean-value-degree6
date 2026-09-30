"""Exact arithmetic in Q(omega), omega = exp(2*pi*i/5), and the exact Taylor coefficients of S_i.

Representation.  An element of the group ring Q[C_5] is a list v[0..4] of Fractions meaning sum_n v[n] omega^n.
Two such vectors define the same element of Q(omega) iff their difference is a constant vector
(the kernel of Q[C_5] -> Q(omega) is spanned by 1 + omega + ... + omega^4).  normalize() removes the
omega^4 coordinate, which gives the unique coordinates in the basis 1, omega, omega^2, omega^3.

Exponential-sum form of S_i.  With b_j = omega^(j-1) exp(eps_j) (eps_1 = 0) and r_ij = b_i / b_j,
    prod_j (1 - t r_ij) = sum_k (-1)^k e_k(r_i1..r_i5) t^k,   int_0^1 t^k dt = 1/(k+1),
so  S_i = sum_{J subset {1..5}} (-1)^|J| / (|J|+1) * prod_{j in J} r_ij
        = sum_J  coef_J * omega^(e_J) * exp(lambda_J . eps),
with e_J = sum_{j in J} (i - j) mod 5 and lambda_J in Z^4 (coordinates eps_2..eps_5):
    lambda_J[k] = |J| [k = i] - [k in J].
Hence the Taylor coefficient of eps^beta in S_i is  sum_J coef_J omega^(e_J) lambda_J^beta / beta!,
an element of Q(omega) that is computed here exactly.
"""
from fractions import Fraction
from itertools import combinations
from math import factorial

N_PTS = 5


def normalize(v):
    """Coordinates in the basis 1, w, w^2, w^3 (w^4 = -1 - w - w^2 - w^3)."""
    return tuple(Fraction(v[n]) - Fraction(v[4]) for n in range(4))


def qeq(u, v):
    return normalize(u) == normalize(v)


def qzero(v):
    return all(c == 0 for c in normalize(v))


def subsets_data(i):
    """List of (coef, e, lam) for S_i, i in 1..5.  coef is a Fraction, e an int mod 5, lam a 4-tuple of ints."""
    out = []
    for size in range(0, 6):
        for J in combinations(range(1, 6), size):
            coef = Fraction((-1) ** size, size + 1)
            e = sum(i - j for j in J) % 5
            lam = tuple(size * (1 if k == i else 0) - (1 if k in J else 0) for k in range(2, 6))
            out.append((coef, e, lam))
    return out


def monomials4(N):
    """All exponent 4-tuples of total degree <= N, ordered by degree."""
    out = []
    for d in range(N + 1):
        def rec(prefix, left, nv):
            if nv == 1:
                out.append(tuple(prefix + [left]))
                return
            for a in range(left, -1, -1):
                rec(prefix + [a], left - a, nv - 1)
        rec([], d, 4)
    return out


def taylor_S(i, N):
    """Exact Taylor coefficients of S_i in eps (eps_2..eps_5) up to total degree N.
    Returns dict beta -> group-ring vector (list of 5 Fractions)."""
    data = subsets_data(i)
    res = {}
    for beta in monomials4(N):
        bf = 1
        for b in beta:
            bf *= factorial(b)
        v = [Fraction(0)] * 5
        for coef, e, lam in data:
            p = 1
            for lk, bk in zip(lam, beta):
                p *= lk ** bk
            if p:
                v[e] += coef * p
        res[beta] = [x / bf for x in v]
    return res


def exact_facts(verbose=True):
    """Prove (exactly) the facts at x = 0 used by the certificate:
    (F1) S_i(0) = 5/6 for i = 1..5;
    (F2) for each k = 2..5, sum_i dS_i/d eps_k (0) = 0 in Q(omega).
    Returns the exact first-order coefficients as well (for the record)."""
    ok = True
    lin = {}
    for i in range(1, 6):
        tc = taylor_S(i, 1)
        s0 = tc[(0, 0, 0, 0)]
        f1 = qeq(s0, [Fraction(5, 6), 0, 0, 0, 0])
        ok &= f1
        for beta, v in tc.items():
            if sum(beta) == 1:
                lin[(i, beta.index(1) + 2)] = normalize(v)
        if verbose:
            print(f"  S_{i}(0) = 5/6 exactly: {f1}")
    for k in range(2, 6):
        tot = [Fraction(0)] * 4
        for i in range(1, 6):
            tot = [a + b for a, b in zip(tot, lin[(i, k)])]
        z = all(t == 0 for t in tot)
        ok &= z
        if verbose:
            print(f"  sum_i dS_i/deps_{k}(0) = 0 exactly in Q(omega): {z}")
    return ok, lin


if __name__ == "__main__":
    ok, lin = exact_facts()
    for key in sorted(lin):
        print(key, lin[key])
    print("ALL EXACT FACTS:", ok)
