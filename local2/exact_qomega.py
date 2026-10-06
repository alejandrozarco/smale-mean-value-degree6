"""Exact arithmetic in Q(omega), omega = exp(2*pi*i/n), n = d - 1, and the exact Taylor coefficients of S_i.

The degree d is set with setup(d) (default d = 6, n = 5, as in SPEC.md); d in {5, 6, 7} is supported
(n = 4, 5, 6; Q(omega) = Q(i), Q(zeta_5), Q(sqrt(-3)) respectively).

Representation.  An element of the group ring Q[C_n] is a list v[0..n-1] of Fractions meaning sum_k v[k] omega^k.
Two such vectors define the same element of Q(omega) iff the polynomials sum_k v[k] x^k agree modulo the
cyclotomic polynomial Phi_n(x) (the minimal polynomial of omega).  normalize() reduces modulo Phi_n, which gives
the unique coordinates in the basis 1, omega, ..., omega^(phi(n)-1).  (For n = 5 this is the same as subtracting
the omega^4 coordinate, since omega^4 = -1 - omega - omega^2 - omega^3.)

Exponential-sum form of S_i.  With b_j = omega^(j-1) exp(eps_j) (eps_1 = 0) and r_ij = b_i / b_j,
    prod_j (1 - t r_ij) = sum_k (-1)^k e_k(r_i1..r_in) t^k,   int_0^1 t^k dt = 1/(k+1),
so  S_i = sum_{J subset {1..n}} (-1)^|J| / (|J|+1) * prod_{j in J} r_ij
        = sum_J  coef_J * omega^(e_J) * exp(lambda_J . eps),
with e_J = sum_{j in J} (i - j) mod n and lambda_J in Z^(n-1) (coordinates eps_2..eps_n):
    lambda_J[k] = |J| [k = i] - [k in J].
Hence the Taylor coefficient of eps^beta in S_i is  sum_J coef_J omega^(e_J) lambda_J^beta / beta!,
an element of Q(omega) that is computed here exactly.
"""
from fractions import Fraction
from itertools import combinations
from math import factorial

D = 6
N_PTS = 5
CYC = None          # coefficients of the cyclotomic polynomial Phi_n, low degree first (monic)


def _pdiv_exact(a, b):
    """Exact quotient of integer polynomials (low degree first), b monic, remainder must be 0."""
    a = list(a)
    q = [0] * (len(a) - len(b) + 1)
    for k in range(len(a) - len(b), -1, -1):
        q[k] = a[k + len(b) - 1]
        for j in range(len(b)):
            a[k + j] -= q[k] * b[j]
    assert all(c == 0 for c in a)
    return q


def cyclotomic(n):
    """Phi_n(x) = (x^n - 1) / prod_{m | n, m < n} Phi_m(x), exact integer coefficients, low degree first."""
    p = [-1] + [0] * (n - 1) + [1]
    for m in range(1, n):
        if n % m == 0:
            p = _pdiv_exact(p, cyclotomic(m))
    return p


def setup(d):
    global D, N_PTS, CYC
    assert d in (5, 6, 7), "supported degrees: d = 5, 6, 7"
    D, N_PTS = d, d - 1
    CYC = cyclotomic(N_PTS)


setup(6)


def c_value():
    """c = n/(n+1) = (d-1)/d."""
    return Fraction(N_PTS, N_PTS + 1)


def normalize(v):
    """Coordinates in the basis 1, w, ..., w^(phi(n)-1): reduce sum_k v[k] x^k modulo Phi_n (monic)."""
    a = [Fraction(x) for x in v]
    m = len(CYC) - 1
    for k in range(len(a) - 1, m - 1, -1):
        q = a[k]
        if q:
            for j in range(m + 1):
                a[k - m + j] -= q * CYC[j]
    return tuple(a[:m])


def qeq(u, v):
    return normalize(u) == normalize(v)


def qzero(v):
    return all(c == 0 for c in normalize(v))


def unit(x):
    """Group-ring vector of the rational x (times omega^0)."""
    return [Fraction(x)] + [Fraction(0)] * (N_PTS - 1)


def subsets_data(i):
    """List of (coef, e, lam) for S_i, i in 1..n.  coef is a Fraction, e an int mod n, lam an (n-1)-tuple of ints."""
    n = N_PTS
    out = []
    for size in range(0, n + 1):
        for J in combinations(range(1, n + 1), size):
            coef = Fraction((-1) ** size, size + 1)
            e = sum(i - j for j in J) % n
            lam = tuple(size * (1 if k == i else 0) - (1 if k in J else 0) for k in range(2, n + 1))
            out.append((coef, e, lam))
    return out


def monomials(N, nv=None):
    """All exponent nv-tuples (default nv = n - 1) of total degree <= N, ordered by degree."""
    nv = N_PTS - 1 if nv is None else nv
    out = []
    for d in range(N + 1):
        def rec(prefix, left, k):
            if k == 1:
                out.append(tuple(prefix + [left]))
                return
            for a in range(left, -1, -1):
                rec(prefix + [a], left - a, k - 1)
        rec([], d, nv)
    return out


monomials4 = monomials      # backward-compatible name (d = 6: 4 complex variables)


def taylor_S(i, N):
    """Exact Taylor coefficients of S_i in eps (eps_2..eps_n) up to total degree N.
    Returns dict beta -> group-ring vector (list of n Fractions)."""
    data = subsets_data(i)
    res = {}
    for beta in monomials(N):
        bf = 1
        for b in beta:
            bf *= factorial(b)
        v = [Fraction(0)] * N_PTS
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
    (F1) S_i(0) = c = (d-1)/d for i = 1..n;
    (F2) for each k = 2..n, sum_i dS_i/d eps_k (0) = 0 in Q(omega).
    Returns the exact first-order coefficients as well (for the record)."""
    n = N_PTS
    c = c_value()
    ok = True
    lin = {}
    zero = tuple([0] * n)
    for i in range(1, n + 1):
        tc = taylor_S(i, 1)
        s0 = tc[zero[:n - 1]]
        f1 = qeq(s0, unit(c))
        ok &= f1
        for beta, v in tc.items():
            if sum(beta) == 1:
                lin[(i, beta.index(1) + 2)] = normalize(v)
        if verbose:
            print(f"  S_{i}(0) = {c} exactly: {f1}")
    for k in range(2, n + 1):
        tot = [Fraction(0)] * len(lin[(1, 2)])
        for i in range(1, n + 1):
            tot = [a + b for a, b in zip(tot, lin[(i, k)])]
        z = all(t == 0 for t in tot)
        ok &= z
        if verbose:
            print(f"  sum_i dS_i/deps_{k}(0) = 0 exactly in Q(omega): {z}")
    return ok, lin


def fmt_sqrt_m3(a):
    """For n = 6: a = (a0, a1) means a0 + a1*omega with omega = (1 + sqrt(-3))/2; return 'p + q sqrt(-3)'."""
    a0, a1 = a
    return f"{a0 + a1 / 2} + ({a1 / 2}) sqrt(-3)"


if __name__ == "__main__":
    import sys
    dd = 6
    for arg in sys.argv[1:]:
        if arg.startswith("--d="):
            dd = int(arg[4:])
    setup(dd)
    print(f"d = {D}, n = {N_PTS}, Phi_n = {CYC} (low degree first)")
    ok, lin = exact_facts()
    for key in sorted(lin):
        if N_PTS == 6:
            print(key, lin[key], f"  = {fmt_sqrt_m3(lin[key])}")
        else:
            print(key, lin[key])
    print("ALL EXACT FACTS:", ok)
