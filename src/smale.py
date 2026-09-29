"""Core numerics for Smale's mean value conjecture (normalized form).

P'(z) = prod_{j=1}^{n} (1 - z/b_j),  n = d-1,   P(0)=0.
S_i   = P(b_i)/b_i = int_0^1 prod_j (1 - t b_i/b_j) dt      (integral form, no cancellation)
F(b)  = min_i |S_i|.     Conjecture: F <= n/d = (d-1)/d.

S_i depends only on ratios b_i/b_j: invariant under b -> lambda*b (lambda in C*),
under permutations, and F under complex conjugation.
"""
import numpy as np

_GL = {m: np.polynomial.legendre.leggauss(m) for m in range(1, 12)}


def gl01(m):
    x, w = _GL[m]
    return (x + 1) / 2, w / 2


def S_all(b):
    """S_i for all i, float/complex numerics (Gauss-Legendre, exact for poly degree n)."""
    b = np.asarray(b, dtype=complex)
    n = len(b)
    t, w = gl01(n // 2 + 2)
    R = b[:, None] / b[None, :]                       # R[i,j] = b_i/b_j
    fac = 1 - t[:, None, None] * R[None, :, :]        # (m, n, n)
    prod = np.prod(fac, axis=2)                       # (m, n)
    return w @ prod


def F(b):
    return np.min(np.abs(S_all(b)))


def equality_config(d):
    n = d - 1
    return np.exp(2j * np.pi * np.arange(n) / n)


def chart_to_b(p, n):
    """p = (x_2..x_n, y_2..y_n) log-polar coordinates, b_1 = 1."""
    x = p[: n - 1]
    y = p[n - 1:]
    return np.concatenate([[1.0 + 0j], np.exp(x + 1j * y)])


def b_to_chart(b):
    b = np.asarray(b, complex) / b[0]
    lb = np.log(b[1:])
    return np.concatenate([lb.real, lb.imag])


def canonical(b):
    """Canonical description modulo scale/rotation/permutation/conjugation:
    (sorted |S_i|, sorted log|b_i| relative to mean)."""
    S = np.abs(S_all(b))
    lm = np.log(np.abs(b))
    lm = np.sort(lm - lm.mean())
    return np.sort(S), lm
