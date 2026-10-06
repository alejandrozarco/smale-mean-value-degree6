"""fpcore.py -- first-principles evaluation of the critical-value quotients.

Given critical points b_1..b_n (all nonzero) of P with P(0)=0, P'(0)=1:
    P'(z) = prod_j (1 - z/b_j),   P(z) = int_0^z P',   V_i = P(b_i)/b_i.
We set b_1 = 1 and b_j = 1/u_j.  With v_j := 1/b_j (v_1 = 1, v_j = u_j):
    P'(z) = sum_k q_k z^k        (q from multiplying the linear factors 1 - v_j z)
    P(z)  = sum_k q_k z^{k+1}/(k+1)        (term-wise integration)
    V_i   = P(b_i)/b_i = sum_k q_k/(k+1) * b_i^k      (Horner in b_i = 1/v_i)
No integral formula over t is used here.
"""
import numpy as np
import mpmath as mp

EPS = 2.0 ** -53


def x_to_u(X):
    """X: (..., D) real -> (..., n-1) complex u_2..u_n."""
    return X[..., 0::2] + 1j * X[..., 1::2]


def crit_values(U):
    """U: (N, n-1) complex, all entries nonzero.  Returns V (N, n) complex and
    A (N, n): A_i = sum_k e_k(|v|)/(k+1) |b_i|^k, an upper bound for the sum of
    the moduli of all terms, used for a floating-point error estimate
    |fl(V_i) - V_i| <~ c_n * EPS * A_i."""
    N, m = U.shape
    n = m + 1
    v = np.concatenate([np.ones((N, 1), complex), U], axis=1)
    q = np.zeros((N, n + 1), complex); q[:, 0] = 1
    qa = np.zeros((N, n + 1)); qa[:, 0] = 1
    av = np.abs(v)
    for j in range(n):
        qn = q.copy(); qn[:, 1:] -= v[:, j:j + 1] * q[:, :-1]; q = qn
        qan = qa.copy(); qan[:, 1:] += av[:, j:j + 1] * qa[:, :-1]; qa = qan
    kk = np.arange(n + 1) + 1.0
    coef = q / kk
    coefa = qa / kk
    b = 1.0 / v
    ab = np.abs(b)
    V = np.repeat(coef[:, n:n + 1], n, axis=1)
    A = np.repeat(coefa[:, n:n + 1], n, axis=1)
    for k in range(n - 1, -1, -1):
        V = V * b + coef[:, k:k + 1]
        A = A * ab + coefa[:, k:k + 1]
    return V, A


def crit_values_mp(u, dps=50):
    """u: sequence of complex (u_2..u_n), evaluated in mpmath at dps digits."""
    with mp.workdps(dps):
        v = [mp.mpc(1)] + [mp.mpc(complex(z).real, complex(z).imag) if not isinstance(z, mp.mpc) else z for z in u]
        n = len(v)
        q = [mp.mpc(1)] + [mp.mpc(0)] * n
        for j in range(n):
            qn = q[:]
            for k in range(1, n + 1):
                qn[k] = q[k] - v[j] * q[k - 1]
            q = qn
        coef = [q[k] / (k + 1) for k in range(n + 1)]
        out = []
        for i in range(n):
            b = 1 / v[i]
            acc = coef[n]
            for k in range(n - 1, -1, -1):
                acc = acc * b + coef[k]
            out.append(acc)
        return out


def crit_values_from_roots_np(b):
    """Independent reference: build P' from its roots with numpy.poly, integrate
    with numpy polyint, evaluate P(b_i)/b_i.  b: 1-D complex array of critical
    points.  Used only to cross-check crit_values on random inputs."""
    b = np.asarray(b, complex)
    c = np.poly(b)                       # monic, highest first: prod (z - b_j)
    c = c / c[-1]                        # normalise so that P'(0) = 1
    P = np.polyint(c)                    # P(0) = 0
    return np.polyval(P, b) / b
