"""Exact check (sympy, Gaussian rationals) that the objective map S: (b_2..b_n) -> (S_1..S_n), b_1 = 1,
has complex rank n-1 = d-2 at a sample point, so that the Zariski closure of its image is a hypersurface
(used in REDUCTION.md, Step 1)."""
import sympy as sp, sys
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = d - 1
t = sp.Symbol('t')
bs = sp.symbols('b2:%d' % (n + 1))
b = [sp.Integer(1)] + list(bs)
S = [sp.integrate(sp.expand(sp.prod([1 - t * b[i] / b[j] for j in range(n)])), (t, 0, 1)) for i in range(n)]
pt = {bs[k]: v for k, v in enumerate([2, 3 + sp.I, -1 + 2 * sp.I, sp.Rational(1, 2) - sp.I, 5, -3 + sp.I][: n - 1])}
J = sp.Matrix([[sp.nsimplify(sp.simplify(sp.diff(Si, x).subs(pt))) for x in bs] for Si in S])
print("d =", d, " Jacobian shape", J.shape, " exact rank =", J.rank(simplify=True))
