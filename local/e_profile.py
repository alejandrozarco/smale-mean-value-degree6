"""Numerical profile of log F on E near the equality point (d=6): for unit u in K = ker J and t > 0,
solve for x = t u + x_N (x_N in K^perp) with |S_1| = ... = |S_n| (Newton), report the deficit
(log c - L(x))/|x|^2 and |x_N|/|x|^2 (curvature of E).  Floating point (numpy), not rigorous: a numerical comparison
for the certified bound of local/local_cert.py.
usage: python3 local/e_profile.py [d=6] [OUT.tsv]   (OUT: one line per sample: t, deficit/|x|^2, |x_N|/|x|^2)"""
import sys, numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from smale import S_all
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = d - 1; nv = n - 1; D = 2 * nv
om = np.exp(2j * np.pi * np.arange(n) / n)
def ell(x):
    e = np.concatenate([[0], x[0::2] + 1j * x[1::2]])
    return np.log(np.abs(S_all(om * np.exp(e))))
def jac(f, x, h=1e-7):
    f0 = f(x); return np.array([(f(x + h * np.eye(D)[k]) - f(x - h * np.eye(D)[k])) / (2 * h) for k in range(D)]).T
J = jac(ell, np.zeros(D))
U_, s, Vt = np.linalg.svd(J)
KN = Vt[:n - 1].T; KK = Vt[n - 1:].T     # K-perp (4), K (4)
logc = np.log((d - 1) / d)
rng = np.random.default_rng(0)
def onE(u, t):
    y = np.zeros(n - 1)
    for it in range(50):
        x = t * u + KN @ y
        l = ell(x); G = (l - l.mean())[:n - 1]
        if np.abs(G).max() < 1e-14: break
        Jx = jac(lambda z: (ell(z) - ell(z).mean())[:n - 1], x)
        y -= np.linalg.solve(Jx @ KN, G)
    return x, ell(x)
out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else None
if out: out.write('# t\tdeficit_over_x2\tcurv\n')
for t in [0.005, 0.01, 0.02, 0.03, 0.04, 0.0527, 0.075, 0.1, 0.15, 0.2]:
    worst = np.inf; curv = 0; wd = None
    for k in range(400):
        u = KK @ rng.normal(size=D - n + 1); u /= np.linalg.norm(u)
        x, l = onE(u, t)
        r2 = x @ x
        dfc = (logc - l.mean()) / r2
        if dfc < worst: worst = dfc
        cv = np.linalg.norm(x - t * u) / r2
        curv = max(curv, cv)
        if out: out.write(f'{t}\t{dfc:.6e}\t{cv:.6e}\n')
    print(f"t={t}: min deficit (log c - log F)/|x|^2 on E = {worst:.5f};  max |x_N|/|x|^2 = {curv:.3f}", flush=True)
