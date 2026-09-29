"""Sanity check of the Taylor models in local_cert.py: at random eps in the complex polydisk ||eps||_inf <= T,
|log S_i(eps) - log c - P_i(eps)| must be <= r_i (and in fact <= r_i (||eps||_inf/T)^{N+1})."""
import sys, numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smale import S_all
from local_cert import build_h
d, T, N = 6, 0.0527, 6
n = d - 1
hs = build_h(d, N, T)
om = np.exp(2j * np.pi * np.arange(n) / n)
P = [{a: complex(float(c.real.mid()), float(c.imag.mid())) for a, c in h.p.items()} for h in hs]
r = [float(h.r.upper()) for h in hs]
rng = np.random.default_rng(3)
worst = 0; worst_s = 0
for _ in range(3000):
    e = T * np.sqrt(rng.uniform(size=n - 1)) * np.exp(2j * np.pi * rng.uniform(size=n - 1))
    s = np.abs(e).max()
    h = np.log(S_all(om * np.exp(np.concatenate([[0], e])))) - np.log(5 / 6)
    for i in range(n):
        p = sum(c * np.prod(e ** np.array(a)) for a, c in P[i].items())
        dev = abs(h[i] - p)
        worst = max(worst, dev / r[i]); worst_s = max(worst_s, dev / (r[i] * (s / T) ** (N + 1)))
print(f"max |h - P| / r = {worst:.3e};  max |h - P| / (r (s/T)^(N+1)) = {worst_s:.3e}  (both must be <= 1)")
