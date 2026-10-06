"""models_v3_check.py — check every interval bound and affine model printed by the `models` mode of a smale_bb
binary against critical values computed from first principles (firstprinciples/fpcore.py: P built from its critical
points, P' multiplied out and integrated; no integral formula), at random points and corners of random domain boxes.
usage: python3 c/regress/models_v3_check.py BIN D NBOX [SEED]
Boxes: half-widths 2^-k (k = 2..12, per coordinate), centres on the matching dyadic grid, |C_k| + H_k <= 1, i.e. the
evaluator's domain (A5).  Half of the boxes are drawn near the d = 7 region where the w-model matters (some |u_j|
between 0.25 and 0.6).  Float evaluation, tolerance 1e-9 (1 + |value|); this is a test, not a proof."""
import os, sys, subprocess
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "firstprinciples"))
from fpcore import crit_values, x_to_u

BIN, d, NBOX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
n = d - 1; D = 2 * (n - 1)
boxes = []
while len(boxes) < NBOX:
    k = rng.integers(2, 13, size=D); H = 2.0 ** (-k.astype(float))
    if len(boxes) % 2:
        # near |u_v| in [0.25, 0.6] for one random v, the others anywhere in the unit disc
        u = (rng.uniform(-1, 1, n - 1) + 1j * rng.uniform(-1, 1, n - 1))
        v = rng.integers(n - 1); u[v] = rng.uniform(0.25, 0.6) * np.exp(1j * rng.uniform(0, 2 * np.pi))
        X = np.empty(D); X[0::2] = u.real; X[1::2] = u.imag
    else:
        X = rng.uniform(-1, 1, D)
    C = (np.floor(X / (2 * H)) * 2 + 1) * H              # odd multiples of H: a cell of the dyadic grid
    if np.all(np.abs(C) + H <= 1):
        boxes.append((C, H))
inp = "".join(" ".join(float.hex(float(c)) for c in C) + " " + " ".join(float.hex(float(h)) for h in H) + "\n" for C, H in boxes)
out = subprocess.run([BIN, "models", str(d)], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
lines = [l for l in out if l.strip()]
assert len(lines) == NBOX * n, (len(lines), NBOX * n)
nmod = nviol = niv = nvalid = 0
worst = 0.0
for b, (C, H) in enumerate(boxes):
    M = [list(map(float, lines[b * n + i].split())) for i in range(n)]
    pts = C + H * rng.uniform(-1, 1, (48, D))
    cor = C + H * (2 * ((np.arange(16)[:, None] >> np.arange(D)[None, :]) & 1) - 1)
    pts = np.vstack([pts, cor])
    u = x_to_u(pts)
    keep = np.all(np.abs(u) > 0, axis=1)
    pts = pts[keep]
    V = crit_values(x_to_u(pts))[0]
    LV = np.log(np.abs(V))
    for i in range(n):
        valid, a, e, lo, hi = M[i][:5]; G = np.array(M[i][5:5 + D])
        tol = 1e-9 * (1 + np.abs(LV[:, i]))
        niv += int(np.sum((LV[:, i] > hi + tol) | (LV[:, i] < lo - tol)))
        if valid:
            nvalid += 1
            err = np.abs(LV[:, i] - (a + (pts - C) @ G))
            nmod += len(err); nviol += int(np.sum(err > e + tol))
            if e > 0: worst = max(worst, float(np.max(err / e)))
print(f"{BIN.split('/')[-1]} d={d}: {NBOX} boxes, {nvalid} valid models, {nmod} affine checks, {nviol} model violations, "
      f"{niv} interval violations, max err/e = {worst:.3f}")
sys.exit(1 if nviol or niv else 0)
