"""eqpoints_v2.py — Euclidean accuracy of the equality points used by smale_bb_v2's exclusion test (mode 2), in both
directions: every computed point is close to an exact equality point and every exact point has a close computed point.
Exact points: u_j = conj(omega^{a_j}), omega = e^{2 pi i/n}, (a_2..a_n) a permutation of (1..n-1).  mpmath, 300 bits.
Certified part: each computed point is matched to its exact point and ||p_hat - p||^2 < (1.1e-16)^2 is proved in arb
(python-flint, 300 bits; cos/sin via arb.cos_pi/sin_pi, rigorous ball arithmetic; the doubles enter exactly).
usage: python3 c/regress/eqpoints_v2.py BIN   (asserts max Euclidean error <= 1.1e-16 and count = (n-1)!)"""
import sys, itertools, math, subprocess
import mpmath as mp
import flint
flint.ctx.prec = 300
BOUND2 = flint.arb(121) / flint.arb(10) ** 34   # (1.1e-16)^2, exact rational in a tight ball
mp.mp.prec = 300
BIN = sys.argv[1]
ok = True
for d in (4, 5, 6, 7):
    n = d - 1
    out = subprocess.run([BIN, 'roots', str(d)], capture_output=True, text=True, check=True).stdout
    got = [[mp.mpf(float.fromhex(x)) for x in l.split()] for l in out.strip().split('\n')]
    ex = []
    for p in itertools.permutations(range(1, n)):
        v = []
        for a in p:
            z = mp.conj(mp.expjpi(mp.mpf(2 * a) / n)); v += [z.real, z.imag]
        ex.append(v)
    dist = lambda x, y: mp.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)))
    fwd = max(min(dist(g, e) for e in ex) for g in got)
    bwd = max(min(dist(g, e) for g in got) for e in ex)
    cnt = len(got) == math.factorial(n - 1)
    # certified: every computed point (as exact doubles) is within 1.1e-16 of the exact point it approximates
    exa = [[(flint.arb(2 * a) / n).cos_pi() if k == 0 else -(flint.arb(2 * a) / n).sin_pi()
            for a in p for k in (0, 1)] for p in itertools.permutations(range(1, n))]
    raw = [[float.fromhex(x) for x in l.split()] for l in out.strip().split('\n')]
    cert = True
    for g in raw:
        best = min(range(len(ex)), key=lambda j: dist([mp.mpf(x) for x in g], ex[j]))
        diffs = [flint.arb(x) - y for x, y in zip(g, exa[best])]
        d2 = sum(t * t for t in diffs)          # not t ** 2: python-flint 0.6 returns nan for (ball containing 0) ** 2
        cert &= bool(d2.is_finite() and d2 < BOUND2)
    good = cnt and fwd <= 1.1e-16 and bwd <= 1.1e-16 and cert
    ok &= good
    print(f'd={d}: {len(got)} points (expected {math.factorial(n - 1)}), max dist computed->exact {mp.nstr(fwd, 4)}, '
          f'exact->computed {mp.nstr(bwd, 4)}, arb-certified < 1.1e-16: {cert}  {"PASS" if good else "FAIL"}')
print('ALL PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
