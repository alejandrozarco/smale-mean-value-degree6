"""prim_check_v2.py: exact checks of prim_test_v2 output (stdin).
  mul/add/sub: Euclidean |exact(a op b) - computed mid| <= computed radius (Fractions).
  mod: lo <= |x+iy| <= hi and |x+iy| <= hi_fast (hi_fast only claimed for |x|,|y| <= 2^500), exact via squares.
  pow: x^e <= powu (exact).
  log: cr_log(x) is the correctly rounded log(x) (mpmath, 200 bits)."""
import sys
from fractions import Fraction as Q
import mpmath as mp
mp.mp.prec = 200
viol = {}; cnt = {}; worst = 0.0; logmax = 0.0
def bad(k, line):
    viol[k] = viol.get(k, 0) + 1
    if viol[k] <= 5: print('VIOL', k, line.strip())
it = iter(sys.stdin)
for line in it:
    v = line.split(); op = v[0]
    cnt[op] = cnt.get(op, 0) + 1
    if op in ('mul', 'add', 'sub'):
        x = [Q(float.fromhex(t)) for t in v[1:]]
        rad = Q(float.fromhex(next(it).split()[1]))
        ar, ai, br, bi, rr, ri = x
        if op == 'mul': er, ei = ar * br - ai * bi, ar * bi + ai * br
        elif op == 'add': er, ei = ar + br, ai + bi
        else: er, ei = ar - br, ai - bi
        e2 = (er - rr) ** 2 + (ei - ri) ** 2
        if e2 > rad * rad: bad(op, line)
        elif rad > 0 and e2 > 0: worst = max(worst, float(e2 / (rad * rad)) ** 0.5)
    elif op == 'mod':
        x, y, lo, hi, hf = [Q(float.fromhex(t)) for t in v[1:]]
        m2 = x * x + y * y
        if not (lo >= 0 and lo * lo <= m2 <= hi * hi): bad('mod', line)
        if abs(x) <= 2 ** 500 and abs(y) <= 2 ** 500 and not (m2 <= hf * hf): bad('mod_fast', line)
    elif op == 'pow':
        z = Q(float.fromhex(v[1])); e = int(v[2]); p = Q(float.fromhex(v[3]))
        if not z ** e <= p: bad('pow', line)
    elif op == 'log':
        x = float.fromhex(v[1]); r = float.fromhex(v[2])
        ex = mp.log(mp.mpf(x))
        cr = float(ex)   # mpmath -> float rounds to nearest (200-bit value; ties impossible for log of a double != 1)
        if r != cr: bad('log', line)
print('checked', cnt, 'violations', viol or 0, 'max ball err/rad %.3g' % worst)
