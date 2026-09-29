"""exact_facts.py — machine check, in EXACT arithmetic, of the two low-order facts used by the local certificate
(local/LOCAL_CERT.md step 2), which were previously proved by hand only:

  (a) S_i(0) = (d-1)/d for every i          (so h_i(0) = log S_i(0) - log c = 0 and L(0) = log c exactly);
  (b) sum_i dS_i/d eps_j (0) = 0 for every j = 2..n   (complex derivative; hence grad L(0) = 0 exactly, because
      L = (1/n) sum_i Re log S_i, d(Re h)/dx = Re h', d(Re h)/dy = -Im h', and h_i' = S_i'/S_i(0) = S_i'/c).

Chart: b_1 = 1, b_j = omega^{j-1} e^{eps_j} (j = 2..n), omega = e^{2 pi i/n}, n = d-1;
S_i = int_0^1 (1-t) prod_{j != i} (1 - t z_j) dt = sum_k (-1)^k e_k(z) / ((k+1)(k+2)),  z_j = omega^{i-j} e^{eps_i - eps_j}.
dz_m/d eps_j = z_m (delta_{ij} - delta_{mj}) (m != i),  d e_k / d z_m = e_{k-1}(z without z_m).
All values at eps = 0 lie in Q(omega) = Q[x]/(Phi_n(x)); we compute with Fraction coefficients modulo the
cyclotomic polynomial Phi_n, so equality to 0 or to (d-1)/d is decided exactly.

usage: python3 local/exact_facts.py [d ...]      (default: 5 6)
"""
import sys
from fractions import Fraction as Q

PHI = {3: [1, 1, 1], 4: [1, 0, 1], 5: [1, 1, 1, 1, 1], 6: [1, -1, 1]}   # coefficients, low degree first (monic)


class Cyc:
    """element of Q(omega_n) as a coefficient list of length deg Phi_n (low degree first)."""
    def __init__(self, n, c):
        self.n, self.c = n, c

    @staticmethod
    def reduce(n, p):
        phi = PHI[n]; m = len(phi) - 1
        p = list(p)
        for k in range(len(p) - 1, m - 1, -1):          # x^m = -(phi_0 + ... + phi_{m-1} x^{m-1})
            a = p[k]
            if a:
                p[k] = 0
                for l in range(m):
                    p[k - m + l] -= a * phi[l]
        return Cyc(n, (p + [Q(0)] * m)[:m])

    @staticmethod
    def const(n, a):
        return Cyc.reduce(n, [Q(a)])

    @staticmethod
    def root(n, k):                                      # omega^k
        k %= n
        return Cyc.reduce(n, [Q(0)] * k + [Q(1)])

    def __add__(self, o):
        return Cyc(self.n, [a + b for a, b in zip(self.c, o.c)])

    def __sub__(self, o):
        return Cyc(self.n, [a - b for a, b in zip(self.c, o.c)])

    def __mul__(self, o):
        if not isinstance(o, Cyc):
            return Cyc(self.n, [a * Q(o) for a in self.c])
        p = [Q(0)] * (len(self.c) + len(o.c))
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    p[i + j] += a * b
        return Cyc.reduce(self.n, p)

    def is_zero(self):
        return all(a == 0 for a in self.c)

    def __eq__(self, o):
        return (self - o).is_zero()


def elem_sym(zs, n):
    """e_0..e_len(zs) of the list zs (Cyc)."""
    E = [Cyc.const(n, 1)] + [Cyc.const(n, 0) for _ in zs]
    for z in zs:
        for k in range(len(zs), 0, -1):
            E[k] = E[k] + E[k - 1] * z
    return E


def check(d):
    n = d - 1
    ck = [Q((-1) ** k, (k + 1) * (k + 2)) for k in range(n)]
    out = {'d': d}
    S0 = []
    for i in range(n):
        zs = [Cyc.root(n, i - j) for j in range(n) if j != i]
        E = elem_sym(zs, n)
        S = Cyc.const(n, 0)
        for k in range(n):
            S = S + E[k] * ck[k]
        S0.append(S == Cyc.const(n, Q(d - 1, d)))
    out['S_i(0) == (d-1)/d for all i'] = all(S0)
    grad_zero = []
    for jj in range(1, n):                     # eps_j, j = 2..n  (index jj = j-1)
        tot = Cyc.const(n, 0)
        for i in range(n):
            idx = [j for j in range(n) if j != i]
            zs = [Cyc.root(n, i - j) for j in idx]
            for mpos, m in enumerate(idx):
                dz = (1 if i == jj else 0) - (1 if m == jj else 0)     # dz_m/d eps_j = z_m (delta_ij - delta_mj)
                if dz == 0:
                    continue
                rest = zs[:mpos] + zs[mpos + 1:]
                Er = elem_sym(rest, n)
                # dS_i/dz_m = sum_k c_k e_{k-1}(rest)
                dS = Cyc.const(n, 0)
                for k in range(1, n):
                    dS = dS + Er[k - 1] * ck[k]
                tot = tot + dS * zs[mpos] * dz
        grad_zero.append(tot.is_zero())
    out['sum_i dS_i/deps_j(0) == 0 for j=2..n'] = all(grad_zero)
    out['ok'] = all(S0) and all(grad_zero)
    return out


if __name__ == '__main__':
    ds = [int(a) for a in sys.argv[1:]] or [5, 6]
    allok = True
    for d in ds:
        r = check(d)
        print(r)
        allok &= r['ok']
    print('EXACT FACTS VERIFIED' if allok else 'EXACT FACTS FAILED')
    sys.exit(0 if allok else 1)
