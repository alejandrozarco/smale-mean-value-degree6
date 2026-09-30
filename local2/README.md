# Independent local certificate near the regular configuration ($`d = 6`$)

Written from `SPEC.md` and the mathematics only. It uses a different method from a box subdivision: a
multiplier identity turns the constrained problem on $`E`$ into an unconstrained one, and an explicit Taylor model
with a rigorous majorant remainder handles that problem.

## Statement proved

Notation as in SPEC §1: $`d = 6`$, $`c = 5/6`$, $`\omega = e^{2\pi i/5}`$, $`b_1 = 1`$, $`b_j = \omega^{j-1} e^{\varepsilon_j}`$,
$`x = (\mathrm{Re}\,\varepsilon_2, \mathrm{Im}\,\varepsilon_2, \ldots, \mathrm{Re}\,\varepsilon_5, \mathrm{Im}\,\varepsilon_5) \in \mathbb{R}^8`$, $`L(x) = (1/5) \sum_i \log\lvert S_i(x)\rvert`$, and $`E`$ is the equal-modulus set.

**Theorem (certified).** For every $`x \in E`$ with $`\lVert x\rVert_2 \le T`$:

| $`T`$ | $`\kappa`$ | statement |
|---|---|---|
| **$`1/19`$** | **$`151/10000 = 0.0151`$** | $`L(x) \le \log(5/6) - 0.0151\lVert x\rVert_2^2`$, so $`\min_i \lvert S_i\rvert \le (5/6) e^{-0.0151\lVert x\rVert^2} \le 5/6`$ |
| $`1/8`$ | $`13/10000 = 0.0013`$ | $`L(x) \le \log(5/6) - 0.0013\lVert x\rVert_2^2`$ (extra run, not needed for the hand-over) |

In fact the certificate proves the stronger form $`\Phi(x) \le -2\kappa\lVert x\rVert^2`$, where $`\Phi = (1/5)\sum_i \lvert S_i/c\rvert^2 - 1 = \lvert S_1/c\rvert^2 - 1`$
on $`E`$. (If some $`\lvert S_i\rvert = 0`$ at a point of $`E`$ then all vanish there: $`\Phi = -1`$, $`L = -\infty`$, and the claim holds trivially. In
general $`\Phi \ge -1`$, and $`\log(1 + \Phi) \le \Phi`$ gives $`L - \log c \le \Phi/2`$ on $`E`$.) The best possible constant as $`x \to 0`$ is 0.01959, which is half the smallest curvature $`-0.039177`$ of $`\Phi`$ on
the tangent space of $`E`$. The certified $`\kappa = 0.0151`$ at $`T = 1/19`$ is 77 % of that value.

**Hand-over (SPEC §3).** Every set-aside box lies in one of the 24 balls $`\sum_j \lvert u_j - p_j\rvert^2 \le (1/20)^2`$. A relabelling
of $`b_2..b_5`$ maps that ball onto the ball around $`p_0 = (\bar\omega, \bar\omega^2, \bar\omega^3, \bar\omega^4)`$. In the $`\varepsilon`$-chart, that ball lies in
$`\lVert x\rVert_2 \le \log(20/19) = 0.051293\ldots \lt  1/19`$. The certificate therefore covers every point of every set-aside box, after
relabelling (details below).

## Method

**0. No logarithms.** On $`E`$ all $`\lvert S_i\rvert`$ are equal, so $`L - \log c = \frac{1}{2}\log(1 + \Phi) \le \Phi/2`$ with
  $`\Phi(x) = (1/5) \sum_i A_i - 1`$,  $`A_i = \lvert S_i/c\rvert^2`$.
$`E`$ is the zero set of $`D_a = A_{a+1} - A_1`$ ($`a = 1..4`$). It therefore suffices to prove $`\Phi \le -2\kappa\lVert x\rVert^2`$ on $`E`$.

**1. Exact facts at 0** (`exact_qomega.py`, exact arithmetic in $`\mathbb{Q}(\omega)`$). Expanding the product and integrating
over $`t`$ gives an exponential sum,

  $`S_i = \sum_{J\subseteq\{1..5\}} (-1)^{\lvert J\rvert}/(\lvert J\rvert+1) \cdot \omega^{e_J} \cdot \exp(\lambda_J\cdot\varepsilon)`$,  with $`\lambda_J \in \mathbb{Z}^4`$.

So every Taylor coefficient of $`S_i`$ in $`\varepsilon`$ is an explicit element of $`\mathbb{Q}(\omega)`$, and each is computed exactly.
Checked exactly:
- (F1) $`S_i(0) = 5/6`$ for all $`i`$.
- (F2) $`\sum_i \partial S_i/\partial \varepsilon_k(0) = 0`$ for $`k = 2..5`$. The derivative matrix is circulant with zero row sums.

From these, $`A_i(0) = 1`$, $`D_a(0) = 0`$, $`\Phi(0) = 0`$ and $`\nabla\Phi(0) = (2/5) \mathrm{Re}\sum_i \partial S_i/c = 0`$. These coefficients are then
set to exact zeros. The code first asserts that their arb enclosures contain 0.

**2. Multiplier identity.** For any functions $`\mu_a`$, the function $`F = \Phi + \sum_a \mu_a D_a`$ equals $`\Phi`$ on $`E`$. We take
$`\mu_a(x) = \Gamma_a\cdot x + \Psi_a(x)`$, with $`\Gamma \in \mathbb{R}^{4\times 8}`$ and $`\Psi_a`$ quadratic forms. If $`F(x) \le -2\kappa\lVert x\rVert^2`$ holds on the whole
ball $`\lVert x\rVert \le T`$, the constraint is no longer needed. Floating point is used only to choose $`\Gamma`$ and $`\Psi`$; they are then
used as exact binary rationals, and nothing about the choice has to be verified.
- $`\Gamma`$ makes the quadratic part $`G_2`$ of $`F`$ negative definite on $`\mathbb{R}^8`$. In an orthonormal basis $`(z, w)`$ of
  $`\ker \Lambda \oplus (\ker \Lambda)^\perp`$, where $`\Lambda`$ is the linear part of $`D`$, $`\Gamma`$ removes the $`z`$–$`w`$ cross block of $`\Phi_2`$ and sets the $`w`$-block
  to $`-\gamma I`$ with $`\gamma = 0.0431`$. What remains is $`\Phi_2`$ restricted to $`\ker \Lambda`$, whose maximum eigenvalue is $`-0.039177`$.
- $`\Psi`$ comes from a weighted least-squares fit. It removes the part of the cubic that is divisible by $`\Lambda x`$, which
  lowers $`\lVert G_3\rVert`$ from 1.64 to 0.10. The part that remains is essentially the intrinsic cubic of $`\Phi`$ along $`E`$.

**3. Explicit Taylor model** (`certify.py`, python-flint arb, 200 bits).
- $`\sigma_i`$ = Taylor polynomial of $`S_i/c`$ of degree $`N = 6`$ in $`\varepsilon`$. Its coefficients are the exact $`\mathbb{Q}(\omega)`$ values of step 1,
  enclosed in acb.
- $`\sigma_i`$ is split as $`U_i + iV_i`$, real polynomials in $`x`$.
- $`\alpha_i = [U_i^2 + V_i^2]_{\le 6}`$.
- The explicit part of $`\Phi`$ is $`\varphi = (1/5)\sum \alpha_i - 1`$, and the explicit part of $`D_a`$ is $`d_a = \alpha_{a+1} - \alpha_1`$.
- $`G = \varphi + \sum_a [\mu_a d_a]_{\le 6}`$. $`G`$ has no terms of degree 0 or 1.

**4. Rigorous remainders.** These are majorant series in $`s = \lVert x\rVert_2`$ with nonnegative coefficients and order $`\ge 7`$.
- Taylor tail of $`S_i`$. By Cauchy–Schwarz, $`\lvert \lambda_J\cdot\varepsilon\rvert \le \lvert \lambda_J\rvert_2 s`$. Hence
  $`\lvert S_i/c - \sigma_i\rvert \le \tau_i(s) = (6/5) \sum_J (1/(\lvert J\rvert+1)) \sum_{m\gt 6} (\lvert \lambda_J\rvert s)^m/m!`$.
  Each exponential tail is bounded by $`r^{7}/7! /(1 - r/8)`$.
- Homogeneous parts on the sphere. For a real form $`p_m`$, $`\sup_{\lvert u\rvert=1} \lvert p_m(u)\rvert \le \lVert p_m\rVert_F := (\sum_\alpha p_\alpha^2 \alpha!/m!)^{1/2}`$.
  This follows from Cauchy–Schwarz and the multinomial identity $`\sum_\alpha (m!/\alpha!) u^{2\alpha} = \lvert u\rvert^{2m}`$. The same bound
  holds for complex forms in $`\varepsilon`$, with $`\lvert \varepsilon\rvert = \lVert x\rVert_2`$.
- Error in $`A_i`$:
  $`\lvert A_i - \alpha_i\rvert \le E_i(s) = \sum_{m+m'\gt 6} n_{i,m} n_{i,m'} s^{m+m'} + 2(\sum_m n_{i,m}s^m) \tau_i + \tau_i^2`$,
  where $`n_{i,m}`$ is the Frobenius norm of the degree-$`m`$ part of $`\sigma_i`$.
- Error in $`F`$:
  $`\lvert F - G\rvert \le \mathrm{err}F(s) = (1/5)\sum E_i + \sum_a (\lvert \Gamma_a\rvert s + \lVert \Psi_a\rVert_F s^2)(E_{a+1} + E_1)`$ + (bound on the truncated part of $`\mu_a d_a`$).

  Every term of $`\mathrm{err}F(s)`$ is a power series in $`s`$ with nonnegative coefficients and order $`\ge 7`$ (the truncation
  degree is 6), so $`\mathrm{err}F(s)/s^2`$ is nondecreasing, and it is evaluated once, at $`s = T`$. Check:
  $`\mathrm{err}F(T)/\mathrm{err}F(T/2) \approx 133 \approx 2^{7}`$.

**5. Final inequality.** For $`\lvert u\rvert = 1`$ and $`0 \lt  s \le T`$:

  $`F(su)/s^2 \le G_2(u) + \sum_{m=3}^{6} \lVert G_m\rVert_F s^{m-2} + \mathrm{err}F(s)/s^2 \le -\lambda + \sum_m \lVert G_m\rVert_F T^{m-2} + \mathrm{err}F(T)/T^2 = -2\kappa`$.

Here $`G_2 \le -\lambda\lvert x\rvert^2`$ is certified by an arb Cholesky factorisation of $`-H_2 - \lambda I`$, with every pivot certainly positive.
Every comparison is a certain arb comparison.

## Certified constants ($`T = 1/19`$, $`N = 6`$; from `run_T1_19.log`)

| quantity | value |
|---|---|
| $`\lambda`$ ($`G_2 \le -\lambda\lVert x\rVert^2`$, Cholesky-certified) | $`42023795/2^{30} \approx 0.0391377`$ (exact value $`0.039137708954\ldots`$) |
| $`\lVert G_3\rVert_F \cdot T`$ | 0.00534424 |
| $`\lVert G_4\rVert_F \cdot T^2`$ | 0.00329886 |
| $`\lVert G_5\rVert_F \cdot T^3`$ | 0.000233792 |
| $`\lVert G_6\rVert_F \cdot T^4`$ | 1.55693e-5 |
| $`\tau_i(T)`$ (Taylor tail of $`S_i/c`$) | $`\le`$ 5.7540e-9 |
| $`E_i(T)`$ (error in $`\lvert S_i/c\rvert^2`$) | $`\le`$ 1.9643e-8 |
| $`\mathrm{err}F(T)/T^2`$ | 2.01052e-5 |
| $`2\kappa = \lambda - \sum - \mathrm{err}F/T^2`$ | $`[0.030225144 \pm 3.4\text{e-}11]`$ |
| $`\kappa`$ (certified lower bound) | $`\gt  151/10000`$ |

At $`T = 1/8`$: $`2\kappa = [0.0026125367 \pm 3.3\text{e-}11]`$ and $`\kappa \gt  13/10000`$. The limiting term is $`\lVert G_4\rVert_F T^2`$.

## Hand-over (`handover.py`)

1. **Relabelling.** Take a permutation $`(a_2..a_5)`$ of $`(1..4)`$. Let $`\tau(1) = 1`$ and $`\tau(k) =`$ the $`j`$ with $`a_j = k-1`$, and
   put $`b'_k = b_{\tau(k)}`$. Then $`u'_k = u_{\tau(k)}`$, and $`\sum_k \lvert u'_k - \bar\omega^{k-1}\rvert^2 = \sum_j \lvert u_j - \bar\omega^{a_j}\rvert^2`$. The ball around $`p`$ is
   therefore mapped onto the ball around $`p_0`$, and the inverse relabelling maps it back.
   Also $`S_k(b') = \int \prod_j(1 - t\, b_{\tau(k)}/b_{\tau(j)})\, dt = S_{\tau(k)}(b)`$, because $`\tau`$ is a bijection of $`\{1..5\}`$. So the $`S_i`$
   are permuted, and $`E`$ and $`\min_i\lvert S_i\rvert`$ are invariant. The index identities are checked for all 24 permutations,
   and a float spot check of the $`S`$-identity is in `crosscheck.py`.
2. **Chart.** Put $`w_j = u_j/p0_j - 1`$, so $`\lvert w_j\rvert = \lvert u_j - p0_j\rvert \le 1/20`$. The choice $`\varepsilon_j = -\mathrm{Log}(1 + w_j)`$ gives
   $`b_j = 1/u_j = \omega^{j-1}e^{\varepsilon_j}`$. Also $`\lvert \varepsilon_j\rvert \le -\log(1 - \lvert w_j\rvert) \le C\lvert w_j\rvert`$, with $`C = -20 \log(19/20) = 1.02587`$, because
   $`-\log(1-q)/q`$ is increasing. Hence $`\lVert x\rVert_2^2 = \sum \lvert \varepsilon_j\rvert^2 \le C^2/400`$, and $`\lVert x\rVert_2 \le \log(20/19) = 0.0512933 \lt  1/19`$ (arb-certified).
3. **Covered set.** The certificate covers $`X = \{b = (1, \omega e^{\varepsilon_2}, \omega^2 e^{\varepsilon_3}, \omega^3 e^{\varepsilon_4}, \omega^4 e^{\varepsilon_5}) : \lVert x\rVert_2 \le 1/19\}`$.
   On $`X \cap E`$, $`\min_i\lvert S_i\rvert \le (5/6)e^{-0.0151\lVert x\rVert^2}`$. Take a set-aside box inside the ball around some $`p`$. Its relabelled
   image lies in the ball around $`p_0`$, which lies in $`X`$ by step 2. $`E`$ and $`\min\lvert S_i\rvert`$ are invariant under relabelling, so
   the bound holds on every set-aside box.

## How to run

`certify.py` reads the floating-point design ($`\Gamma`$, $`\Psi`$) from `design_T*.json` when that file exists, so the certified
constants do not depend on the numpy/BLAS build. Delete the file to regenerate the design.


```
bash local2/run_all.sh        # writes local2/run_T1_19.log
```
The script runs, at nice 15 on a single thread, in this order: the main certificate ($`T = 1/19`$), the hand-over,
the extra certificate ($`T = 1/8`$), and a non-rigorous float cross-check. **Run time: about 12 s wall in total; each
certificate takes under 1 s** (peak memory about 55 MB). Single runs:
`python3 certify.py 1/19 6`, `python3 handover.py`, `python3 crosscheck.py 1/19 6`.

Files:
- `exact_qomega.py`: $`\mathbb{Q}(\omega)`$ arithmetic, the exact Taylor coefficients, and facts F1–F2.
- `certify.py`: the certificate.
- `handover.py`: SPEC §3.
- `crosscheck.py`: an independent float evaluator of $`S_i`$, used to compare the explicit polynomial $`G`$ with
  $`F`$ computed directly. The differences scale like $`s^7`$ and stay 2–3 orders below $`\mathrm{err}F`$.
- `design_T*.json`: the multipliers $`\Gamma`$ and $`\Psi`$ that were used, as exact float hex values.
- `run_T1_19.log`: the log of a full run.

## What is and is not certified

- Certified: facts F1 and F2 (exact rational arithmetic in $`\mathbb{Q}(\omega)`$), every Taylor coefficient (exact, then
  enclosed in arb), every remainder bound, the Cholesky step, the final inequalities, and the hand-over constant.
- Floating point only chooses $`\gamma`$, $`\Gamma`$, $`\Psi`$ and $`\lambda`$. Every one of these is used as an exact number afterwards.
- `crosscheck.py` is a sanity check only and is not part of the proof.
- Trusted: python-flint 0.6 / FLINT arb ball arithmetic being correct, and the short analytic arguments written
  above: the exponential-sum form, Cauchy–Schwarz, the majorant monotonicity, and the relabelling identity.
