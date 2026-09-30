# Local certificate at the equality configuration ($`d = 6`$; also $`d = 5`$)

The code is `local/local_cert.py`, which uses python-flint arb/acb balls at 128 bits, together with
`local/exact_facts.py`, which does exact arithmetic in $`\mathbb{Q}(\omega)`$. The run logs are
`runs/local_cert_v3_d6_T0.0527.log` and `runs/local_cert_v3_d5_T0.0527.log`; each takes about 1 minute. The logs
`runs/local_cert_v3_d6_T{0.01,0.02,0.03,0.04}.log` are the same computation at smaller radii (used in the figure).

## Statement

Let $`d = 6`$, $`n = 5`$, $`c = 5/6`$. Chart: $`b_1 = 1`$ and $`b_j = \omega^{j-1} e^{\varepsilon_j}`$ ($`j = 2, \dots, 5`$,
$`\omega = e^{2\pi i/5}`$), with $`\varepsilon \in \mathbb{C}^4`$,
$`x = (\mathrm{Re}\varepsilon_2, \mathrm{Im}\varepsilon_2, \dots) \in \mathbb{R}^8`$ and $`t = \lVert x\rVert_2`$.
Put $`\ell_i = \log\lvert S_i\rvert`$, $`L = \frac1n \sum_i \ell_i`$, and
$`E = \{\lvert S_1\rvert = \dots = \lvert S_5\rvert\}`$.

**Claim.** For every $`x \in E`$ with $`0 \lt  \lVert x\rVert_2 \le T := 0.0527`$,

```math
\log F(x) = L(x) \le \log\tfrac56 - 0.00757\,\lVert x\rVert_2^2 \lt  \log\tfrac56 .
```

At $`x = 0`$, $`F = 5/6`$. $`T`$ is the exact decimal, and the code uses an arb ball that contains it.

**Hand-over to the branch and bound.** `c/smale_bb_v2` (excl_mode 2) sets a box aside only when it has verified that
the box lies in a closed Euclidean $`u`$-ball $`\sum_j \lvert u_j - p_j\rvert^2 \le r^2`$ ($`r = 0.05`$) around one of the 24
equality points $`p`$. Relabelling maps each such ball to the ball around
$`p = (\bar\omega, \bar\omega^2, \bar\omega^3, \bar\omega^4)`$, where $`u_j = \omega^{-(j-1)} e^{-\varepsilon_j}`$. With
$`w_j = u_j/p_j - 1`$ we have $`\lvert\varepsilon_j\rvert = \lvert\log(1 + w_j)\rvert \le \lvert w_j\rvert/(1 - \lvert w_j\rvert)`$,
so $`\lVert\varepsilon\rVert_2 \le r/(1 - r) = 1/19 \lt  T`$. `c/check_done.py` checks $`r/(1-r) \le 0.0527`$ exactly for
the IEEE value of $`r`$.

## Proof structure (all constants computed in arb)

1. **Taylor models.** $`h_i = \log S_i - \log c`$ is holomorphic on the polydisk $`\lVert\varepsilon\rVert_\infty \le T`$.
   Taylor-model arithmetic (degree $`N = 6`$ in 4 complex variables) gives $`h_i = P_i + R_i`$. Here $`P_i`$ is the
   degree-$`\le 6`$ Taylor polynomial with acb ball coefficients, and $`\lvert R_i\rvert \le r_i`$ on the polydisk
   ($`r_i \le 2.9\cdot10^{-6}`$). The code uses $`\lvert R_i\rvert \le r_i (t/T)^2`$ for $`t \le T`$.
2. **Exact low-order facts** (`local/exact_facts.py`, asserted by `local_cert.py`). In $`\mathbb{Q}(\omega)`$:
   $`S_i(0) = (d-1)/d`$ for all $`i`$, and $`\sum_i \partial S_i/\partial\varepsilon_j(0) = 0`$ for $`j = 2, \dots, n`$. Hence
   $`h_i(0) = 0`$, $`L(0) = \log c`$ and $`\nabla L(0) = 0`$ exactly (checked for $`d = 4, \dots, 7`$: `python3 local/exact_facts.py 4 5 6 7`). With $`J = D\ell(0)`$
   (a $`5\times 8`$ matrix whose rows sum to $`0`$), $`G_i := \ell_i - L = (Jx)_i + G_i^{\ge 2}(x)`$.
3. **Tube around $`K = \ker J`$.** On $`E`$ all $`G_i = 0`$, so $`y := \lVert Jx\rVert = \lVert G^{\ge 2}(x)\rVert`$. Write
   $`x = \Pi x + w`$ with $`\Pi = VV^{\mathsf T}`$ (formed in arb) and $`\lVert w\rVert \le \beta y + \gamma t`$, where
   $`\beta = \lVert J^+\rVert \le 1.7929`$ ($`J^+`$ a floating-point pseudo-inverse; its defect is accounted for by $`\gamma`$) and $`\gamma \le 2.04\cdot10^{-15}`$. With:
   - $`q = \sup_{\lVert v\rVert \le 1} \lVert G^2(\Pi v)\rVert \le 0.23673`$ (a rigorous branch-and-bound over the unit
     sphere in $`K`$);
   - $`b_2 = \bigl(\sum_i \lVert A_i\rVert^2\bigr)^{1/2} \le 1.2650`$;
   - the degree $`3..6`$ parts and remainders of $`G`$,

   the fixed-point inequality
   $`y \le q t^2 + 2 b_2 \lVert\Pi\rVert\, t(\beta y + \gamma t) + b_2 (\beta y + \gamma t)^2 + g_3(t)`$,
   started from $`y \le \lVert J\rVert t`$, gives $`\lVert Jx\rVert \le \kappa t^2`$ with $`\kappa = 0.4787`$ for all $`x \in E`$
   with $`\lVert x\rVert \le T`$.
4. **Quadratic part.** For $`\sigma = 5`$, $`H_\sigma := \mathrm{Hess} L(0) - \sigma J^{\mathsf T} J \preceq -\mu I`$ with
   $`\mu = 0.03192459`$ (arb $`LDL^{\mathsf T}`$ of $`-H_\sigma - \mu I`$). Hence
   $`\tfrac12 x^{\mathsf T} \mathrm{Hess} L(0)\, x \le -\tfrac{\mu}{2} t^2 + \tfrac{\sigma}{2}\kappa^2 t^4`$.
5. **Cubic part.**
   - $`\sup_{\lVert z\rVert = 1} \lvert L^3(Vz)\rvert \le 0.0325473`$, by a rigorous sphere branch-and-bound with exact
     rational sphere-intersection tests.
   - $`\lvert L^3(a + w) - L^3(a)\rvert \le C_3\bigl((\lVert a\rVert + \lVert w\rVert)^3 - \lVert a\rVert^3\bigr)`$ with
     $`C_3 \le 0.5192`$.
   - Degrees $`4..6`$ are bounded by $`\sum_\alpha \lvert C_\alpha\rvert\, m_\alpha t^k`$, with
     $`m_\alpha = \prod_j (\alpha_j/k)^{\alpha_j/2}`$ computed in arb.
   - The remainder contributes $`\frac1n \sum_i r_i (t/T)^2`$.
6. **Conclusion.** $`L(x) - \log c \le t^2 g(t)`$, where every term of $`g`$ is increasing in $`t`$:
   - quadratic: $`-0.015962`$;
   - penalty: $`0.001591`$;
   - cubic on $`K`$: $`0.001715`$;
   - cubic off $`K`$: $`0.003883`$;
   - degrees 4–6: $`0.000782`$;
   - remainder: $`0.000417`$.

   So $`g(T) \le -0.0075740986`$ (an arb upper bound; the displayed terms are rounded).

For $`d = 5`$ the same code gives $`g(T) \le -0.02284`$ ($`\mu = 0.050387`$).

## Numerical comparison (not rigorous)

`local/e_profile.py` projects points onto $`E`$ by Newton's method, in float64, from 400 random $`K`$-directions per radius.
The sampled minimum of $`(\log c - \log F)/\lVert x\rVert^2`$ on $`E`$ is $`0.0198`$–$`0.0215`$ for $`t = 0.005, \dots, 0.2`$
(`runs/e_profile_d6.tsv`). `local/check_tm.py` samples $`\lvert h_i - P_i\rvert`$ at random points of the polydisk and
compares it with $`r_i`$.
