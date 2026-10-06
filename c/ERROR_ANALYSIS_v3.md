# Floating-point error analysis of `c/smale_bb_v3.c`: addendum to `c/ERROR_ANALYSIS_v2.md`

`c/smale_bb_v3.c` is `c/smale_bb_v2.c` plus change A9: a second model of $`S_i`$ ($`i \ge 2`$) in the variable $`w = 1/u_i`$.
Apart from A9, the two sources differ in version strings and in A10, a run-mode option that restricts a run to a
subset of the task ids (no effect on any computed bound; see W.9). Every statement of `c/ERROR_ANALYSIS_v2.md` (A.1–A.8)
therefore applies to v3 unchanged. That covers the trust base, the domain, ball arithmetic, the T-form enclosures and
models, the tests, geometry and exclusion. This addendum covers A9 only: the functions `taylor_S` and `wmodel`, and
the eleven lines at the end of `build_models` that combine the two models. It uses the notation of v2's analysis:
$`U = 2^{-53}`$, $`\eta_0 = 2^{-1075}`$, $`\mathrm{up}(x) = \mathrm{fl}(\mathrm{fl}(x\cdot\mathsf{INFL}) + \mathsf{TINY})`$, facts (F1)–(F6), and lemmas L1–L5, R, S, LB and A.

## W.0 The identity
For $`i \ge 2`$, $`b_1 = 1`$, $`b_j = 1/u_j`$ and $`w := 1/u_i`$:

```math
S_i = \int_0^1 \prod_j \Bigl(1 - t\,\frac{b_i}{b_j}\Bigr)dt = \int_0^1 (1-t)\,(1 - t w) \prod_{j \ne 1,i} (1 - t\,u_j\,w)\,dt ,
```

because $`b_i/b_1 = w`$, $`b_i/b_j = u_j w`$, and the factor $`j = i`$ is $`(1-t)`$. Multiplying by $`u_i^{\,n-1} = w^{-(n-1)}`$ gives v2's
$`T_i = \int_0^1 (1-t)(u_i - t)\prod_{j\ne 1,i}(u_i - t u_j)\,dt`$. As a polynomial, $`S_i`$ has degree $`\le n-1`$ in $`w`$ and degree $`\le 1`$
in each $`u_j`$, $`j \ne 1, i`$. These are the degree bounds v2 uses for $`T_i`$ in $`u_i`$ and the $`u_j`$. The identity was
checked exactly with sympy for $`d = 4..7`$ (`firstprinciples/derive.py`, item 2, with $`V_i = P(b_i)/b_i`$).

## W.1 Domain of the w-model and the radius $`\rho_w`$
Write $`v = i - 1`$, $`m = m_v`$ (exact doubles, the box centre), $`\delta = u_i - m`$, and use the computed values
$`r = \rho_v \ge \lvert\delta\rvert`$ on $`B`$ (L3) and $`\mathrm{am} = \mathrm{absm}_{\mathrm{lo}}[v] \le \lvert m\rvert`$ (L3).

**Lemma W1.** The w-model is computed only if $`\mathrm{am} \gt  2r`$ (an exact comparison of doubles). Then on $`B`$:
- $`\lvert u_i\rvert \ge \lvert m\rvert - r \gt  \lvert m\rvert/2 \gt  0`$, so $`w = 1/u_i`$ is defined;
- $`\lvert w - 1/m\rvert = \lvert\delta\rvert/(\lvert u_i\rvert\lvert m\rvert) \le r/((\mathrm{am} - r)\,\mathrm{am}) =: \rho_w^{\ast}`$;
- the computed $`\mathrm{rw} = \mathrm{up}(\mathrm{fl}(r/\mathrm{fl}(\mathrm{am}\cdot\mathrm{fl}(\mathrm{am}-r))))`$ satisfies $`\mathrm{rw} \ge \rho_w^{\ast}`$.

*Proof.* The first two items are exact real inequalities. The bound $`\rho/((\lvert m\rvert - \rho)\lvert m\rvert)`$ increases in $`\rho`$ and
decreases in $`\lvert m\rvert`$, so it may be evaluated at $`r \ge \rho`$ and $`\mathrm{am} \le \lvert m\rvert`$.

Magnitudes (A.2, A5): $`H_k \ge 2^{-40}`$ gives $`r \ge 2^{-40}`$. With $`\mathrm{am} \gt  2r`$, also $`\mathrm{am} - r \gt  r \ge 2^{-40}`$. Further,
$`\mathrm{am} \le \lvert m\rvert \le \sqrt2`$, so $`\mathrm{am}(\mathrm{am} - r) \in [2^{-80}, 2]`$ and $`\rho_w^{\ast} \in [2^{-41}, 2^{39}]`$. All intermediate values are
normal. So each of the three roundings has relative error $`\le U`$, and the computed quotient is
$`\ge \rho_w^{\ast}(1-U)/(1+U)^2 \ge \rho_w^{\ast}(1 - 3U)`$. By (F6), $`\mathrm{up}(\cdot)`$ restores $`\ge \rho_w^{\ast}`$, since $`\mathsf{INFL} - 1 = 2^{-40} \gg 3U`$. ∎

## W.2 The centre ball $`w_0 \ni 1/m`$
**Lemma W2.** Let $`\mathrm{m2} = \mathrm{fl}(\mathrm{fl}(m_{re}^2) + \mathrm{fl}(m_{im}^2))`$, $`w_{0,re} = \mathrm{fl}(m_{re}/\mathrm{m2})`$, $`w_{0,im} = \mathrm{fl}(-m_{im}/\mathrm{m2})`$ and
$`w_{0,\mathrm{rad}} = \mathrm{fl}(\mathrm{fl}(8U\cdot\mathrm{fl}(\lvert w_{0,re}\rvert + \lvert w_{0,im}\rvert)) + \mathsf{TINY})`$. Then $`\lvert (w_{0,re} + i\,w_{0,im}) - 1/m\rvert \le w_{0,\mathrm{rad}}`$.

*Proof.* Each nonzero component of $`m`$ is $`\ge 2^{-40}`$ in magnitude (A5: an integer multiple of $`H \ge 2^{-40}`$), and $`\lvert m\rvert \le \sqrt2`$. So
both squares are $`0`$ or normal, and $`\mathrm{m2} = \lvert m\rvert^2(1 + \theta)`$ with $`\lvert\theta\rvert \le 2U + U^2`$. The exact value is
$`1/m = (m_{re} - i\,m_{im})/\lvert m\rvert^2`$. Each computed component equals the exact one times $`(1+\delta')/(1+\theta)`$ with
$`\lvert\delta'\rvert \le U`$, so its relative error is $`\le 3.1U`$. The error vector therefore has modulus
$`\le 3.1U(\lvert m_{re}\rvert + \lvert m_{im}\rvert)/\lvert m\rvert^2 \le 3.2U(\lvert w_{0,re}\rvert + \lvert w_{0,im}\rvert)`$.

The radius computation has two roundings and one addition of $`\mathsf{TINY}`$. It returns at least
$`8U(\lvert w_{0,re}\rvert + \lvert w_{0,im}\rvert)(1-U)^2 \gt  3.2U(\ldots)`$. ∎

## W.3 Taylor coefficients (`taylor_S`)
**Lemma TS.** The balls $`\mathrm{coef}[a][mk]`$ contain the exact Taylor coefficients of $`S_i`$ at the point
$`(w, u_j) = (1/m, m_j)`$, in the variables $`\delta_w = w - 1/m`$ (power $`a`$) and $`\delta_j`$ (multilinear, bit mask $`mk`$).

*Proof.* `taylor_S` multiplies the factors

```math
1 - t w = 1 - t w_0 - t\,\delta_w \qquad\text{and}\qquad 1 - t u_j w = 1 - t\,m_j w_0 - t\,m_j\,\delta_w - t\,w_0\,\delta_j - t\,\delta_j\delta_w
```

as polynomials in $`(t, \delta_w, \delta_j)`$. Here $`w_0`$ is the ball of W2, and $`m_j w_0`$ is formed by `cb_mul`. It then
integrates in $`t`$ with the moment balls $`1/((k+1)(k+2))`$ (`P->mom`, A.3). Each coefficient is a polynomial in $`w_0`$,
the $`m_j`$ and the moments, evaluated in ball arithmetic. By L1 and the inclusion property, the result contains the value
of that polynomial at the exact point $`w_0 = 1/m`$. This is v2's Lemma T, with one more input ball.

Magnitudes, for L1's no-overflow hypothesis:
- $`\lvert w_0\rvert \le 1/\mathrm{am}\,(1 + 4U) \lt  2^{40}`$.
- Every factor's coefficients are bounded in modulus by $`1 + \lvert m_j w_0\rvert + \lvert m_j\rvert + \lvert w_0\rvert + 1 \lt  2^{42}`$.
- With at most $`n - 1 \le 5`$ factors, all coefficient midpoints and radii stay below $`2^{215}`$, far below the $`2^{500}`$
  limit of L2.
- Underflow in a radius is covered by $`\mathsf{TINY}`$ per operation (L1).

Array bounds: the $`w`$-degree is at most $`1 + (n-2) = n - 1 \le 5 \lt  \mathsf{NMAX}`$, and the $`t`$-degree at most $`n - 1 \le 5 \lt  \mathsf{TDEG}`$. ∎

## W.4 Enclosures (interval part of `wmodel`)
**Lemma WR.** Lemma R of v2 holds for $`S_i`$ in the w-form, with the radius vector $`\mathrm{rr}`$ (where $`\mathrm{rr}[v] = \mathrm{rw}`$ and
$`\mathrm{rr}[k] = \rho_k`$ otherwise) in place of $`\rho`$. On $`B`$, $`\lvert\delta_w\rvert \le \rho_w^{\ast} \le \mathrm{rw}`$ (W1) and $`\lvert\delta_j\rvert \le \rho_j`$.

*Proof.* The code is v2's loop with $`\mathrm{rr}`$ in place of $`\rho`$; $`\mathrm{ra} = \mathrm{powu}(\mathrm{rw}, a)`$ (L4). L4 requires
$`\mathrm{rw} \ge 2^{-60}`$, and $`\mathrm{rw} \ge 2^{-41}`$ by W1. Each product rp has at most five factors $`\mathrm{rw} \le 2^{39}`$ and at most four other radii $`\le 2`$ (since $`\lvert C_k\rvert + H_k \le 1`$).
So rp $`\le 2^{195}\cdot 2^{4} \lt  2^{200}`$, and there is no overflow. The
summand count is v2's, since the layout is identical, and v2's rounding argument applies verbatim. ∎

**Lemma WS/WLB.**
- $`\mathrm{supT} \ge \sup_B \lvert S_i\rvert`$ and $`\mathrm{infT} \le \inf_B \lvert S_i\rvert`$. This is v2's Lemma S, same code.
- $`\mathrm{hi} = \mathrm{fl}(\mathrm{ls} + \mathrm{fl}(\mathsf{REL}\cdot\mathrm{fl}(1 + \lvert\mathrm{ls}\rvert))) \ge \log\sup_B\lvert S_i\rvert`$, where $`\mathrm{ls} = \mathrm{cr\_log}(\mathrm{supT})`$: the log loses $`\le U\lvert\mathrm{ls}\rvert`$
  (L5), and the allowance $`\ge 890U(1 + \lvert\mathrm{ls}\rvert)`$ covers it together with the final rounding.
- $`\mathrm{lo}`$ is symmetric when $`\mathrm{infT} \gt  0`$; otherwise $`\mathrm{lo} = -\infty`$.
- These are bounds for $`\log\lvert S_i\rvert`$ itself. Unlike v2's LB, there is no $`(n-1)\log\lvert u_i\rvert`$ term, because the w-model
  is only formed when $`u_i`$ is bounded away from 0 on $`B`$ (W1). ∎

## W.5 The affine model (`wmodel`, model part)
Assume the model is marked valid: grad_ok ($`\mathrm{c0n2} \ge 2^{-900}`$), $`\mathrm{A0lo} \gt  0`$, $`q \lt  1/2`$, and all quantities finite.

**Exact model.** As in v2's A.6, with $`\gamma_\alpha = \hat c_\alpha/\hat c_0`$:

```math
\log\lvert S_i\rvert = \log\lvert\hat c_0\rvert + \mathrm{Re}\Bigl(\gamma_w\delta_w + \sum_{j}\gamma_j\delta_j\Bigr) \pm \Bigl(\frac{\mathrm{rT}}{\mathrm{A0lo}} + \frac{q^2}{2(1-q)}\Bigr)\quad\text{on } B .
```

The own variable is substituted exactly:

```math
\delta_w = \frac{1}{m+\delta} - \frac1m = s\,\delta + \eta,\qquad s = -\frac{1}{m^2},\qquad \eta = \frac{\delta^2}{m^2\,(m+\delta)},\qquad
\lvert\eta\rvert \le \frac{r^2}{\mathrm{am}^2(\mathrm{am} - r)} =: \eta^{\ast}.
```

So $`\mathrm{Re}(\gamma_w\delta_w) = \mathrm{Re}(\gamma_w s\,\delta) \pm \lvert\gamma_w\rvert\,\eta^{\ast}`$, and the exact affine model has gradient $`\gamma_w s`$ in the variable $`\delta`$
and $`\gamma_j`$ in $`\delta_j`$. Its error is

```math
e^{\ast} = \mathrm{rT}/\mathrm{A0lo} + q^2/(2(1-q)) + \lvert\gamma_w\rvert\,\eta^{\ast} .
```

**Lemma AW (the computed model is valid).** For all $`\delta \in B - m`$: $`\bigl\lvert\log\lvert S_i\rvert - a - G\cdot(x,y)\bigr\rvert \le e`$, where $`(a, G, e)`$ are the
computed quantities. Error sources, all charged to $`e`$:

1. **$`\gamma_j`$, $`j \ne i`$.** Exactly as in v2's A.6, item 1: $`\lvert\gamma_c - \gamma\rvert \le 11U(\lvert\mathrm{gr}_c\rvert + \lvert\mathrm{gi}_c\rvert) + 2^{-170}`$, charged
   $`(\mathsf{REL}(\lvert\mathrm{gr}_c\rvert + \lvert\mathrm{gi}_c\rvert) + 2^{-160})\,\rho_j`$.

2. **$`s = -1/m^2`$.** The code computes $`\mathrm{m4} = \mathrm{fl}(\mathrm{m2}^2)`$, $`\mathrm{sre} = \mathrm{fl}(-\mathrm{fl}(\mathrm{fl}(m_{re}^2) - \mathrm{fl}(m_{im}^2))/\mathrm{m4})`$ and
   $`\mathrm{sim} = \mathrm{fl}(\mathrm{fl}(2 m_{re} m_{im})/\mathrm{m4})`$. The exact values are $`s_{re} = -(m_{re}^2 - m_{im}^2)/\lvert m\rvert^4`$ and $`s_{im} = 2m_{re}m_{im}/\lvert m\rvert^4`$.
   - The numerator of $`s_{re}`$ may cancel. Its absolute error is $`\le 3U(m_{re}^2 + m_{im}^2) = 3U\lvert m\rvert^2`$.
   - The numerator of $`s_{im}`$ has relative error $`\le U`$; $`2\times`$ is exact.
   - Division by $`\mathrm{m4} = \lvert m\rvert^4(1 + \theta_4)`$, $`\lvert\theta_4\rvert \le 5.1U`$, adds a relative $`\le 6.2U`$.

   Hence $`\lvert s_c - s\rvert_1 \le 3U/\lvert m\rvert^2 + 7.3U(\lvert s_{re}\rvert + \lvert s_{im}\rvert)(1 + O(U)) \le 11U\,\lvert s\rvert_1`$, using $`\lvert s\rvert_1 \ge \lvert s\rvert = 1/\lvert m\rvert^2`$.
   All quantities are normal, since $`\lvert m\rvert \ge 2^{-40}`$.

3. **Own-variable gradient $`t = \gamma_w s`$.** The code computes $`\mathrm{tr} = \mathrm{fl}(\mathrm{fl}(\mathrm{gr}\cdot\mathrm{sre}) - \mathrm{fl}(\mathrm{gi}\cdot\mathrm{sim}))`$ and
   $`\mathrm{ti} = \mathrm{fl}(\mathrm{fl}(\mathrm{gr}\cdot\mathrm{sim}) + \mathrm{fl}(\mathrm{gi}\cdot\mathrm{sre}))`$. With $`\gamma_c = \mathrm{gr} + i\,\mathrm{gi}`$:

   ```math
   \lvert t_c - \gamma_w s\rvert \le \lvert t_c - \gamma_c s_c\rvert + \lvert\gamma_c\rvert\,\lvert s_c - s\rvert + \lvert\gamma_c - \gamma_w\rvert\,\lvert s\rvert .
   ```

   - The first term is $`\le 2.1U\,\lvert\gamma_c\rvert_1\lvert s_c\rvert_1`$, from the products and one sum per component.
   - The second is $`\le 11U\,\lvert\gamma_c\rvert_1\lvert s\rvert_1`$ by item 2.
   - The third is $`\le (11U\lvert\gamma_c\rvert_1 + 2^{-170})\lvert s\rvert`$ by item 1.

   For complex numbers $`\lvert z\rvert_1 \le \sqrt2\,\lvert z\rvert`$, so $`\lvert\gamma_c\rvert_1\lvert s_c\rvert_1 \le 2\lvert\gamma_c s_c\rvert \le 2(1 + 5U)\lvert t_c\rvert_1`$. That is the cancellation bound.
   Altogether, $`\lvert t_c - \gamma_w s\rvert \le 50U\,(\lvert\mathrm{tr}\rvert + \lvert\mathrm{ti}\rvert) + 2^{-170}(\lvert\mathrm{sre}\rvert + \lvert\mathrm{sim}\rvert)(1 + 12U) + 8\eta_0`$.
   The last term covers underflow in the four products and two sums. The products need not be normal; each underflow
   error is $`\le \eta_0`$, which the relative argument above does not include. Also, the computed $`\lvert\mathrm{sre}\rvert + \lvert\mathrm{sim}\rvert`$ bounds
   $`\lvert s\rvert = \lvert 1/m^2\rvert`$ only up to a factor $`1 + 12U`$ (e.g. $`m = 3/4`$). The source comment at the charge says "$`\le`$" without this
   factor. Both points are covered by the existing slack, $`2^{-160} \gg 2^{-170}(1 + 12U) + 8\eta_0`$.

   The code charges $`(\mathsf{REL}(\lvert\mathrm{tr}\rvert + \lvert\mathrm{ti}\rvert) + 2^{-160}(1 + \lvert\mathrm{sre}\rvert + \lvert\mathrm{sim}\rvert))\cdot r \ge \lvert(t_c - \gamma_w s)\,\delta\rvert`$, with $`\mathsf{REL} \ge 900U`$.
   The extra factor $`(1 + \lvert\mathrm{sre}\rvert + \lvert\mathrm{sim}\rvert)`$ on the $`2^{-160}`$ term was added in v3 while this analysis was written: $`\lvert s\rvert`$ can be
   as large as $`2^{80}`$. The sums $`\mathrm{gre}[v] \mathrel{+}= \mathrm{tr}`$ and $`\mathrm{gim}[v] \mathrel{+}= \mathrm{ti}`$ add one rounding each. The own-variable term is the
   first and only term added to slot $`v`$, since `others` excludes $`v`$, so these roundings are exact.

4. **$`\lvert\gamma_w\rvert`$.** $`\mathrm{gw} = \mathrm{up}(\mathrm{cmod\_hi\_fast}(\mathrm{gr}, \mathrm{gi}) + \mathsf{REL}(\lvert\mathrm{gr}\rvert + \lvert\mathrm{gi}\rvert) + 2^{-160}) \ge \lvert\gamma_c\rvert + \lvert\gamma_c - \gamma_w\rvert \ge \lvert\gamma_w\rvert`$ (L2, item 1).
   There is exactly one coefficient with $`a = 1`$ and $`mk = 0`$, so gw is set once.

5. **$`\eta^{\ast}`$.** $`\mathrm{eta\_up} = \mathrm{up}(\mathrm{fl}(\mathrm{fl}(r\cdot r)/\mathrm{fl}(\mathrm{fl}(\mathrm{am}\cdot\mathrm{am})\cdot\mathrm{fl}(\mathrm{am} - r)))) \ge \eta^{\ast}`$: five roundings, all normal (W1), covered
   by (F6). The charge is $`\mathrm{up}(\mathrm{gw}\cdot\mathrm{eta\_up}) \ge \lvert\gamma_w\rvert\eta^{\ast}`$.

6. **Constant term.** $`a = \mathrm{la} = \mathrm{cr\_log}(\mathrm{fl}(\sqrt{\mathrm{c0n2}}))`$, with error $`\le 2.2U + U\lvert\mathrm{la}\rvert`$ (A.6, item 3). It is charged
   $`\mathsf{REL}(1 + \lvert\mathrm{la}\rvert)`$ in the final step.

7. **Remainder terms and accumulation.**
   - $`\mathrm{up}(\mathrm{rT}/\mathrm{A0lo})`$ and $`\mathrm{up}(q^2/(2(1-q))(1+4U))`$ are as in A.6, item 4.
   - ei is a sum of at most $`2 + 2(n-2) + 1`$ nonnegative terms.
   - The final $`\mathrm{up}(\mathrm{ei}\,(1 + 10^{-12}) + \mathsf{REL}(1 + \lvert\mathrm{la}\rvert))`$ dominates the exact sum (A.6, item 5).

So $`e \ge e^{\ast}`$ plus every computed-versus-exact discrepancy in $`a`$ and $`G`$ over $`B`$. ∎

**Finiteness (A2).** The model is declared invalid unless $`\mathrm{c0n2} \ge 2^{-900}`$ and $`a`$, $`e`$ and every gradient entry are finite.
Before the $`q`$ test, the w-form coefficients can be as large as about $`2^{210}`$ (W.3; L1 and L2 need only no overflow
and arguments $`\le 2^{500}`$). So $`\lvert\gamma\rvert`$ can reach about $`2^{660}`$, and tr, gw or ei can become $`\infty`$ or NaN. Every such case
fails `isfin` or the comparison $`q \lt  1/2`$, which is false on NaN, so no non-finite quantity can prune a box. When the
model is valid, $`q \lt  1/2`$ forces $`\lvert\gamma_\alpha\rvert\rho^\alpha \lt  1/2`$, in particular $`\lvert\gamma_w\rvert \lt  1/(2\,\mathrm{rw}) \le 2^{40}`$.

## W.6 Combination (end of `build_models`)
For $`i \ge 2`$, after v2's bounds and model are computed:
- $`\mathrm{hi}_i \leftarrow \min(\mathrm{hi}_i, \mathrm{hi}^w)`$ and $`\mathrm{lo}_i \leftarrow \max(\mathrm{lo}_i, \mathrm{lo}^w)`$;
- if the w-model is valid and v2's is invalid, or has a larger $`e`$, the w-model replaces it.

Every bound and every model is valid on its own: v2's by A.5/A.6, the w-model's by W.4/W.5. So the combination is
valid. v2's lower bound $`\mathrm{lo}_i`$ is valid at the nondegenerate points of $`B`$; the w-bound only exists when $`u_i \ne 0`$ on $`B`$.
The tests in `eval_box` (A.7) use only these quantities, so Theorem A holds for v3 unchanged.

## W.7 Mechanical cross-checks
1. `c/regress/models_v3_check.py`: every bound and model printed by the `models` mode is compared with critical values
   computed from first principles (`firstprinciples/fpcore.py`) at random points and corners of 4000 random domain boxes
   per degree. Half of the boxes are drawn where some $`\lvert u_j\rvert \in [0.25, 0.6]`$.

   | binary | $`d = 5`$ | $`d = 6`$ | $`d = 7`$ | violations | max err/e |
   |---|---|---|---|---|---|
   | v3 | 848 496 | 980 908 | 1 074 988 | 0 | 0.932 |
   | v2 | 731 724 | 793 456 | 824 780 | 0 | 0.932 |

   The counts are affine checks; v3 has 15–30 % more valid models.
2. The $`d = 5`$ v3 run is certificate-complete (`c/check_done.py`): 489 000 boxes, against 2 490 060 for v2.
   `firstprinciples/sample_tree.py` on its exported tree evaluated all 199 483 F/E/L leaves at 15.3 M points, with 0
   violations.

## W.9 A10 (run partition) and known limitations
A10 only decides which task ids a process computes. Records, the configuration and the checks are unchanged.
`c/merge_done.py` validates every record in full before removing duplicates, and `c/check_done.py` then requires every
task id exactly once.

Known limitation, from the review: the run program accepts malformed `SMALE_PART` strings, e.g. an empty residue list
or trailing characters. Such a run computes a different subset of tasks than intended. That cannot produce a false
completion, because the merged checkpoint must still contain every task id, which `check_done.py` verifies. The source
was not changed for this, since a change would alter the configuration hash of the running $`d = 7`$ computation.

`check_done.py` binds the records to the source through the configuration hash, and to one binary through
`--build`/`--bin`. A run split over machines has one build record per machine (`BUILD*.txt`); the source hash, which is
part of every record's configuration, is the common binding.

## W.8 Review status
See `AI_DISCLOSURE.md`.
