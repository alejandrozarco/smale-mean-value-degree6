# Appendix A — Floating-point error analysis of `c/smale_bb_v2.c`

Source: `c/smale_bb_v2.c`, sha256 `8ddcec1767e0b64b404843cc0e0d99684b6526ef9c6d34288129a89fe0b9487f`, built by
`c/build_v2.sh` (`-O2 -ffp-contract=off -fno-fast-math`, Apple clang 21, arm64). This appendix states and proves the
enclosure property of every numerical step used by the production mode `run`. Line numbers refer to that source.

**Main claim (Theorem A).** Let $`B = \prod_k [C_k - H_k, C_k + H_k]`$ be a box that satisfies the domain hypotheses (A5,
§A.2), and let $`d \le 7`$ and $`c = (d-1)/d`$. If `eval_box` returns
* 1 (F): then on $`B`$ there is an $`i`$ with $`\lvert S_i\rvert \lt  c`$;
* 2 (E): then there is a pair $`i \ne j`$ with $`\lvert S_i\rvert \ne \lvert S_j\rvert`$ at every point of $`B`$ where both are defined ($`u_i, u_j \ne 0`$), so
  $`B \cap E`$ contains no nondegenerate point;
* 4 (L): then $`(1/n) \sum_i \log\lvert S_i\rvert \lt  \log c`$ at every point of $`B`$.

If `geom_discard` returns 1, then $`\lvert u_v\rvert \gt  1`$ on $`B`$ for some $`v`$. If it returns 2, then on $`B`$ either $`\lvert u_v\rvert \lt  \lvert u_{v+1}\rvert`$ for
some $`v`$, or $`\mathrm{Im}\ u_2 \lt  0`$. If `excluded` returns 1 (mode 2), then $`B`$ lies in the closed Euclidean ball of radius $`r_{\mathrm{excl}}`$
around an exact equality point. Here "on $`B`$" means at every point of $`B`$ where the quantity is defined, i.e. $`u_i \ne 0`$ for
$`S_i`$ with $`i \ge 2`$. In particular the E-interval test may fire on a box that meets $`\{u_i = 0\}`$: $`\mathrm{lo}_i`$ stays finite there
(Lemma LB), and it bounds $`\log\lvert S_i\rvert`$ from below on the punctured box. This is exactly the quantification of the Prop R
corollary: its tests are over the nondegenerate points of the box. Every inequality holds with a positive margin that is uniform on $`B`$; this is what the corollary of
Proposition R (REDUCTION.md) uses at the degenerate boundary.

The mathematical content of the tests (Taylor models, the E-restriction, the uniform-weight L test) is described in
README.md and REDUCTION.md. This appendix covers only the passage from exact real arithmetic to the computed doubles.

---
## A.1 Trust base and floating-point model

**Assumptions** (outside this appendix):
1. The hardware implements IEEE-754 binary64 $`+, -, \times, \div, \sqrt{\cdot}`$ correctly rounded to nearest-even, with **gradual
   underflow** (no flush-to-zero), in the main thread and in every worker thread. Neither the program nor CORE-MATH
   changes the floating-point environment: there is no fesetenv / fesetround / FPCR write. The macOS arm64 default
   FPCR = 0 (round-to-nearest, FZ = 0, no default-NaN mode) was confirmed in the main thread and in three pthreads on
   the production machine (2026-09-29; a subnormal product rounded correctly in each). On arm64, FLT_EVAL_METHOD = 0,
   so there is no extended precision.
2. The compiler evaluates every expression as written: no contraction (`-ffp-contract=off`), no reassociation (no
   `-ffast-math`). `BUILD.txt` records the command and shows that every scalar fused multiply-add (`fmadd`/`fmsub`/`fnmadd`/`fnmsub`)
   in the binary lies inside `_cr_log`. The production binary (sha256 e99fdb88…) was also checked for vector
   `fmla`/`fmls` (none outside `_cr_log`). `build_v2.sh` now matches both kinds for future builds.
3. The only libm functions linked are `frexp` and `ldexp`. `fabs` and `floor` compile to exact instructions, and
   `sqrt` compiles to the IEEE `fsqrt` instruction: the production binary has 7 `fsqrt` and no call to libm
   sqrt/hypot/pow/log/fma. `fabs`, `frexp`, `floor` are exact. `ldexp(x, k)` is exact when the result is normal. When the result is subnormal, L3 needs
   only an error $`\lt  2^{-1074}`$ (faithful rounding); Apple's libm rounds correctly.
4. Compiler semantics used implicitly:
   - `isfin(x) = (x − x == 0)` is false exactly for $`\pm\infty`$ and NaN. This needs `-fno-fast-math` (in particular no
     finite-math assumption), which the build uses; checked under the build flags.
   - int → double conversions of small integers are exact, and literal expressions folded at compile time ($`1 + 8U`$,
     $`2(a+1)U`$, $`4U`$, $`(\mathrm{double})(d-1)/d`$, $`1.0/((k+1)(k+2))`$) are folded with IEEE round-to-nearest semantics.
5. The **compiled** CORE-MATH `cr_log` (third_party/core-math-log/log.c + dint.h, sha256 in BUILD.txt; upstream
   commit in UPSTREAM_COMMIT; authors P. Zimmermann and T. Hubrecht) is correctly rounded for every positive finite
   binary64 input. It is never called on 0, negative or NaN inputs (L5). $`+\infty`$ can occur only when $`\mathrm{supT} = \infty`$, and that
   case is excluded by the `supT < INFINITY` guard. This rests on the CORE-MATH authors' correctness argument for this implementation, and on its
   compilation being faithful: explicit `__builtin_fma` and integer operations are compiled as written, with
   `-ffp-contract=off` and no fast-math. We only tested it: 150k inputs against mpmath at
   200 bits, 0 misroundings. The analysis below needs much less than correct rounding: every log
   enters with an allowance $`\mathsf{REL}(1 + \lvert\log\rvert)`$, $`\mathsf{REL} \ge 900U`$, so a log that is wrong by up to ~50 ulp would still be
   covered.

**Notation.** $`U = 2^{-53}`$ (unit roundoff), $`\eta_0 = 2^{-1075}`$ (half the smallest subnormal). $`\mathrm{fl}(\cdot)`$ is the correctly rounded
result.

**(F1) Relative error.** If $`x`$ is real and $`\mathrm{fl}(x)`$ is normal, then $`\lvert \mathrm{fl}(x) - x\rvert \le U \lvert \mathrm{fl}(x)\rvert`$ and $`\lvert \mathrm{fl}(x) - x\rvert \le U\lvert x\rvert`$. If
$`\mathrm{fl}(x)`$ is subnormal or zero, then $`\lvert \mathrm{fl}(x) - x\rvert \le \eta_0`$. Hence, for $`x \ge 0`$ with $`\lvert x\rvert`$ below the overflow threshold (true for every operation in this analysis; see (M)):
$`\mathrm{fl}(x) \ge x(1-U) - \eta_0`$ and $`\mathrm{fl}(x) \le x(1+U) + \eta_0`$.

**(F2) Exact sums in the subnormal range.** If $`a, b`$ are doubles and $`\mathrm{fl}(a \pm b)`$ is subnormal, then $`a \pm b`$ is exact.
So the error term $`\eta_0`$ in (F1) is needed only for $`\times`$ and $`\div`$.

**(F3) Monotonicity.** $`x \le y`$ implies $`\mathrm{fl}(x) \le \mathrm{fl}(y)`$. Consequently, for a double $`y`$:
$`\mathrm{fl}(x) \lt  y \Rightarrow x \lt  y`$, and $`\mathrm{fl}(x) \gt  y \Rightarrow x \gt  y`$. (If $`x \ge y`$ then $`\mathrm{fl}(x) \ge \mathrm{fl}(y) = y`$.)

**(F4) Exact constants.** The source literals give exactly $`U = 2^{-53}`$ (`1.1102230246251565e-16`) and
$`\mathsf{INFL} = 1 + 2^{-40}`$ (`1.0 + 9.094947017729282e-13`); both were checked with exact rationals. The numbers
$`1 + kU`$ ($`k`$ even, $`k \le 18`$) and $`1 - kU`$ ($`k \le 18`$) are exact; the code uses $`1 + kU`$ only with even $`k`$. The products $`2U\cdot x`$ and
$`4U\cdot x`$ are exact whenever the product is normal. $`\mathsf{TINY} = 10^{-300} \gt  2^{-997} \gt  2^{77}\eta_0`$. $`\mathsf{REL} = 10^{-13}`$
$`\ge 900U`$. $`\mathsf{MARG} = 10^{-12} \ge 9000U`$. None of REL, MARG, TINY needs to be exact.

**(F5) Chains of nonnegative operations.** Let a nonnegative real $`y`$ be given by an expression tree whose leaves are
nonnegative doubles and whose inner nodes are $`+`$, $`\times`$, or $`\div`$ by a positive leaf. Let $`y_c`$ be the computed value, with each
node rounded. Define the rounding depth $`k`$ by: $`k(\text{leaf}) = 0`$, $`k(a + b) = \max(k(a), k(b)) + 1`$, $`k(a \times b) = k(a) + k(b) + 1`$,
$`k(a \div \text{leaf}) = k(a) + 1`$. Suppose no rounded $`\times`$ or $`\div`$ underflows, except possibly nodes whose result enters only further
additions. Then $`y_c \ge y(1-U)^{k} - N_u \eta_0`$, where $`N_u`$ is the number of such underflowing nodes.
*Proof:* induction with (F1), (F2), and monotonicity of $`+`$ and $`\times`$ on nonnegatives. At a product, the factors
$`(1-U)^{k(a)}`$ and $`(1-U)^{k(b)}`$ multiply. An underflow error that enters only additions is not amplified. ∎
The depth is at most the total number of operations. Every chain in the code has rounding depth $`\le 110`$, far below the
limit 4094 of (F6).

**(F6) `up`.** $`\mathrm{up}(x) = \mathrm{fl}(\mathrm{fl}(x\cdot\mathsf{INFL}) + \mathsf{TINY})`$. For $`x \ge 0`$, $`\mathrm{up}(x) \ge (x\cdot\mathsf{INFL}(1-U) - \eta_0 + \mathsf{TINY})(1-U)`$. Consequently: if
$`x \ge y(1-U)^{k} - N \eta_0`$ with $`k \le 4094`$ and $`N \le 2^{70}`$, then $`\mathrm{up}(x) \ge y`$. Reason: by Bernoulli, $`\mathsf{INFL}(1-U)^{k+2} \ge (1 + 2^{-40})(1 - (k+2)2^{-53}) = 1 + 2^{-40} - (k+2)2^{-53} - (k+2)2^{-93} \gt  1`$ for
$`k + 2 \le 4096`$, and $`\mathsf{TINY}(1-U) \gt  (N+1)\eta_0 \cdot \mathsf{INFL}`$. For two nested `up`: $`\mathsf{INFL}^2(1-U)^{k+4} \ge (1 + 2^{-39})(1 - (k+4)2^{-53}) \gt  1`$
for $`k + 4 \le 8192`$. (Every chain in this appendix has $`k \le 110`$.)

## A.2 Domain and magnitudes

**(A5) Domain.** Every box satisfies $`H_k = 2^{-e_k}`$ with $`0 \le e_k \le 40`$, $`C_k = m_k H_k`$ with $`m_k \in \mathbb{Z}`$, and
$`\lvert C_k\rvert + H_k \le 1`$. `in_domain` checks this for every top box with exact operations: frexp, and $`C/H`$ with $`H`$ a power of
two in the normal range, which is exact. Bisection halves $`H`$ (exact) and sets $`C \pm H/2`$ (exact, since it is a multiple
of $`2^{-41}`$ with modulus $`\le 1`$). It also asserts $`H \ge 2^{-40}`$ on every split. So:
* every $`C_k`$ and $`H_k`$ is a dyadic rational with at most 41 fractional bits, and $`\lvert C_k\rvert \le 1`$. Every centre
  $`m_v = C_{2v} + iC_{2v+1}`$ has $`\lvert m_v\rvert \le \sqrt{2}`$. The box $`B`$ is exactly the real set $`\prod[C_k - H_k, C_k + H_k]`$, with no rounding;
* $`\lvert C_k\rvert \pm H_k`$ is computed exactly: it is a multiple of $`2^{-41}`$ of modulus $`\le 2`$ (used in `modrange`). For out-of-domain
  input, `in_domain` still rejects: if $`\lvert c\rvert + h \gt  1`$ exactly, then by monotonicity $`\mathrm{fl}(\lvert c\rvert + h) \ge \mathrm{fl}`$ of the next double
  above $`1`$, which is $`\gt  1`$;
* $`\rho_v \ge H_{2v} \ge 2^{-40}`$. A monomial has at most $`n - 1 + (\mathrm{nv} - 1) \le 9`$ radius factors ($`d = 7`$), so every such
  product is $`\ge 2^{-360}`$ (no underflow);
* $`\mathrm{mre}^2 + \mathrm{mim}^2`$ is either $`0`$ or $`\ge 2^{-82}`$.

**(M) Magnitudes.** $`T_i`$ is a product of $`n - 1 \le 5`$ linear factors (the factor $`(1-t)`$ is absorbed into the moments),
integrated against moments $`\le 1/2`$.

*Majorant recurrence.* For a ball $`\beta = \langle z, r\rangle`$ put $`N(\beta) = \lvert \mathrm{Re}\ z\rvert + \lvert \mathrm{Im}\ z\rvert + r`$, and write $`\varepsilon = 2^{-38}`$.

**Claim (one operation).** As long as all quantities are $`\lt  2^{100}`$, so that nothing overflows:
- (i) $`N(\mathrm{cb\_mul}(a, b)) \le N(a)N(b)(1 + \varepsilon) + 2\cdot\mathsf{TINY}`$;
- (ii) $`N(\mathrm{cb\_addto}(q, a)) \le (N(q) + N(a))(1 + \varepsilon) + 2\cdot\mathsf{TINY}`$, and the same for cb_subto.

*Proof of (i).* By (F1), $`\lvert \mathrm{p1}\rvert \le \lvert a_{\mathrm{re}} b_{\mathrm{re}}\rvert(1+U) + \eta_0`$, and similarly for $`\mathrm{p2}, \mathrm{p3}, \mathrm{p4}`$. So
$`\lvert \mathrm{re}\rvert + \lvert \mathrm{im}\rvert \le (\lvert \mathrm{p1}\rvert + \lvert \mathrm{p2}\rvert + \lvert \mathrm{p3}\rvert + \lvert \mathrm{p4}\rvert)(1+U) \le \lVert a\rVert_1 \lVert b\rVert_1 (1+U)^2 + 5\eta_0`$, where $`\lVert x\rVert_1 = \lvert \mathrm{Re}\ x\rvert + \lvert \mathrm{Im}\ x\rvert`$.
The six partials sum to $`S \le 2\lVert a\rVert_1 \lVert b\rVert_1 (1+U)^2 + 10\eta_0`$. The pre-inflation radius is a nonnegative expression of depth 7,
so it is at most $`(\lVert a\rVert_1 r_b + r_a \lVert b\rVert_1 + r_a r_b + U\cdot S)(1+U)^{7} + 4\eta_0`$. The code returns $`\mathrm{fl}(\mathrm{fl}(\cdot\ \mathsf{INFL}) + \mathsf{TINY})`$, which is at
most that times $`\mathsf{INFL}(1+U)^2`$, plus $`\mathsf{TINY}(1+U)`$. Adding the midpoint bound, and using
$`\lVert a\rVert_1 \lVert b\rVert_1 + \lVert a\rVert_1 r_b + r_a \lVert b\rVert_1 + r_a r_b = N(a)N(b)`$, $`3U + (1+U)^{9}\cdot\mathsf{INFL} - 1 \lt  \varepsilon`$, and $`20\eta_0 \lt  \mathsf{TINY}`$, gives (i).

*Proof of (ii).* Same argument, with depth 3.

**Induction.** Let $`N^{\ast}`$ be the exact majorant: the same recurrence evaluated in exact arithmetic with every rounding
and TINY removed. Every computed ball on a chain of $`h`$ operations then satisfies
$`N \le N^{\ast}(1 + \varepsilon)^h + 2\cdot\mathsf{TINY}\cdot\sum_j \prod(\text{later factors}) \le N^{\ast}(1+\varepsilon)^h + 2^{-900}`$. The later factors are at most $`2^{13}`$ each and there are fewer than $`2^{13}`$ operations, so the
TINY terms total $`\lt  2^{13} \cdot 2^{13} \cdot 2^{65} \cdot 2\cdot\mathsf{TINY} \lt  2^{-900}`$ (the $`2^{65}`$ bounds products of at most 5 factors $`\le 2^{13}`$).

In `taylor_T`, every chain has $`h \le 5\cdot(1\text{ multiplication} + 4\text{ accumulations}) + 8 \lt  2^{6}`$. So $`(1+\varepsilon)^h \lt  1 + 2^{-31}`$.

$`N^{\ast}`$ is sub-multiplicative ($`\lVert zw\rVert_1 \le \lVert z\rVert_1 \lVert w\rVert_1`$). Each linear factor has $`N^{\ast} \le 2 + 2 + 1 + 1 = 6`$, because $`\lvert \mathrm{Re}\ m\rvert + \lvert \mathrm{Im}\ m\rvert \le 2`$.
So every intermediate slot has $`N^{\ast} \le 6^{n-1} \le 6^{5} \lt  2^{13}`$. After integration against moments $`\le 1/2`$, every
coefficient slot has $`N^{\ast} \le \tfrac12\cdot 6^{5} \lt  2^{12}`$.

**Hence** every intermediate ball has $`N \lt  2^{13}`$ and every Taylor ball satisfies $`\lvert \hat c_\alpha\rvert + r_\alpha \lt  2^{12}`$, for $`d \le 7`$. This
also justifies the $`2^{100}`$ premise used in the Claim.

*Derived quantities.* $`\rho_v \le \lvert (1, 1)\rvert(1 + 5U) + 2^{-1074} \lt  1.5`$. There are $`\le 9`$ radius factors and $`\le 96`$ summands, so
$`\mathrm{linmaj}, \mathrm{rT} \lt  96\cdot 2^{12}\cdot 1.5^{9}\cdot(1 + 2^{-38}) \lt  2^{26}`$. $`\mathrm{up}(\mathrm{up}(\cdot))`$ adds a factor $`\lt  1 + 2^{-38}`$. So $`\mathrm{supT} \lt  2^{27}`$, $`\mathrm{A0} \lt  2^{12}`$ and $`\lvert \mathrm{infT}\rvert \lt  2^{27}`$. The log arguments are $`\le 2^{27}`$ and $`\mathrm{phi} \le 2 + 1.5`$. So no overflow occurs in `taylor_T`, `cmod_*`,
`powu`, rT, linmaj, supT, infT, or the logs. Quantities that can be large are guarded:
* $`\gamma = c/\mathrm{c0}`$ is only formed when $`\lvert \mathrm{c0}\rvert^2 \ge 2^{-900}`$, so $`\lvert\gamma\rvert \lt  2^{14} \cdot 2^{450}`$;
* $`p`$ and $`q`$ may be large or $`+\infty`$ (division by a tiny or zero lower bound). The model is then marked invalid, because
  it requires $`p, q \lt  1/2`$;
* the final `isfin` check rejects any non-finite model quantity.

*Exceptional values: case analysis.* A strict comparison with $`\pm\infty`$ is not automatically false ($`\infty \gt  1`$ is true), so
safety follows from where the non-finite values can appear:
* $`\mathrm{hi}_i \in \mathbb{R} \cup \{+\infty\}`$. $`+\infty`$ arises when $`\mathrm{plo} \le 0`$ or $`\mathrm{supT} = \infty`$, and NaN is mapped to $`+\infty`$. hi appears only on the
  **smaller** side of a pruning comparison (F: $`\mathrm{hi} + \mathrm{margin} \lt  L_c`$; E: $`\ldots \gt  \mathrm{mhi} + \mathrm{margin}`$, with $`\mathrm{mhi} = \min_j \mathrm{hi}_j`$). $`+\infty`$
  there makes the comparison false.
* $`\mathrm{lo}_i \in \mathbb{R} \cup \{-\infty\}`$. $`-\infty`$ arises when $`\mathrm{infT} \le 0`$, and NaN is mapped to $`-\infty`$. lo appears only on the **larger** side
  (E: $`\mathrm{mlo} - \mathrm{margin} \gt  \ldots`$, with $`\mathrm{mlo} = \max_i \mathrm{lo}_i`$). $`-\infty`$ there makes it false. lo is never $`+\infty`$: $`\mathrm{li} = \log \mathrm{infT} \le \log \mathrm{A0lo} \lt  12\cdot\log 2`$.
* Affine-model quantities ($`a`$, $`G`$, $`e`$) enter tests only when the model is valid, and validity requires all of them to
  be finite (`isfin`).
* Geometric quantities (mn, mx, far2) are finite in the domain.
So a non-finite value can only prevent pruning, never cause it.

## A.3 Complex ball arithmetic (L1)

A ball $`\langle z, r\rangle = \{w \in \mathbb{C} : \lvert w - z\rvert \le r\}`$, with $`z = \mathrm{re} + i\cdot\mathrm{im}`$ and $`\mathrm{re}, \mathrm{im}, r`$ doubles. **Invariant:** every ball has either
$`r = 0`$ and an exact midpoint (inputs, the constants $`1, \pm m_v`$), or $`r \ge \mathsf{TINY}`$.

**Lemma L1a (`cb_mul`).** If $`\alpha \in \langle a, r_a\rangle`$ and $`\beta \in \langle b, r_b\rangle`$, then $`\alpha\beta \in`$ `cb_mul`$`(\langle a, r_a\rangle, \langle b, r_b\rangle)`$.
*Proof.* Let $`\mathrm{p1} = \mathrm{fl}(a_{\mathrm{re}} b_{\mathrm{re}})`$, $`\mathrm{p2} = \mathrm{fl}(a_{\mathrm{im}} b_{\mathrm{im}})`$, $`\mathrm{p3} = \mathrm{fl}(a_{\mathrm{re}} b_{\mathrm{im}})`$, $`\mathrm{p4} = \mathrm{fl}(a_{\mathrm{im}} b_{\mathrm{re}})`$, $`\mathrm{re} = \mathrm{fl}(\mathrm{p1} - \mathrm{p2})`$,
$`\mathrm{im} = \mathrm{fl}(\mathrm{p3} + \mathrm{p4})`$. By (F1)/(F2), $`\lvert \mathrm{re} - \mathrm{Re}(ab)\rvert \le U(\lvert \mathrm{p1}\rvert + \lvert \mathrm{p2}\rvert + \lvert \mathrm{re}\rvert) + 2\eta_0`$, and the same for im with $`\mathrm{p3}, \mathrm{p4}`$. Using
$`\lvert x + iy\rvert \le \lvert x\rvert + \lvert y\rvert`$, the midpoint error is at most $`U\cdot S + 4\eta_0`$, where $`S`$ is the sum of the six absolute values in the
code. Moreover $`\lvert \alpha\beta - ab\rvert \le \lvert a\rvert r_b + r_a \lvert b\rvert + r_a r_b`$, and $`\lvert a\rvert \le \lvert a_{\mathrm{re}}\rvert + \lvert a_{\mathrm{im}}\rvert`$. The required radius is therefore
$`R = \lVert a\rVert_1 r_b + r_a \lVert b\rVert_1 + r_a r_b + U\cdot S + 4\eta_0`$.

The computed pre-inflation radius is a nonnegative expression tree. Its rounding depth (F5) is 7: the three
products have depth $`\le 2`$, the six-term sum has depth 5 and $`U\times`$ adds 1, and the three outer additions (left to right)
give $`\max(\max(2, 2) + 1 + 1, 6) + 1 = 7`$. There are at most 4 underflowing products (three radius products, and
$`U\cdot\mathrm{sum}`$), none of them amplified. By (F5), $`\mathrm{rad}_c \ge (R - 4\eta_0)(1-U)^{7} - 4\eta_0`$. The code returns
$`\mathrm{fl}(\mathrm{fl}(\mathrm{rad}_c\cdot\mathsf{INFL}) + \mathsf{TINY}) = \mathrm{up}(\mathrm{rad}_c) \ge R`$ by (F6). ∎

**Lemma L1b (`cb_addto`, `cb_subto`).** If $`\alpha \in \langle q, r_q\rangle`$ and $`\beta \in \langle a, r_a\rangle`$, then $`\alpha \pm \beta`$ lies in the result. Midpoint
error: $`\le U(\lvert \mathrm{re}\rvert + \lvert \mathrm{im}\rvert)`$ by (F1), or 0 in the subnormal case (F2). Required: $`R = r_q + r_a + U(\lvert \mathrm{re}\rvert + \lvert \mathrm{im}\rvert)`$. The chain
has rounding depth 3 (two inner sums, $`U\times`$, outer sum) and one possible $`U\times`$ underflow, and $`\mathrm{up}(\cdot) \ge R`$ by (F6). ∎

**Lemma L1c (moments).** $`\mu_k = \int_0^1 (1-t)t^k\ dt = 1/((k+1)(k+2))`$. The denominator is an exact integer. The midpoint
$`\mathrm{fl}(1/\cdot)`$ has error $`\le U\cdot\mathrm{mid}`$, and the radius $`2U\cdot\mathrm{mid}`$ is exact and $`\ge U\cdot\mathrm{mid}`$. ∎

**Lemma T (`taylor_T`).** Write $`\delta_v = u_v - m_v`$. For $`i = 1`$ (code index 0),
$`S_1 = \int (1-t) \prod_v (1 - t m_v - t \delta_v)\ dt`$. For $`i \ge 2`$,

```math
T_i = \int (1-t)(m_i - t + \delta_i) \prod_{v\ne i}(m_i - t m_v + \delta_i - t \delta_v)\ dt.
```

The code expands these products as polynomials in $`(t, \delta_i, \delta_{\mathrm{others}})`$, with others entering multilinearly (bit masks).
The state $`A[t][a][\mathrm{mask}]`$ holds the coefficient of $`t^t \delta_i^a \prod_{q\in\mathrm{mask}} \delta_{\mathrm{others}[q]}`$. For $`i \ge 2`$ it starts from the factor
$`(m_i - t + \delta_i)`$:
- $`A[0][0][0] = m_i`$;
- $`A[1][0][0] = -1`$;
- $`A[0][1][0] = 1`$.

Each further factor $`(m_i - t m_v + \delta_i - t \delta_v)`$, where $`v = \mathrm{others}[q]`$ with $`\mathrm{bit} = 2^q`$, maps every state entry $`p`$ at
$`(t, a, \mathrm{mask})`$ to four contributions. Each contribution is a sum over all entries that map there:
- $`(t, a, \mathrm{mask})`$ receives $`p\cdot m_i`$;
- $`(t+1, a, \mathrm{mask})`$ receives $`p\cdot(-m_v)`$;
- $`(t, a+1, \mathrm{mask})`$ receives $`p`$;
- $`(t+1, a, \mathrm{mask}\vert\mathrm{bit})`$ receives $`-p`$.

For $`i = 1`$ the factor is $`(1 - t m_v - t \delta_v)`$, which gives the rules $`(t, \mathrm{mask}) \mathrel{+}= p`$, $`(t+1, \mathrm{mask}) \mathrel{+}= p\cdot(-m_v)`$ and
$`(t+1, \mathrm{mask}\vert\mathrm{bit}) \mathrel{-}= p`$. Finally $`\mathrm{coef}[a][\mathrm{mask}] = \sum_t A[t][a][\mathrm{mask}]\cdot\mu_t`$, with $`\mu_t = \int_0^1 (1-t)t^t\ dt`$. Expanding the
products term by term gives exactly these rules, so $`\mathrm{coef}[a][\mathrm{mask}] = c_{a,\mathrm{mask}}`$.

The products with the exact constant ($`m_i`$ or $`-m_v`$) are formed with `cb_mul`. The $`\delta`$-terms shift a ball to another slot,
possibly negated, and all contributions to a slot are accumulated with `cb_addto` / `cb_subto`. These are rounded
additions, covered by L1b. The result is then integrated with the moment balls. Every elementary step is an exact
polynomial operation whose ball version satisfies L1a–L1c. So each output `coef[a][mask]` contains the exact Taylor
coefficient $`c_{a,\mathrm{mask}}`$ of $`T_i`$ (resp. $`S_1`$) at the centre $`m`$. Skipping exact-zero balls (`p == 0, rad == 0`) is exact.
The algebra was also cross-validated against an independent Python implementation and python-flint. ∎

## A.4 Moduli, powers, logarithms (L2–L5)

**L2 (`cmod_hi_fast`, used for $`\lvert c\rvert`$ of Taylor coefficients; $`\lvert x\rvert, \lvert y\rvert \le 2^{14} \ll 2^{500}`$).**
Claim: $`\sqrt{x^2 + y^2} \le \mathrm{fl}(\mathrm{fl}(\mathrm{fl}(\sqrt{s})\cdot(1+8U)) + 2^{-530})`$, where $`s = \mathrm{fl}(\mathrm{fl}(x^2) + \mathrm{fl}(y^2))`$.
*Proof.* By (F1), $`s \ge (x^2 + y^2)(1-U)^2 - 2\eta_0(1-U) \ge (x^2 + y^2)(1-U)^2 - 2^{-1073}`$. Then
$`\sqrt{s} \ge \lvert z\rvert(1-U) - 2^{-536.5}`$ (using $`\sqrt{A - B} \ge \sqrt{A} - \sqrt{B}`$), and $`r = \mathrm{fl}(\sqrt{s}) \ge \sqrt{s}(1-U)`$. Next,
$`\mathrm{fl}(r(1+8U)) \ge r(1+8U)(1-U) - \eta_0`$, and the final sum is $`\ge (\text{that} + 2^{-530})(1-U)`$. Collecting terms:
$`\ge \lvert z\rvert(1-U)^{4}(1+8U) - 2^{-536.5}(1+8U) - \eta_0 + 2^{-530}(1-U) \ge \lvert z\rvert`$, because $`(1-U)^{4}(1+8U) \ge 1`$ and
$`2^{-530}(1-U) \gt  2^{-536.4} + \eta_0`$. ∎

**L3 (`cmod_bounds`, finite $`x, y`$ with $`\max(\lvert x\rvert, \lvert y\rvert) \le 2^{1000}`$).** Claim: $`\mathrm{lo} \le \lvert x + iy\rvert \le \mathrm{hi}`$. (Not for arbitrary finite
inputs: near the overflow threshold $`\mathrm{ldexp}(r(1-8U), e)`$ can overflow, so $`\mathrm{lo} = +\infty`$. Every caller passes values $`\lt  2^{13}`$:
box coordinates, half-widths, the Taylor constant term (M), and $`\lvert C\rvert \pm H \le 2`$.)
*Proof.* Let $`M = \max(\lvert x\rvert, \lvert y\rvert) \gt  0`$ and $`M = f\cdot 2^e`$ with $`f \in [1/2, 1)`$ (frexp, exact). Then $`\mathrm{xs} = x\cdot 2^{-e}`$, $`\mathrm{ys} = y\cdot 2^{-e}`$ are
exact, except that the smaller one may underflow (error $`\le \eta_0`$). The larger one lies in $`[1/2, 1)`$, so
$`\mathrm{xs}^2 + \mathrm{ys}^2 \in [1/4, 2]`$. Two products and one sum give $`s = (\mathrm{xs}^2 + \mathrm{ys}^2)(1 + \theta_1)`$ with $`\lvert\theta_1\rvert \le 2U + 2U^2 + 12\eta_0`$. Then
$`r = \mathrm{fl}(\sqrt{s}) = \lvert z\rvert 2^{-e}(1 + \theta)`$, $`\lvert\theta\rvert \le U + \lvert\theta_1\rvert/2 + U^2 \le 2.5U`$.

hi: $`\mathrm{fl}(r(1+8U)) \ge r(1+8U)(1-U) \ge \lvert z\rvert 2^{-e}(1 + 4U)`$. If $`\mathrm{ldexp}(\cdot, e)`$ is normal it is exact. If it is subnormal, it is
rounded down by $`\lt  2^{-1074}`$ ($`\le \eta_0`$ if correctly rounded), and then $`+2^{-1074}`$ restores it, exactly by (F2). If the result is normal, $`\mathrm{fl}(y + 2^{-1074}) \ge y`$ by
(F3), and $`y \ge \lvert z\rvert(1 + 4U) \gt  \lvert z\rvert`$.
lo: symmetric with $`(1 - 8U)`$, $`-2^{-1074}`$, and $`\max(0, \cdot)`$. ∎

**L4 (`powu`, $`x \ge 0`$, $`0 \le a \le 8`$).** Claim: $`x^a \le \mathrm{fl}(\mathrm{fl}(x\cdots x)\cdot(1 + 2(a+1)U)) + \mathsf{TINY}`$.
*Proof.* For $`a = 0`$ (called with $`a = 0`$ in `build_models`), $`r_c = 1`$ and the result is $`\mathrm{fl}(\mathrm{fl}(1\cdot(1+2U)) + \mathsf{TINY})`$
$`= 1 + 2U \ge 1 = x^0`$. For $`a \ge 1`$: in every call $`x = \rho \ge 2^{-40}`$ (the source header says the weaker $`x \ge 2^{-60}`$), so the
repeated products do not underflow ($`x^8 \ge 2^{-320}`$). The chain has $`a - 1`$ roundings (the first product $`1\cdot x`$ is exact), so
by (F5) $`r_c \ge x^a(1-U)^{a-1} \ge x^a(1-U)^a`$. The
constant $`1 + 2(a+1)U`$ is exact. Then $`\mathrm{fl}(r_c(1 + 2(a+1)U)) \ge x^a(1-U)^{a+1}(1 + 2(a+1)U) \ge x^a`$, because
$`(1-U)^{a+1} \ge 1 - (a+1)U`$ and $`(1 - kU)(1 + 2kU) \ge 1`$ for $`kU \le 1/2`$. The final addition $`\mathrm{fl}(y + \mathsf{TINY}) \ge y`$ by (F3),
since $`y`$ is a double. ∎

**L5 (`cr_log`).** Let $`x`$ be positive and finite, and $`\ell = \mathrm{cr\_log}(x)`$. Then $`\lvert \ell - \log x\rvert \le \tfrac12 \mathrm{ulp}(\log x) \le U\lvert\log x\rvert`$,
and so $`\lvert \ell - \log x\rvert \le U\lvert\ell\rvert(1 + 2U)`$. ($`\log x`$ is never subnormal for a double $`x`$: the smallest nonzero $`\lvert\log x\rvert`$ is
$`\approx 2^{-53}`$.) Every argument passed to cr_log in the code is positive and finite:
* $`\mathrm{supT} \ge \mathsf{TINY}`$ (up);
* $`\mathrm{infT} \gt  0`$ and $`\mathrm{plo} \gt  0`$ are tested;
* $`\mathrm{phi} \gt  0`$;
* $`\sqrt{\mathrm{c0n2}} \ge 2^{-450}`$;
* $`\mathrm{absm} \ge 2^{-41}`$;
* $`\mathrm{fl}((d-1)/d) \gt  0`$.

## A.5 Enclosures for $`\lvert S_i\rvert`$ (`build_models`, interval part)

Fix $`i`$, and write $`\delta \in B - m`$. So $`\delta_v`$ ranges over the rectangle $`[-H_{2v}, H_{2v}] \times [-H_{2v+1}, H_{2v+1}]`$, which lies
in the disc of radius $`\lvert (H_{2v}, H_{2v+1})\rvert \le \rho_v`$ (L3). Let $`c_\alpha`$ be the exact coefficients, with $`\alpha`$ a multi-index of
degree $`\lvert\alpha\rvert`$, and $`\langle \hat c_\alpha, r_\alpha\rangle`$ the computed balls (Lemma T).

**Lemma R (remainder and linear majorants).** Define
$`\mathrm{rT}^{\ast} = r_0 + \sum_{\lvert\alpha\rvert=1} r_\alpha \rho^\alpha + \sum_{\lvert\alpha\rvert\ge 2} (\lvert \hat c_\alpha\rvert + r_\alpha) \rho^\alpha`$ and $`\mathrm{linmaj}^{\ast} = \sum_{\lvert\alpha\rvert=1} \lvert \hat c_\alpha\rvert \rho^\alpha`$.
Then the computed $`\mathrm{rT} \ge \mathrm{rT}^{\ast}`$ and $`\mathrm{linmaj} \ge \mathrm{linmaj}^{\ast}`$, and for all $`\delta`$,

```math
\Bigl\lvert T(m+\delta) - \hat c_0 - \sum_{\lvert\alpha\rvert=1} \hat c_\alpha \delta^\alpha\Bigr\rvert \le \mathrm{rT}^{\ast}.
```

*Proof.* The inequality follows from $`\lvert\delta^\alpha\rvert \le \rho^\alpha`$ and $`\lvert c_\alpha - \hat c_\alpha\rvert \le r_\alpha`$.

Computation of the terms:
* $`\mathrm{ra} = \mathrm{powu}(\rho_i, a) \ge \rho_i^a`$ (L4).
* rp multiplies in $`\le \mathrm{nb} \le 5`$ factors $`\rho`$ (no underflow, A.2), with a chain of $`\le 5`$ roundings.
* $`\mathrm{ac} \ge \lvert \hat c_\alpha\rvert`$ (L2), and $`\mathrm{ac} \ge 2^{-530}`$, so $`\mathrm{ac}\cdot\mathrm{rp}`$ and $`(\mathrm{ac} + r)\cdot\mathrm{rp}`$ do not underflow.
* $`r_\alpha\cdot\mathrm{rp}`$ may underflow, but it enters only additions.
* For $`d \le 7`$, rT has at most $`N \le 1 + (6\cdot 2^{4} - 1) = 96`$ summands ($`d = 6`$: 40). The source comment "≤ 40 terms" is the
  $`d = 6`$ count.

Each summand has rounding depth $`\le 7`$ ($`\mathrm{rp} \le 5`$, times ac or $`\mathrm{ac} + r \le 1 + 1`$). Recursive summation adds $`\le 96`$ more. By (F5),
$`\mathrm{rT}_c \ge \mathrm{rT}^{\ast}(1-U)^{103} - 96\eta_0`$, and $`\mathrm{up}(\mathrm{up}(\cdot)) \ge \mathrm{rT}^{\ast}`$ by (F6). The same argument applies to linmaj ($`\le 6`$ summands). ∎

**Lemma S (supT, infT).** Let $`\mathrm{A0lo} \le \lvert \hat c_0\rvert \le \mathrm{A0hi}`$ (L3). Then
$`\mathrm{supT} = \mathrm{up}(\mathrm{A0hi} + \mathrm{linmaj} + \mathrm{rT}) \ge \mathrm{A0hi} + \mathrm{linmaj} + \mathrm{rT} \ge \sup_B \lvert T\rvert`$, and $`\mathrm{infT} \le \mathrm{A0lo} - \mathrm{linmaj} - \mathrm{rT} \le \inf_B \lvert T\rvert`$.

*Proof of the infT bound.* Write $`\Sigma = \mathrm{A0lo} + \mathrm{linmaj} + \mathrm{rT}`$ (reals, all $`\ge 0`$) and $`X = \mathrm{A0lo} - \mathrm{linmaj} - \mathrm{rT}`$. The code
computes $`\mathrm{x2} = \mathrm{fl}(\mathrm{fl}(\mathrm{A0lo} - \mathrm{linmaj}) - \mathrm{rT})`$ with $`\mathrm{x2} \le X + 2U\Sigma + U^2\Sigma`$ (F1). It also computes
$`y = \mathrm{fl}(4U\cdot\mathrm{fl}(\mathrm{fl}(\mathrm{A0lo} + \mathrm{linmaj}) + \mathrm{rT}))`$. Then $`z = \mathrm{fl}(\mathrm{x2} - y)`$ and $`w = \mathrm{fl}(z - \mathsf{TINY})`$.

($`\Sigma \ge \mathrm{rT} \ge \mathsf{TINY} \gt  0`$ always: rT includes the coefficient radius $`r_0 \ge \mathsf{TINY}`$ after up.)
Case 1: $`y`$ is normal. Then $`y \ge 4U\Sigma(1-U)^2`$ and $`\lvert z\rvert \le \Sigma(1 + 7U)`$, so
$`z \le \mathrm{x2} - y + U\lvert z\rvert \le X - U\Sigma(1 - 17U) \lt  X`$, and $`w \le z`$ by (F3).
Case 2: $`y`$ is subnormal (so $`\Sigma \lt  2^{-970}`$). Then $`\mathrm{x2} \le X + 2^{-1021}`$ and $`z \le \mathrm{x2} + U\lvert z\rvert \le X + 2^{-1020}`$. So
$`w \le (z - \mathsf{TINY}) + U\lvert z - \mathsf{TINY}\rvert \lt  X`$, since $`\mathsf{TINY} \gt  2^{-997}`$.
In both cases $`\mathrm{infT} = w \le X`$. ∎

**Lemma LB (log bounds, `hi`/`lo`).** For $`i \ge 2`$ (code $`i \gt  0`$), $`\lvert u_i\rvert \in [\lvert m_i\rvert - \rho_i, \lvert m_i\rvert + \rho_i]`$ on $`B`$.
* $`\mathrm{plo} = \mathrm{fl}(\mathrm{fl}(\mathrm{absm}_{\mathrm{lo}} - \rho)\cdot(1 - 4U)) \le (\mathrm{absm}_{\mathrm{lo}} - \rho)(1+U)^2(1 - 4U) \lt  \mathrm{absm}_{\mathrm{lo}} - \rho \le \min_B \lvert u_i\rvert`$ when
  $`t := \mathrm{absm}_{\mathrm{lo}} - \rho \gt  0`$. The difference $`t`$ is never subnormal: $`\mathrm{absm}_{\mathrm{lo}}`$ and $`\rho`$ are doubles $`\ge 2^{-41}`$ (A5), hence multiples of
  $`2^{-93}`$, so $`t \le 0`$ or $`t \ge 2^{-93}`$.
* Symmetrically, $`\mathrm{phi} \ge \max_B \lvert u_i\rvert`$.

Exact bounds: $`\log\lvert S_i\rvert \le \log \mathrm{supT} - (n-1) \log \mathrm{plo}`$ and $`\log\lvert S_i\rvert \ge \log \mathrm{infT} - (n-1) \log \mathrm{phi}`$.

Computed: $`\mathrm{ls} = \mathrm{cr\_log}(\mathrm{supT})`$ and $`\mathrm{dlo} = \mathrm{fl}((n-1)\cdot\mathrm{cr\_log}(\mathrm{plo}))`$. The expression
$`\mathrm{hi} = \mathrm{fl}(\mathrm{fl}(\mathrm{ls} - \mathrm{dlo}) + \mathrm{fl}(\mathsf{REL}\cdot\mathrm{fl}(\mathrm{fl}(1 + \lvert \mathrm{ls}\rvert) + \mathrm{dmag})))`$ differs from $`\mathrm{ls} - \mathrm{dlo}`$ by the following roundings:
* at most $`U\lvert \mathrm{ls}\rvert`$ from the log (L5),
* $`(n-1)U\lvert\log \mathrm{plo}\rvert + U\lvert \mathrm{dlo}\rvert`$ from dlo,
* $`U(\lvert \mathrm{ls}\rvert + \lvert \mathrm{dlo}\rvert)`$ from the subtraction.

The total is $`\le 10U(1 + \lvert \mathrm{ls}\rvert + \mathrm{dmag})`$, where $`\mathrm{dmag} \ge \lvert \mathrm{dlo}\rvert(1-3U)`$. The added allowance (the code evaluates
$`(1 + \lvert \mathrm{ls}\rvert) + \mathrm{dmag}`$) is $`\mathsf{REL}(1 + \lvert \mathrm{ls}\rvert + \mathrm{dmag})(1-U)^{3} \ge 890U(1 + \lvert \mathrm{ls}\rvert + \mathrm{dmag})`$. The final addition can lower the result by
at most $`U\lvert \mathrm{hi}\rvert \le U(\lvert \mathrm{ls}\rvert + \mathrm{dmag} + 1)(1 + O(U))`$. So $`\mathrm{hi} \ge \log \mathrm{supT} - (n-1) \log \mathrm{plo} \ge \sup_B \log\lvert S_i\rvert`$. For $`i = 1`$ there is
no $`u`$-term. The lo bound is symmetric.

Non-finite cases, as coded: $`\mathrm{hi} = +\infty`$ if $`\mathrm{plo} \le 0`$ or $`\mathrm{supT} = \infty`$; $`\mathrm{lo} = -\infty`$ if $`\mathrm{infT} \le 0`$. **lo can remain finite when $`\mathrm{plo} \le 0`$**
($`u_i`$ may vanish on $`B`$), **provided $`\mathrm{infT} \gt  0`$**: then $`\mathrm{lo}_i = \log \mathrm{infT} - (n-1) \log \mathrm{phi} - \mathsf{REL}(1 + \lvert \mathrm{li}\rvert + \mathrm{dmag})`$ with $`\mathrm{dmag} = (n-1)\lvert\log \mathrm{phi}\rvert`$.
This is a valid lower bound for $`\log\lvert S_i\rvert`$ at every point of $`B`$ with $`u_i \ne 0`$ (there $`\lvert S_i\rvert = \lvert T_i\rvert/\lvert u_i\rvert^{n-1}`$
$`\ge \mathrm{infT}/\mathrm{phi}^{n-1}`$; $`\lvert S_i\rvert \to \infty`$ as $`u_i \to 0`$). So an E-interval verdict that uses $`\mathrm{lo}_i`$ on such a box holds on the
nondegenerate points, as stated in Theorem A. ∎

## A.6 Affine models of $`\log\lvert S_i\rvert`$ (`build_models`, affine part)

Assume the model is marked valid. That requires:
* grad_ok ($`\mathrm{c0n2} := \mathrm{fl}(\hat c_{0\mathrm{re}}^2 + \hat c_{0\mathrm{im}}^2) \ge 2^{-900}`$);
* $`\mathrm{A0lo} \gt  0`$ and $`q \lt  1/2`$;
* for $`i \ge 2`$: $`p \lt  1/2`$ and $`\mathrm{m2} \ge 2^{-900}`$;
* all quantities finite.

**Exact model.** On $`B`$, $`T(m+\delta) = \hat c_0(1 + w)`$ with $`w = \sum_{\lvert\alpha\rvert=1} \gamma_\alpha \delta^\alpha + E/\hat c_0`$, where $`\gamma_\alpha = \hat c_\alpha/\hat c_0`$ and $`\lvert E\rvert \le \mathrm{rT}^{\ast}`$ (Lemma
R). So $`\lvert w\rvert \le (\mathrm{linmaj} + \mathrm{rT})/\lvert \hat c_0\rvert \le (\mathrm{linmaj} + \mathrm{rT})/\mathrm{A0lo} \le q \lt  1/2`$. Taking the principal branch,
$`\log\lvert T\rvert = \log\lvert \hat c_0\rvert + \mathrm{Re}(\gamma\cdot\delta) + \mathrm{Re}(E/\hat c_0) + \mathrm{Re}(\log(1+w) - w)`$, with
$`\lvert \mathrm{Re}(E/\hat c_0)\rvert \le \mathrm{rT}/\mathrm{A0lo}`$ and $`\lvert \log(1+w) - w\rvert \le \sum_{k\ge 2}\lvert w\rvert^k/k \le q^2/(2(1-q))`$ ($`x \mapsto x^2/(2(1-x))`$ is increasing on $`[0,1)`$).

Likewise, for $`i \ge 2`$, $`\log\lvert u_i\rvert = \log\lvert m\rvert + \mathrm{Re}(\delta/m) + R'`$ with $`\lvert R'\rvert \le p^2/(2(1-p))`$, where $`p \ge \rho/\lvert m\rvert`$. With
$`\mathrm{Re}(\gamma\delta) = \gamma_{\mathrm{re}} x - \gamma_{\mathrm{im}} y`$ and $`\mathrm{Re}(\delta/m) = (m_{\mathrm{re}} x + m_{\mathrm{im}} y)/\lvert m\rvert^2`$, the code's gradient row ($`G[2v] = \mathrm{gre}`$, $`G[2v+1] = -\mathrm{gim}`$)
is exactly the real gradient of $`\mathrm{Re}(\gamma\cdot\delta) - (n-1)\mathrm{Re}(\delta/m)`$. Hence the exact statement:
$`\log\lvert S_i(m+\delta)\rvert = \hat a + G^{\ast}\cdot(x,y) + \varepsilon`$ with $`\lvert\varepsilon\rvert \le e^{\ast}`$, where $`\hat a = \log\lvert \hat c_0\rvert - (n-1)\log\lvert m\rvert`$, $`G^{\ast}`$ is the exact gradient, and

```math
e^{\ast} = \mathrm{rT}/\mathrm{A0lo} + q^2/(2(1-q)) + (n-1)p^2/(2(1-p)).
```

**Lemma A (computed model is valid).** For all $`\delta \in B - m`$: $`\bigl\lvert \log\lvert S_i(m+\delta)\rvert - a - G\cdot(x,y)\bigr\rvert \le e`$, where $`(a, G, e)`$ are the
computed quantities. The error sources, all charged to $`e`$:

1. **$`\gamma`$.** $`N = \hat c_{\mathrm{re}}\ \mathrm{c0re} + \hat c_{\mathrm{im}}\ \mathrm{c0im}`$. By Cauchy–Schwarz, $`\lvert \hat c_{\mathrm{re}}\ \mathrm{c0re}\rvert + \lvert \hat c_{\mathrm{im}}\ \mathrm{c0im}\rvert \le \lvert \hat c\rvert \lvert \hat c_0\rvert`$, so
   $`\lvert \mathrm{fl}(N) - N\rvert \le 2.1U\lvert \hat c\rvert \lvert \hat c_0\rvert + 2\eta_0`$. Also $`\mathrm{c0n2} = \lvert \hat c_0\rvert^2(1 + \theta)`$ with $`\lvert\theta\rvert \le 2.1U`$: component squares may underflow, but
   the error $`2\eta_0 \ll U\cdot 2^{-900}`$. Hence $`\lvert \mathrm{gr}_c - \mathrm{Re}\ \gamma\rvert \le 4.5U\lvert\gamma\rvert + U\lvert \mathrm{gr}_c\rvert + 2^{-172}`$, and the same holds for gi. Since
   $`\lvert\gamma\rvert \le (\lvert \mathrm{gr}_c\rvert + \lvert \mathrm{gi}_c\rvert)(1 + 12U) + 2^{-171}`$, $`\lvert \gamma_c - \gamma\rvert \le 11U(\lvert \mathrm{gr}_c\rvert + \lvert \mathrm{gi}_c\rvert) + 2^{-170}`$. The code charges
   $`(\mathsf{REL}(\lvert \mathrm{gr}_c\rvert + \lvert \mathrm{gi}_c\rvert) + 2^{-160})\cdot\rho_{\mathrm{var}} \ge \lvert (\gamma_c - \gamma)\delta\rvert`$ ($`\mathsf{REL} \ge 900U`$). For $`i \ge 2`$, adding the $`u`$-gradient into $`\mathrm{gre}[v]`$
   costs one more rounding of $`\le U(\lvert \mathrm{gr}_T\rvert + \lvert \mathrm{gr}_u\rvert)`$, which is covered by the two REL charges.
2. **$`u`$-gradient.** $`\mathrm{gr} = \mathrm{fl}(\mathrm{fl}(-(n-1)m_{\mathrm{re}})/\mathrm{m2})`$ has relative error $`\le 4.2U`$ ($`\mathrm{m2}`$ relative error $`\le 2.1U`$); the same holds
   for gi. It is charged $`\mathsf{REL}(\lvert \mathrm{gr}\rvert + \lvert \mathrm{gi}\rvert)\rho_v + 2^{-160}\rho_v`$.
3. **Constant term.** $`\sqrt{\mathrm{c0n2}} = \lvert \hat c_0\rvert(1 + \theta')`$ with $`\lvert\theta'\rvert \le 2.1U`$, so $`\lvert \mathrm{la} - \log\lvert \hat c_0\rvert\rvert \le 2.2U + U\lvert \mathrm{la}\rvert`$ (L5). Likewise
   $`\lvert \mathrm{lm} - \log\lvert m\rvert\rvert \le 2.2U + U\lvert \mathrm{lm}\rvert`$. The operations $`\mathrm{ai} \mathrel{-}= \mathrm{fl}((n-1)\cdot\mathrm{lm})`$ add $`\le U((n-1)\lvert \mathrm{lm}\rvert + \lvert \mathrm{ai}\rvert)`$. The total is
   $`\le 20U(1 + \mathrm{amag} + \lvert \mathrm{ai}\rvert)`$, charged $`\mathsf{REL}(1 + \mathrm{amag} + \lvert \mathrm{ai}\rvert)`$.
4. **Remainder terms.** $`\mathrm{up}(\mathrm{rT}/\mathrm{A0lo}) \ge \mathrm{rT}/\mathrm{A0lo}`$ by (F5)+(F6). Now take $`x = q^2/(2(1-q))`$, computed with roundings
   $`q\cdot q`$, $`1 - q`$ (up $`\le U`$, which lowers the value), $`2\times`$ (exact), $`\div`$, and $`\times(1+4U)`$. $`q`$ is only known to be $`\gtrsim 2^{-584}`$, so $`q\cdot q`$,
   the $`\div`$ and the $`\times`$ may underflow. Each underflow error $`\eta_0`$ is then multiplied by at most $`(1+4U)/(2(1-q)) \le 1 + 4U`$. So
   $`x_c \ge x(1-U)^{4} - 3\eta_0(1+4U)`$, and $`\mathrm{up}(\cdot) \ge x`$ by (F6) with $`k = 4`$, $`N = 3`$. The term $`(n-1)p^2/(2(1-p))`$ has one more
   rounding ($`(n-1)\cdot p`$): $`x_c \ge x(1-U)^{5} - 4\eta_0(1+4U)`$, and again $`\mathrm{up}(\cdot) \ge x`$. (The factor $`(1+4U)`$ in the code is not needed
   for this argument.) $`q \ge (\mathrm{linmaj} + \mathrm{rT})/\mathrm{A0lo}`$
   and $`p \ge \rho/\lvert m\rvert`$ hold because each is $`\mathrm{up}(\cdot)`$ of a quotient whose numerator is an upper bound and whose denominator is
   a lower bound ($`\mathrm{absm}_{\mathrm{lo}}(1-4U)`$ rounded is $`\lt  \mathrm{absm}_{\mathrm{lo}} \le \lvert m\rvert`$).
5. **Accumulation.** ei is a sum of $`\le 2 + 2\mathrm{nv}`$ nonnegative terms ($`\le 12`$ roundings). The final
   $`\mathrm{up}(\mathrm{ei}\cdot(1 + 10^{-12}) + \ldots)`$ makes the result $`\ge`$ the exact sum, since $`10^{-12} \gg 12U`$.

Therefore $`e \ge e^{\ast} +`$ (all computed-vs-exact discrepancies in $`a`$ and $`G`$, integrated over $`B`$), which is the claim. ∎

## A.7 Tests (`eval_box`)

Let $`L_c =`$ `logc_lo` $`= \mathrm{fl}(\mathrm{cr\_log}(\mathrm{fl}((d-1)/d)) - 10^{-14})`$. $`\mathrm{fl}((d-1)/d)`$ has relative error $`\le U`$, and $`\lvert\log c\rvert \le \log(4/3)`$
$`\lt  0.29`$ for $`d \ge 4`$, so $`\lvert \mathrm{cr\_log}(\mathrm{fl}(c)) - \log c\rvert \le 1.3U`$. So $`L_c \le \log c`$.

**F-interval.** $`\mathrm{fl}(\mathrm{hi} + \mathrm{fl}(\mathsf{MARG}(1 + \lvert \mathrm{hi}\rvert))) \lt  L_c \Rightarrow \mathrm{hi} \lt  L_c`$ by (F3), since the margin is $`\ge 0`$. Then
$`\log\lvert S_i\rvert \lt  \log c`$ on $`B`$ (Lemma LB).

**F-affine.** Let $`s^{\ast} = a + \sum_k\lvert G_k\rvert H_k + e`$ be the exact real value with the computed $`a, G, e`$. Then $`\sup_B \log\lvert S_i\rvert \le s^{\ast}`$
(Lemma A; $`\lvert G\cdot(x,y)\rvert \le \sum\lvert G_k\rvert H_k`$ on the rectangle). Also $`s_c`$ is a $`(D+2)`$-term recursive sum of products, so
$`\lvert s_c - s^{\ast}\rvert \le (D+3)U\cdot\mathrm{mag}^{\ast}(1 + U) + D\eta_0`$, where $`\mathrm{mag}^{\ast}`$ is the same sum with $`\lvert a\rvert`$, $`D + 3 \le 13`$, and $`D\eta_0`$ covers products
$`\lvert G_k\rvert\cdot H_k`$ that underflow. By (F3), $`\mathrm{fl}(s_c + \mathsf{MARG}\cdot(1 + \mathrm{mag}_c)) \lt  L_c`$ gives $`s_c + \mathsf{MARG}(1 + \mathrm{mag}_c)(1 - 2U) \lt  L_c`$. So
$`s^{\ast} \le s_c + 14U\cdot\mathrm{mag}_c + D\eta_0 \lt  L_c`$, since $`\mathsf{MARG}\cdot 1`$ absorbs $`D\eta_0`$.

**E-interval.** $`\mathrm{fl}(\mathrm{mlo} - m_1) \gt  \mathrm{fl}(\mathrm{mhi} + m_2)`$ with $`m_1, m_2 \ge 0`$ implies $`\mathrm{mlo} \gt  \mathrm{mhi}`$ (F3). So for some $`i, j`$,
$`\inf \log\lvert S_i\rvert \ge \mathrm{lo}_i \gt  \mathrm{hi}_j \ge \sup \log\lvert S_j\rvert`$ on $`B`$. Moreover $`i \ne j`$. $`B`$ has interior points with all $`u_v \ne 0`$ ($`H_k \gt  0`$), and
at such a point every finite $`\mathrm{lo}_i`$ satisfies $`\mathrm{lo}_i \le \log\lvert S_i\rvert \le \mathrm{hi}_i`$ (Lemma LB; if $`\mathrm{lo}_i`$ is finite then $`\mathrm{infT} \gt  0`$, so
$`T_i \ne 0`$). Hence $`\mathrm{lo}_i \le \mathrm{hi}_i`$ for every $`i`$, and $`\mathrm{mlo} \gt  \mathrm{mhi}`$ forces the index of the maximum lo to differ from the index of
the minimum hi.

**E-affine.** Let $`\mathrm{lhs}^{\ast} = \lvert a_i - a_j\rvert`$ and $`\mathrm{rhs}^{\ast} = e_i + e_j + \sum_k\lvert G_{ik} - G_{jk}\rvert H_k`$. The rounding errors satisfy
$`\lvert \mathrm{lhs}_c - \mathrm{lhs}^{\ast}\rvert \le U\cdot\mathrm{lhs}_c`$ and $`\lvert \mathrm{rhs}_c - \mathrm{rhs}^{\ast}\rvert \le (D+4)U\cdot\mathrm{rhs}_c(1+U) + D\eta_0`$ ($`D\eta_0`$ for underflowing products
$`\lvert G_{ik} - G_{jk}\rvert\cdot H_k`$). The test is
$`\mathrm{fl}(\mathrm{lhs}_c - \mathrm{rhs}_c) \gt  \mathsf{MARG}(1 + \lvert a_i\rvert + \lvert a_j\rvert + \mathrm{rhs}_c)`$ (computed). It implies $`\mathrm{lhs}^{\ast} - \mathrm{rhs}^{\ast} \gt  0`$, because the margin exceeds
$`2U\cdot\mathrm{lhs}_c + 15U\cdot\mathrm{rhs}_c + U\lvert \mathrm{lhs}_c - \mathrm{rhs}_c\rvert + D\eta_0`$. Then on $`B`$, $`\lvert l_i - l_j\rvert \ge \mathrm{lhs}^{\ast} - \mathrm{rhs}^{\ast} \gt  0`$, i.e. $`\lvert S_i\rvert \ne \lvert S_j\rvert`$.

**L (uniform weights).** On $`B \cap E`$ all $`\log\lvert S_i\rvert`$ are equal, so $`(1/n)\sum_i \log\lvert S_i\rvert = \log F`$. Let
$`s^{\ast} = \sum_i (a_i + e_i) + \sum_k \lvert \sum_i G_{ik}\rvert H_k`$. Then $`\sum_i \log\lvert S_i\rvert \le s^{\ast}`$ on $`B`$. The inner sums $`g_k = \sum_i G_{ik}`$ carry error
$`\le nU\cdot\sum_i\lvert G_{ik}\rvert`$, and the outer sum has $`\le 2n + D`$ terms. The products $`\lvert g_k\rvert\cdot H_k`$ and $`\mathrm{ga}\cdot H_k`$ can underflow ($`\le 2D\eta_0`$ in
total). So $`\lvert s_c - s^{\ast}\rvert \le (3n + D + 2)U\cdot\mathrm{mag}_c(1 + U) + 2D\eta_0 \le 31U\cdot\mathrm{mag}_c + 2D\eta_0`$ for $`d \le 7`$ ($`3n + D + 2 = 30`$ at $`d = 7`$). The right-hand side is
$`Y_c = \mathrm{fl}(\mathrm{fl}(n\cdot L_c) - \mathsf{MARG}) \le nL_c + 2.1U\cdot n\lvert L_c\rvert - \mathsf{MARG}(1 - U)`$. The test $`\mathrm{fl}(s_c + \mathsf{MARG}(1 + \mathrm{mag}_c)) \lt  Y_c`$ gives, by (F3),
$`s_c + \mathsf{MARG}(1 + \mathrm{mag}_c)(1 - 2U) \lt  Y_c`$, hence

```math
s^{\ast} \lt  nL_c + 2.1U\cdot n\lvert L_c\rvert - \mathsf{MARG}(1 - U) - \mathsf{MARG}(1 + \mathrm{mag}_c)(1 - 2U) + 31U\cdot\mathrm{mag}_c + 2D\eta_0 \lt  nL_c \le n \log c.
```

So $`\log F \lt  \log c`$ on $`B \cap E`$. More precisely, $`(1/n)\sum \log\lvert S_i\rvert \lt  \log c`$ on
all of $`B`$, which is the uniform strict form the Prop R corollary needs.

**Uniformity.** Each implication above gives a strict inequality whose gap is bounded below on $`B`$. Indeed $`B`$ is
compact, and each bound is a fixed real computed from $`(a, G, e, \mathrm{hi}, \mathrm{lo})`$ with slack $`\ge \mathsf{MARG}`$. This is the "uniform
strict test" of the Prop R corollary.

## A.8 Geometry and exclusion

**Outside / symmetry (`modrange`, `geom_discard`).** In the domain, $`\mathrm{dx} = \max(0, \lvert C_x\rvert - H_x)`$ and the analogous $`\mathrm{dy}`$,
$`\mathrm{ex}`$, $`\mathrm{ey}`$ are exact (A.2), and $`\min_B \lvert u_v\rvert = \lvert (\mathrm{dx}, \mathrm{dy})\rvert`$, $`\max_B \lvert u_v\rvert = \lvert (\mathrm{ex}, \mathrm{ey})\rvert`$. By L3, $`\mathrm{mn}_c \le \min \lvert u_v\rvert`$ and
$`\mathrm{mx}_c \ge \max \lvert u_v\rvert`$. Multiplying by $`(1 \mp 4U)`$ and rounding moves each further in the safe direction (F3).
* `mn > 1 + 10^−12` $`\Rightarrow \lvert u_v\rvert \gt  1`$ on $`B`$.
* `fl(mx(1 + 10^−12)) < fl(mn′(1 − 10^−12))` $`\Rightarrow \mathrm{mx} \lt  \mathrm{mn}'`$ (F3) $`\Rightarrow \lvert u_v\rvert \lt  \lvert u_{v+1}\rvert`$ on $`B`$.
(Here and below, $`1 \pm 10^{-12}`$, $`10^{-13}`$ denote the doubles the literals round to; only $`1 + 10^{-12} \gt  1 \gt  1 - 10^{-12}`$ and
$`10^{-13} \gt  0.99\cdot 10^{-13}`$ are used.)
* `C[1] + H[1] < 0` is exact $`\Rightarrow \mathrm{Im}\ u_2 \lt  0`$ on $`B`$.

**Exclusion, mode 2 (`excluded`).**
Let $`\hat p`$ be the computed equality point. `c/regress/eqpoints_v2.py` **certifies** $`\lVert \hat p - p\rVert^2 \lt  (1.1\cdot 10^{-16})^2`$ in arb
(python-flint ball arithmetic, 300 bits; exact points via arb cos_pi/sin_pi; the doubles enter exactly) for every
output point, $`d = 4..7`$ (runs/regress_v2/eqpoints_v2.txt: ALL PASS). It also compares with mpmath in both directions,
which is provenance only, since soundness needs only the computed → exact direction. Each computed point is within $`1.1\cdot 10^{-16}`$ (Euclidean) of an exact
equality point, each exact point has a computed point that close, and the counts are $`(n-1)!`$. The maximum Euclidean
errors are $`d = 4`$: $`7.1\cdot 10^{-17}`$, $`d = 5`$: $`0`$, $`d = 6`$: $`8.0\cdot 10^{-17}`$, $`d = 7`$: $`1.0\cdot 10^{-16}`$. We use $`\lVert \hat p - p\rVert \le 1.6\cdot 10^{-16}`$. (Even the
weaker per-coordinate bound $`10^{-15}`$ asserted by c/regress/regress_v2.py, i.e. $`\le 3.2\cdot 10^{-15}`$ in norm, would leave
$`9.6\cdot 10^{-14}`$ of slack below.)

For $`u \in B`$, $`\lvert u - \hat p\rvert_2 \le (\sum_k (\lvert C_k - \hat p_k\rvert + H_k)^2)^{1/2} =: \Phi`$. Each $`\mathrm{fx}_c = \mathrm{fl}(\mathrm{fl}(\lvert \mathrm{fl}(C - \hat p)\rvert) + H)`$
$`\ge (\lvert C - \hat p\rvert + H)(1-U)^2 - \eta_0`$. far2 is a sum of $`2\mathrm{nv} \le 10`$ squares, so $`\mathrm{far2}_c \ge \Phi^2(1-U)^{15} - 20\eta_0`$, and
$`\mathrm{fl}(\sqrt{\mathrm{far2}_c}) \ge \Phi(1 - 8.5U) - 2^{-535}`$. The test is $`\mathrm{fl}(Y) \le r_{\mathrm{excl}}`$ with $`Y = \mathrm{fl}(\sqrt{\mathrm{far2}_c}\cdot(1 + 10^{-12})) + 10^{-13}`$, the last
sum taken exactly. If $`Y \gt  r_{\mathrm{excl}} + \tfrac12\mathrm{ulp}(r_{\mathrm{excl}})`$, then $`\mathrm{fl}(Y) \gt  r_{\mathrm{excl}}`$. So the test implies $`Y \le r_{\mathrm{excl}}(1 + U)`$, i.e.
$`\Phi(1 - 8.5U)(1 + 10^{-12})(1 - U) - 2^{-534} + 10^{-13} \le r_{\mathrm{excl}}(1 + U)`$. Since $`(1 - 8.5U)(1 + 10^{-12})(1 - U) \gt  1`$, we get
$`\Phi \le r_{\mathrm{excl}}(1 + U) - 10^{-13} + 2^{-534}`$. Therefore $`\sup_B \lvert u - p\rvert \le \Phi + 1.6\cdot 10^{-16} \lt  r_{\mathrm{excl}} - 9\cdot 10^{-14} \lt  r_{\mathrm{excl}}`$, and $`B`$ lies
in the closed ball. Here $`r_{\mathrm{excl}}`$ is the double parsed from "0.05", which is $`1/20 + 1/360287970189639680`$. Since
$`\sup_B \lvert u - p\rvert \lt  r_{\mathrm{excl}} - 9\cdot 10^{-14} \lt  0.05`$, $`B`$ also lies in the closed ball of exact radius 0.05 used by
local/LOCAL_CERT.md and the Prop R corollary.

## A.9 What is not covered here

* Non-floating-point control flow of the driver. One known limitation: in `run`, `long total = s^D` would overflow for
  s = 2^20, D = 10, far outside the configurations used (s = 4, D = 6 or 8). Grid coverage is re-checked independently
  in exact arithmetic by c/check_done.py.

* The algebraic identities S_i = T_i/u_i^{n−1} and the chart and symmetry reduction: README.md.
* Proposition R and its corollary: REDUCTION.md (two audits).
* The local certificate on the excluded balls: local/LOCAL_CERT.md (arb).
* Completeness of the cover, i.e. every grid box is a task or exactly discardable, all tasks done, 0 unresolved:
  c/check_done.py (exact-rational grid recomputation).
* Correct rounding of CORE-MATH `cr_log`: its authors' proof (assumption 4).

## A.10 Mechanical cross-checks of this appendix

These are tests, not proofs. They are listed so a reader can see which lemmas were also exercised numerically.
* L1a/b, L2–L5: c/regress/prim_test_v2.c + prim_check_v2.py (adversarial inputs in all binades, exact rationals /
  mpmath).
* A.8 equality points: c/regress/eqpoints_v2.py (both directions, d = 4..7, mpmath 300 bits): ALL PASS, max 1.004e−16.

