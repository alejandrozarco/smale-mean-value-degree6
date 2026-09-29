# Appendix A — Floating-point error analysis of `c/smale_bb_v2.c`

Source: `c/smale_bb_v2.c`, sha256 `8ddcec1767e0b64b404843cc0e0d99684b6526ef9c6d34288129a89fe0b9487f`, built by
`c/build_v2.sh` (`-O2 -ffp-contract=off -fno-fast-math`, Apple clang 21, arm64). This appendix states and proves the
enclosure property of every numerical step used by the production mode `run`. Line numbers refer to that source.

**Main claim (Theorem A).** Let B = Π_k [C_k − H_k, C_k + H_k] be a box that satisfies the domain hypotheses (A5,
§A.2), and let d ≤ 7 and c = (d−1)/d. If `eval_box` returns
* 1 (F): then on B there is an i with |S_i| < c;
* 2 (E): then there is a pair i ≠ j with |S_i| ≠ |S_j| at every point of B where both are defined (u_i, u_j ≠ 0), so
  B ∩ E contains no nondegenerate point;
* 4 (L): then (1/n) Σ_i log|S_i| < log c at every point of B.

If `geom_discard` returns 1, then |u_v| > 1 on B for some v. If it returns 2, then on B either |u_v| < |u_{v+1}| for
some v, or Im u_2 < 0. If `excluded` returns 1 (mode 2), then B lies in the closed Euclidean ball of radius r_excl
around an exact equality point. Here "on B" means at every point of B where the quantity is defined, i.e. u_i ≠ 0 for
S_i with i ≥ 2. In particular the E-interval test may fire on a box that meets {u_i = 0}: lo_i stays finite there
(Lemma LB), and it bounds log|S_i| from below on the punctured box. This is exactly the quantification of the Prop R
corollary: its tests are over the nondegenerate points of the box. Every inequality holds with a positive margin that is uniform on B; this is what the corollary of
Proposition R (REDUCTION.md) uses at the degenerate boundary.

The mathematical content of the tests (Taylor models, the E-restriction, the uniform-weight L test) is described in
README.md and REDUCTION.md. This appendix covers only the passage from exact real arithmetic to the computed doubles.

---
## A.1 Trust base and floating-point model

**Assumptions** (outside this appendix):
1. The hardware implements IEEE-754 binary64 +, −, ×, ÷, √ correctly rounded to nearest-even, with **gradual
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
   only an error < 2^−1074 (faithful rounding); Apple's libm rounds correctly.
4. Compiler semantics used implicitly:
   - `isfin(x) = (x − x == 0)` is false exactly for ±∞ and NaN. This needs `-fno-fast-math` (in particular no
     finite-math assumption), which the build uses; checked under the build flags.
   - int → double conversions of small integers are exact, and literal expressions folded at compile time (1 + 8U,
     2(a+1)U, 4U, (double)(d−1)/d, 1.0/((k+1)(k+2))) are folded with IEEE round-to-nearest semantics.
5. The **compiled** CORE-MATH `cr_log` (third_party/core-math-log/log.c + dint.h, sha256 in BUILD.txt; upstream
   commit in UPSTREAM_COMMIT; authors P. Zimmermann and T. Hubrecht) is correctly rounded for every positive finite
   binary64 input. It is never called on 0, negative or NaN inputs (L5). +∞ can occur only when supT = ∞, and that
   case is excluded by the `supT < INFINITY` guard. This rests on the CORE-MATH authors' correctness argument for this implementation, and on its
   compilation being faithful: explicit `__builtin_fma` and integer operations are compiled as written, with
   `-ffp-contract=off` and no fast-math. We only tested it: 150k inputs against mpmath at
   200 bits, 0 misroundings. The analysis below needs much less than correct rounding: every log
   enters with an allowance REL(1 + |log|), REL ≥ 900U, so a log that is wrong by up to ~50 ulp would still be
   covered.

**Notation.** U = 2^−53 (unit roundoff), η₀ = 2^−1075 (half the smallest subnormal). fl(·) is the correctly rounded
result.

**(F1) Relative error.** If x is real and fl(x) is normal, then |fl(x) − x| ≤ U|fl(x)| and |fl(x) − x| ≤ U|x|. If
fl(x) is subnormal or zero, then |fl(x) − x| ≤ η₀. Hence, for x ≥ 0 with |x| below the overflow threshold (true for every operation in this analysis; see (M)):
fl(x) ≥ x(1−U) − η₀ and fl(x) ≤ x(1+U) + η₀.

**(F2) Exact sums in the subnormal range.** If a, b are doubles and fl(a ± b) is subnormal, then a ± b is exact.
So the error term η₀ in (F1) is needed only for × and ÷.

**(F3) Monotonicity.** x ≤ y implies fl(x) ≤ fl(y). Consequently, for a double y:
fl(x) < y ⇒ x < y, and fl(x) > y ⇒ x > y. (If x ≥ y then fl(x) ≥ fl(y) = y.)

**(F4) Exact constants.** The source literals give exactly U = 2^−53 (`1.1102230246251565e-16`) and
INFL = 1 + 2^−40 (`1.0 + 9.094947017729282e-13`); both were checked with exact rationals. The numbers
1 + kU (k even, k ≤ 18) and 1 − kU (k ≤ 18) are exact; the code uses 1 + kU only with even k. The products 2U·x and
4U·x are exact whenever the product is normal. TINY = 10^−300 > 2^−997 > 2^77 η₀. REL = 10^−13
≥ 900U. MARG = 10^−12 ≥ 9000U. None of REL, MARG, TINY needs to be exact.

**(F5) Chains of nonnegative operations.** Let a nonnegative real y be given by an expression tree whose leaves are
nonnegative doubles and whose inner nodes are +, ×, or ÷ by a positive leaf. Let y_c be the computed value, with each
node rounded. Define the rounding depth k by: k(leaf) = 0, k(a + b) = max(k(a), k(b)) + 1, k(a × b) = k(a) + k(b) + 1,
k(a ÷ leaf) = k(a) + 1. Suppose no rounded × or ÷ underflows, except possibly nodes whose result enters only further
additions. Then y_c ≥ y(1−U)^k − N_u η₀, where N_u is the number of such underflowing nodes.
*Proof:* induction with (F1), (F2), and monotonicity of + and × on nonnegatives. At a product, the factors
(1−U)^{k(a)} and (1−U)^{k(b)} multiply. An underflow error that enters only additions is not amplified. ∎
The depth is at most the total number of operations. Every chain in the code has rounding depth ≤ 110, far below the
limit 4094 of (F6).

**(F6) `up`.** up(x) = fl(fl(x·INFL) + TINY). For x ≥ 0, up(x) ≥ (x·INFL(1−U) − η₀ + TINY)(1−U). Consequently: if
x ≥ y(1−U)^k − N η₀ with k ≤ 4094 and N ≤ 2^70, then up(x) ≥ y. Reason: by Bernoulli, INFL(1−U)^{k+2} ≥ (1 + 2^−40)(1 − (k+2)2^−53) = 1 + 2^−40 − (k+2)2^−53 − (k+2)2^−93 > 1 for
k + 2 ≤ 4096, and TINY(1−U) > (N+1)η₀ · INFL. For two nested `up`: INFL²(1−U)^{k+4} ≥ (1 + 2^−39)(1 − (k+4)2^−53) > 1
for k + 4 ≤ 8192. (Every chain in this appendix has k ≤ 110.)

## A.2 Domain and magnitudes

**(A5) Domain.** Every box satisfies H_k = 2^{−e_k} with 0 ≤ e_k ≤ 40, C_k = m_k H_k with m_k ∈ Z, and
|C_k| + H_k ≤ 1. `in_domain` checks this for every top box with exact operations: frexp, and C/H with H a power of
two in the normal range, which is exact. Bisection halves H (exact) and sets C ± H/2 (exact, since it is a multiple
of 2^−41 with modulus ≤ 1). It also asserts H ≥ 2^−40 on every split. So:
* every C_k and H_k is a dyadic rational with at most 41 fractional bits, and |C_k| ≤ 1. Every centre m_v = C_{2v} +
  iC_{2v+1} has |m_v| ≤ √2. The box B is exactly the real set Π[C_k − H_k, C_k + H_k], with no rounding;
* |C_k| ± H_k is computed exactly: it is a multiple of 2^−41 of modulus ≤ 2 (used in `modrange`). For out-of-domain
  input, `in_domain` still rejects: if |c| + h > 1 exactly, then by monotonicity fl(|c| + h) ≥ fl of the next double
  above 1, which is > 1;
* ρ_v ≥ H_{2v} ≥ 2^−40. A monomial has at most n − 1 + (nv − 1) ≤ 9 radius factors (d = 7), so every such
  product is ≥ 2^−360 (no underflow);
* mre² + mim² is either 0 or ≥ 2^−82.

**(M) Magnitudes.** T_i is a product of n − 1 ≤ 5 linear factors (the factor (1−t) is absorbed into the moments),
integrated against moments ≤ 1/2.

*Majorant recurrence.* For a ball β = ⟨z, r⟩ put N(β) = |Re z| + |Im z| + r, and write ε = 2^−38.

**Claim (one operation).** As long as all quantities are < 2^100, so that nothing overflows:
- (i) N(cb_mul(a, b)) ≤ N(a)N(b)(1 + ε) + 2·TINY;
- (ii) N(cb_addto(q, a)) ≤ (N(q) + N(a))(1 + ε) + 2·TINY, and the same for cb_subto.

*Proof of (i).* By (F1), |p1| ≤ |a_re b_re|(1+U) + η₀, and similarly for p2, p3, p4. So
|re| + |im| ≤ (|p1| + |p2| + |p3| + |p4|)(1+U) ≤ ‖a‖₁‖b‖₁(1+U)² + 5η₀, where ‖x‖₁ = |Re x| + |Im x|.
The six partials sum to S ≤ 2‖a‖₁‖b‖₁(1+U)² + 10η₀. The pre-inflation radius is a nonnegative expression of depth 7,
so it is at most (‖a‖₁r_b + r_a‖b‖₁ + r_ar_b + U·S)(1+U)^7 + 4η₀. The code returns fl(fl(·INFL) + TINY), which is at
most that times INFL(1+U)², plus TINY(1+U). Adding the midpoint bound, and using ‖a‖₁‖b‖₁ + ‖a‖₁r_b + r_a‖b‖₁ + r_ar_b =
N(a)N(b), 3U + (1+U)^9·INFL − 1 < ε, and 20η₀ < TINY, gives (i).

*Proof of (ii).* Same argument, with depth 3.

**Induction.** Let N* be the exact majorant: the same recurrence evaluated in exact arithmetic with every rounding
and TINY removed. Every computed ball on a chain of h operations then satisfies N ≤ N*(1 + ε)^h + 2·TINY·Σ_j Π(later
factors) ≤ N*(1+ε)^h + 2^−950. The later factors are at most 2^13 each and there are fewer than 2^13 operations, so the
TINY terms total < 2^13 · 2^13 · 2^65 · 2·TINY < 2^−950 (the 2^65 bounds products of at most 5 factors ≤ 2^13).

In `taylor_T`, every chain has h ≤ 5·(1 multiplication + 4 accumulations) + 8 < 2^6. So (1+ε)^h < 1 + 2^−31.

N* is sub-multiplicative (‖zw‖₁ ≤ ‖z‖₁‖w‖₁). Each linear factor has N* ≤ 2 + 2 + 1 + 1 = 6, because |Re m| + |Im m| ≤ 2.
So every intermediate slot has N* ≤ 6^{n−1} ≤ 6^5 < 2^13. After integration against moments ≤ 1/2, every
coefficient slot has N* ≤ ½·6^5 < 2^12.

**Hence** every intermediate ball has N < 2^13 and every Taylor ball satisfies |ĉ_α| + r_α < 2^12, for d ≤ 7. This
also justifies the 2^100 premise used in the Claim.

*Derived quantities.* ρ_v ≤ |(1, 1)|(1 + 5U) + 2^−1074 < 1.5. There are ≤ 9 radius factors and ≤ 96 summands, so
linmaj, rT < 96·2^12·1.5^9·(1 + 2^−38) < 2^26. up(up(·)) adds a factor < 1 + 2^−38. So supT < 2^27, A0 < 2^12 and |infT| < 2^27. The log arguments are ≤ 2^27 and phi ≤ 2 + 1.5. So no overflow occurs in `taylor_T`, `cmod_*`,
`powu`, rT, linmaj, supT, infT, or the logs. Quantities that can be large are guarded:
* γ = c/c0 is only formed when |c0|² ≥ 2^−900, so |γ| < 2^14 · 2^450;
* p and q may be large or +∞ (division by a tiny or zero lower bound). The model is then marked invalid, because
  it requires p, q < 1/2;
* the final `isfin` check rejects any non-finite model quantity.

*Exceptional values: case analysis.* A strict comparison with ±∞ is not automatically false (∞ > 1 is true), so
safety follows from where the non-finite values can appear:
* hi_i ∈ R ∪ {+∞}. +∞ arises when plo ≤ 0 or supT = ∞, and NaN is mapped to +∞. hi appears only on the
  **smaller** side of a pruning comparison (F: hi + margin < L_c; E: … > mhi + margin, with mhi = min_j hi_j). +∞
  there makes the comparison false.
* lo_i ∈ R ∪ {−∞}. −∞ arises when infT ≤ 0, and NaN is mapped to −∞. lo appears only on the **larger** side
  (E: mlo − margin > …, with mlo = max_i lo_i). −∞ there makes it false. lo is never +∞: li = log infT ≤ log A0lo < 12·log 2.
* Affine-model quantities (a, G, e) enter tests only when the model is valid, and validity requires all of them to
  be finite (`isfin`).
* Geometric quantities (mn, mx, far2) are finite in the domain.
So a non-finite value can only prevent pruning, never cause it.

## A.3 Complex ball arithmetic (L1)

A ball ⟨z, r⟩ = {w ∈ C : |w − z| ≤ r}, with z = re + i·im and re, im, r doubles. **Invariant:** every ball has either
r = 0 and an exact midpoint (inputs, the constants 1, ±m_v), or r ≥ TINY.

**Lemma L1a (`cb_mul`).** If α ∈ ⟨a, r_a⟩ and β ∈ ⟨b, r_b⟩, then αβ ∈ `cb_mul`(⟨a,r_a⟩, ⟨b,r_b⟩).
*Proof.* Let p1 = fl(a_re b_re), p2 = fl(a_im b_im), p3 = fl(a_re b_im), p4 = fl(a_im b_re), re = fl(p1 − p2),
im = fl(p3 + p4). By (F1)/(F2), |re − Re(ab)| ≤ U(|p1| + |p2| + |re|) + 2η₀, and the same for im with p3, p4. Using
|x + iy| ≤ |x| + |y|, the midpoint error is at most U·S + 4η₀, where S is the sum of the six absolute values in the
code. Moreover |αβ − ab| ≤ |a| r_b + r_a |b| + r_a r_b, and |a| ≤ |a_re| + |a_im|. The required radius is therefore
R = ‖a‖₁ r_b + r_a ‖b‖₁ + r_a r_b + U·S + 4η₀.

The computed pre-inflation radius is a nonnegative expression tree. Its rounding depth (F5) is 7: the three
products have depth ≤ 2, the six-term sum has depth 5 and U× adds 1, and the three outer additions (left to right)
give max(max(2, 2) + 1 + 1, 6) + 1 = 7. There are at most 4 underflowing products (three radius products, and
U·sum), none of them amplified. By (F5), rad_c ≥ (R − 4η₀)(1−U)^7 − 4η₀. The code returns
fl(fl(rad_c·INFL) + TINY) = up(rad_c) ≥ R by (F6). ∎

**Lemma L1b (`cb_addto`, `cb_subto`).** If α ∈ ⟨q, r_q⟩ and β ∈ ⟨a, r_a⟩, then α ± β lies in the result. Midpoint
error: ≤ U(|re| + |im|) by (F1), or 0 in the subnormal case (F2). Required: R = r_q + r_a + U(|re| + |im|). The chain
has rounding depth 3 (two inner sums, U×, outer sum) and one possible U× underflow, and up(·) ≥ R by (F6). ∎

**Lemma L1c (moments).** μ_k = ∫₀¹(1−t)t^k dt = 1/((k+1)(k+2)). The denominator is an exact integer. The midpoint
fl(1/·) has error ≤ U·mid, and the radius 2U·mid is exact and ≥ U·mid. ∎

**Lemma T (`taylor_T`).** Write δ_v = u_v − m_v. For i = 1 (code index 0),
S_1 = ∫(1−t) Π_v (1 − t m_v − t δ_v) dt. For i ≥ 2,
T_i = ∫(1−t)(m_i − t + δ_i) Π_{v≠i}(m_i − t m_v + δ_i − t δ_v) dt.
The code expands these products as polynomials in (t, δ_i, δ_others), with others entering multilinearly (bit masks).
The state A[t][a][mask] holds the coefficient of t^t δ_i^a Π_{q∈mask} δ_{others[q]}. For i ≥ 2 it starts from the factor
(m_i − t + δ_i):
- A[0][0][0] = m_i;
- A[1][0][0] = −1;
- A[0][1][0] = 1.

Each further factor (m_i − t m_v + δ_i − t δ_v), where v = others[q] with bit = 2^q, maps every state entry p at
(t, a, mask) to four contributions. Each contribution is a sum over all entries that map there:
- (t, a, mask) receives p·m_i;
- (t+1, a, mask) receives p·(−m_v);
- (t, a+1, mask) receives p;
- (t+1, a, mask|bit) receives −p.

For i = 1 the factor is (1 − t m_v − t δ_v), which gives the rules (t, mask) += p, (t+1, mask) += p·(−m_v) and
(t+1, mask|bit) −= p. Finally coef[a][mask] = Σ_t A[t][a][mask]·μ_t, with μ_t = ∫₀¹(1−t)t^t dt. Expanding the
products term by term gives exactly these rules, so coef[a][mask] = c_{a,mask}.

The products with the exact constant (m_i or −m_v) are formed with `cb_mul`. The δ-terms shift a ball to another slot,
possibly negated, and all contributions to a slot are accumulated with `cb_addto` / `cb_subto`. These are rounded
additions, covered by L1b. The result is then integrated with the moment balls. Every elementary step is an exact
polynomial operation whose ball version satisfies L1a–L1c. So each output `coef[a][mask]` contains the exact Taylor
coefficient c_{a,mask} of T_i (resp. S_1) at the centre m. Skipping exact-zero balls (`p == 0, rad == 0`) is exact.
The algebra was also cross-validated against an independent Python implementation and python-flint. ∎

## A.4 Moduli, powers, logarithms (L2–L5)

**L2 (`cmod_hi_fast`, used for |c| of Taylor coefficients; |x|, |y| ≤ 2^14 ≪ 2^500).**
Claim: sqrt(x² + y²) ≤ fl(fl(fl(√s)·(1+8U)) + 2^−530), where s = fl(fl(x²) + fl(y²)).
*Proof.* By (F1), s ≥ (x² + y²)(1−U)² − 2η₀(1−U) ≥ (x² + y²)(1−U)² − 2^−1073. Then
√s ≥ |z|(1−U) − 2^−536.5 (using √(A − B) ≥ √A − √B), and r = fl(√s) ≥ √s(1−U). Next,
fl(r(1+8U)) ≥ r(1+8U)(1−U) − η₀, and the final sum is ≥ (that + 2^−530)(1−U). Collecting terms:
≥ |z|(1−U)^4(1+8U) − 2^−536.5(1+8U) − η₀ + 2^−530(1−U) ≥ |z|, because (1−U)^4(1+8U) ≥ 1 and
2^−530(1−U) > 2^−536.4 + η₀. ∎

**L3 (`cmod_bounds`, finite x, y with max(|x|, |y|) ≤ 2^1000).** Claim: lo ≤ |x + iy| ≤ hi. (Not for arbitrary finite
inputs: near the overflow threshold ldexp(r(1−8U), e) can overflow, so lo = +∞. Every caller passes values < 2^13:
box coordinates, half-widths, the Taylor constant term (M), and |C| ± H ≤ 2.)
*Proof.* Let M = max(|x|, |y|) > 0 and M = f·2^e with f ∈ [1/2, 1) (frexp, exact). Then xs = x·2^−e, ys = y·2^−e are
exact, except that the smaller one may underflow (error ≤ η₀). The larger one lies in [1/2, 1), so
xs² + ys² ∈ [1/4, 2]. Two products and one sum give s = (xs² + ys²)(1 + θ₁) with |θ₁| ≤ 2U + 2U² + 12η₀. Then
r = fl(√s) = |z|2^−e(1 + θ), |θ| ≤ U + |θ₁|/2 + U² ≤ 2.5U.

hi: fl(r(1+8U)) ≥ r(1+8U)(1−U) ≥ |z|2^−e(1 + 4U). If ldexp(·, e) is normal it is exact. If it is subnormal, it is
rounded down by < 2^−1074 (≤ η₀ if correctly rounded), and then +2^−1074 restores it, exactly by (F2). If the result is normal, fl(y + 2^−1074) ≥ y by
(F3), and y ≥ |z|(1 + 4U) > |z|.
lo: symmetric with (1 − 8U), −2^−1074, and max(0, ·). ∎

**L4 (`powu`, x ≥ 0, 0 ≤ a ≤ 8).** Claim: x^a ≤ fl(fl(x⋯x)·(1 + 2(a+1)U)) + TINY.
*Proof.* For a = 0 (called with a = 0 in `build_models`), r_c = 1 and the result is fl(fl(1·(1+2U)) + TINY) =
1 + 2U ≥ 1 = x^0. For a ≥ 1: in every call x = ρ ≥ 2^−40 (the source header states x ≥ 2^−60), so the
repeated products do not underflow (x^8 ≥ 2^−320). The chain has a − 1 roundings (the first product 1·x is exact), so
by (F5) r_c ≥ x^a(1−U)^{a−1} ≥ x^a(1−U)^a. The
constant 1 + 2(a+1)U is exact. Then fl(r_c(1 + 2(a+1)U)) ≥ x^a(1−U)^{a+1}(1 + 2(a+1)U) ≥ x^a, because
(1−U)^{a+1} ≥ 1 − (a+1)U and (1 − kU)(1 + 2kU) ≥ 1 for kU ≤ 1/2. The final addition fl(y + TINY) ≥ y by (F3),
since y is a double. ∎

**L5 (`cr_log`).** Let x be positive and finite, and ℓ = cr_log(x). Then |ℓ − log x| ≤ ½ ulp(log x) ≤ U|log x|,
and so |ℓ − log x| ≤ U|ℓ|(1 + 2U). (log x is never subnormal for a double x: the smallest nonzero |log x| is
≈ 2^−53.) Every argument passed to cr_log in the code is positive and finite:
* supT ≥ TINY (up);
* infT > 0 and plo > 0 are tested;
* phi > 0;
* sqrt(c0n2) ≥ 2^−450;
* absm ≥ 2^−41;
* fl((d−1)/d) > 0.

## A.5 Enclosures for |S_i| (`build_models`, interval part)

Fix i, and write δ ∈ B − m. So δ_v ranges over the rectangle [−H_{2v}, H_{2v}] × [−H_{2v+1}, H_{2v+1}], which lies
in the disc of radius |(H_{2v}, H_{2v+1})| ≤ ρ_v (L3). Let c_α be the exact coefficients, with α a multi-index of
degree |α|, and ⟨ĉ_α, r_α⟩ the computed balls (Lemma T).

**Lemma R (remainder and linear majorants).** Define
rT* = r_0 + Σ_{|α|=1} r_α ρ^α + Σ_{|α|≥2} (|ĉ_α| + r_α) ρ^α and linmaj* = Σ_{|α|=1} |ĉ_α| ρ^α.
Then the computed rT ≥ rT* and linmaj ≥ linmaj*, and for all δ,
|T(m+δ) − ĉ_0 − Σ_{|α|=1} ĉ_α δ^α| ≤ rT*.
*Proof.* The inequality follows from |δ^α| ≤ ρ^α and |c_α − ĉ_α| ≤ r_α.

Computation of the terms:
* ra = powu(ρ_i, a) ≥ ρ_i^a (L4).
* rp multiplies in ≤ nb ≤ 5 factors ρ (no underflow, A.2), with a chain of ≤ 5 roundings.
* ac ≥ |ĉ_α| (L2), and ac ≥ 2^−530, so ac·rp and (ac + r)·rp do not underflow.
* r_α·rp may underflow, but it enters only additions.
* For d ≤ 7, rT has at most N ≤ 1 + (6·2^4 − 1) = 96 summands (d = 6: 40). The source comment "≤ 40 terms" is the
  d = 6 count.

Each summand has rounding depth ≤ 7 (rp ≤ 5, times ac or ac + r ≤ 1 + 1). Recursive summation adds ≤ 96 more. By (F5),
rT_c ≥ rT*(1−U)^{103} − 96η₀, and up(up(·)) ≥ rT* by (F6). The same argument applies to linmaj (≤ 6 summands). ∎

**Lemma S (supT, infT).** Let A0lo ≤ |ĉ_0| ≤ A0hi (L3). Then
supT = up(A0hi + linmaj + rT) ≥ A0hi + linmaj + rT ≥ sup_B |T|, and infT ≤ A0lo − linmaj − rT ≤ inf_B |T|.

*Proof of the infT bound.* Write Σ = A0lo + linmaj + rT (reals, all ≥ 0) and X = A0lo − linmaj − rT. The code
computes x2 = fl(fl(A0lo − linmaj) − rT) with x2 ≤ X + 2UΣ + U²Σ (F1). It also computes y = fl(4U·fl(fl(A0lo +
linmaj) + rT)). Then z = fl(x2 − y) and w = fl(z − TINY).

(Σ ≥ rT ≥ TINY > 0 always: rT includes the coefficient radius r_0 ≥ TINY after up.)
Case 1: y is normal. Then y ≥ 4UΣ(1−U)^2 and |z| ≤ Σ(1 + 7U), so
z ≤ x2 − y + U|z| ≤ X − UΣ(1 − 17U) < X, and w ≤ z by (F3).
Case 2: y is subnormal (so Σ < 2^−970). Then x2 ≤ X + 2^−1021 and z ≤ x2 + U|z| ≤ X + 2^−1020. So
w ≤ (z − TINY) + U|z − TINY| < X, since TINY > 2^−997.
In both cases infT = w ≤ X. ∎

**Lemma LB (log bounds, `hi`/`lo`).** For i ≥ 2 (code i > 0), |u_i| ∈ [|m_i| − ρ_i, |m_i| + ρ_i] on B.
* plo = fl(fl(absm_lo − ρ)·(1 − 4U)) ≤ (absm_lo − ρ)(1+U)²(1 − 4U) < absm_lo − ρ ≤ min_B |u_i| when t :=
  absm_lo − ρ > 0. The difference t is never subnormal: absm_lo and ρ are doubles ≥ 2^−41 (A5), hence multiples of
  2^−93, so t ≤ 0 or t ≥ 2^−93.
* Symmetrically, phi ≥ max_B |u_i|.

Exact bounds: log|S_i| ≤ log supT − (n−1) log plo and log|S_i| ≥ log infT − (n−1) log phi.

Computed: ls = cr_log(supT) and dlo = fl((n−1)·cr_log(plo)). The expression hi = fl(fl(ls − dlo) + fl(REL·fl(fl(1 +
|ls|) + dmag))) differs from ls − dlo by the following roundings:
* at most U|ls| from the log (L5),
* (n−1)U|log plo| + U|dlo| from dlo,
* U(|ls| + |dlo|) from the subtraction.

The total is ≤ 10U(1 + |ls| + dmag), where dmag ≥ |dlo|(1−3U). The added allowance (the code evaluates
(1 + |ls|) + dmag) is REL(1 + |ls| + dmag)(1−U)^3 ≥ 890U(1 + |ls| + dmag). The final addition can lower the result by
at most U|hi| ≤ U(|ls| + dmag + 1)(1 + O(U)). So hi ≥ log supT − (n−1) log plo ≥ sup_B log|S_i|. For i = 1 there is
no u-term. The lo bound is symmetric.

Non-finite cases, as coded: hi = +∞ if plo ≤ 0 or supT = ∞; lo = −∞ if infT ≤ 0. **lo can remain finite when plo ≤ 0**
(u_i may vanish on B), **provided infT > 0**: then lo_i = log infT − (n−1) log phi − REL(1 + |li| + dmag) with dmag = (n−1)|log phi|.
This is a valid lower bound for log|S_i| at every point of B with u_i ≠ 0 (there |S_i| = |T_i|/|u_i|^{n−1} ≥
infT/phi^{n−1}; |S_i| → ∞ as u_i → 0). So an E-interval verdict that uses lo_i on such a box holds on the
nondegenerate points, as stated in Theorem A. ∎

## A.6 Affine models of log|S_i| (`build_models`, affine part)

Assume the model is marked valid. That requires:
* grad_ok (c0n2 := fl(ĉ0re² + ĉ0im²) ≥ 2^−900);
* A0lo > 0 and q < 1/2;
* for i ≥ 2: p < 1/2 and m2 ≥ 2^−900;
* all quantities finite.

**Exact model.** On B, T(m+δ) = ĉ_0(1 + w) with w = Σ_{|α|=1} γ_α δ^α + E/ĉ_0, where γ_α = ĉ_α/ĉ_0 and |E| ≤ rT* (Lemma
R). So |w| ≤ (linmaj + rT)/|ĉ_0| ≤ (linmaj + rT)/A0lo ≤ q < 1/2. Taking the principal branch,
log|T| = log|ĉ_0| + Re(γ·δ) + Re(E/ĉ_0) + Re(log(1+w) − w), with
|Re(E/ĉ_0)| ≤ rT/A0lo and |log(1+w) − w| ≤ Σ_{k≥2}|w|^k/k ≤ q²/(2(1−q)) (x ↦ x²/(2(1−x)) is increasing on [0,1)).

Likewise, for i ≥ 2, log|u_i| = log|m| + Re(δ/m) + R′ with |R′| ≤ p²/(2(1−p)), where p ≥ ρ/|m|. With
Re(γδ) = γ_re x − γ_im y and Re(δ/m) = (m_re x + m_im y)/|m|², the code's gradient row (G[2v] = gre, G[2v+1] = −gim)
is exactly the real gradient of Re(γ·δ) − (n−1)Re(δ/m). Hence the exact statement:
log|S_i(m+δ)| = â + G*·(x,y) + ε with |ε| ≤ e*, where â = log|ĉ_0| − (n−1)log|m|, G* is the exact gradient, and
e* = rT/A0lo + q²/(2(1−q)) + (n−1)p²/(2(1−p)).

**Lemma A (computed model is valid).** For all δ ∈ B − m: |log|S_i(m+δ)| − a − G·(x,y)| ≤ e, where (a, G, e) are the
computed quantities. The error sources, all charged to e:

1. **γ.** N = ĉ_re c0re + ĉ_im c0im. By Cauchy–Schwarz, |ĉ_re c0re| + |ĉ_im c0im| ≤ |ĉ||ĉ_0|, so
   |fl(N) − N| ≤ 2.1U|ĉ||ĉ_0| + 2η₀. Also c0n2 = |ĉ_0|²(1 + θ) with |θ| ≤ 2.1U: component squares may underflow, but
   the error 2η₀ ≪ U·2^−900. Hence |gr_c − Re γ| ≤ 4.5U|γ| + U|gr_c| + 2^−172, and the same holds for gi. Since
   |γ| ≤ (|gr_c| + |gi_c|)(1 + 12U) + 2^−171, |γ_c − γ| ≤ 11U(|gr_c| + |gi_c|) + 2^−170. The code charges
   (REL(|gr_c| + |gi_c|) + 2^−160)·ρ_var ≥ |(γ_c − γ)δ| (REL ≥ 900U). For i ≥ 2, adding the u-gradient into gre[v]
   costs one more rounding of ≤ U(|gr_T| + |gr_u|), which is covered by the two REL charges.
2. **u-gradient.** gr = fl(fl(−(n−1)m_re)/m2) has relative error ≤ 4.2U (m2 relative error ≤ 2.1U); the same holds
   for gi. It is charged REL(|gr| + |gi|)ρ_v + 2^−160ρ_v.
3. **Constant term.** sqrt(c0n2) = |ĉ_0|(1 + θ′) with |θ′| ≤ 2.1U, so |la − log|ĉ_0|| ≤ 2.2U + U|la| (L5). Likewise
   |lm − log|m|| ≤ 2.2U + U|lm|. The operations ai −= fl((n−1)·lm) add ≤ U((n−1)|lm| + |ai|). The total is
   ≤ 20U(1 + amag + |ai|), charged REL(1 + amag + |ai|).
4. **Remainder terms.** up(rT/A0lo) ≥ rT/A0lo by (F5)+(F6). Now take x = q²/(2(1−q)), computed with roundings
   q·q, 1 − q (up ≤ U, which lowers the value), 2× (exact), ÷, and ×(1+4U). q is only known to be ≳ 2^−584, so q·q,
   the ÷ and the × may underflow. Each underflow error η₀ is then multiplied by at most (1+4U)/(2(1−q)) ≤ 1 + 4U. So
   x_c ≥ x(1−U)^4 − 3η₀(1+4U), and up(·) ≥ x by (F6) with k = 4, N = 3. The term (n−1)p²/(2(1−p)) has one more
   rounding ((n−1)·p): x_c ≥ x(1−U)^5 − 4η₀(1+4U), and again up(·) ≥ x. (The factor (1+4U) in the code is not needed
   for this argument.) q ≥ (linmaj + rT)/A0lo
   and p ≥ ρ/|m| hold because each is up(·) of a quotient whose numerator is an upper bound and whose denominator is
   a lower bound (absm_lo(1−4U) rounded is < absm_lo ≤ |m|).
5. **Accumulation.** ei is a sum of ≤ 2 + 2nv nonnegative terms (≤ 12 roundings). The final
   up(ei·(1 + 10^−12) + ...) makes the result ≥ the exact sum, since 10^−12 ≫ 12U.

Therefore e ≥ e* + (all computed-vs-exact discrepancies in a and G, integrated over B), which is the claim. ∎

## A.7 Tests (`eval_box`)

Let L_c = `logc_lo` = fl(cr_log(fl((d−1)/d)) − 10^−14). fl((d−1)/d) has relative error ≤ U, and |log c| ≤ log(4/3)
< 0.29 for d ≥ 4, so |cr_log(fl(c)) − log c| ≤ 1.3U. So L_c ≤ log c.

**F-interval.** fl(hi + fl(MARG(1 + |hi|))) < L_c ⇒ hi < L_c by (F3), since the margin is ≥ 0. Then
log|S_i| < log c on B (Lemma LB).

**F-affine.** Let s* = a + Σ_k|G_k|H_k + e be the exact real value with the computed a, G, e. Then sup_B log|S_i| ≤ s*
(Lemma A; |G·(x,y)| ≤ Σ|G_k|H_k on the rectangle). Also s_c is a (D+2)-term recursive sum of products, so
|s_c − s*| ≤ (D+3)U·mag*(1 + U) + Dη₀, where mag* is the same sum with |a|, D + 3 ≤ 13, and Dη₀ covers products
|G_k|·H_k that underflow. By (F3), fl(s_c + MARG·(1 + mag_c)) < L_c gives s_c + MARG(1 + mag_c)(1 − 2U) < L_c. So
s* ≤ s_c + 14U·mag_c + Dη₀ < L_c, since MARG·1 absorbs Dη₀.

**E-interval.** fl(mlo − m₁) > fl(mhi + m₂) with m₁, m₂ ≥ 0 implies mlo > mhi (F3). So for some i, j,
inf log|S_i| ≥ lo_i > hi_j ≥ sup log|S_j| on B. Moreover i ≠ j. B has interior points with all u_v ≠ 0 (H_k > 0), and
at such a point every finite lo_i satisfies lo_i ≤ log|S_i| ≤ hi_i (Lemma LB; if lo_i is finite then infT > 0, so
T_i ≠ 0). Hence lo_i ≤ hi_i for every i, and mlo > mhi forces the index of the maximum lo to differ from the index of
the minimum hi.

**E-affine.** Let lhs* = |a_i − a_j| and rhs* = e_i + e_j + Σ_k|G_ik − G_jk|H_k. The rounding errors satisfy
|lhs_c − lhs*| ≤ U·lhs_c and |rhs_c − rhs*| ≤ (D+4)U·rhs_c(1+U) + Dη₀ (Dη₀ for underflowing products
|G_ik − G_jk|·H_k). The test is
fl(lhs_c − rhs_c) > MARG(1 + |a_i| + |a_j| + rhs_c) (computed). It implies lhs* − rhs* > 0, because the margin exceeds
2U·lhs_c + 15U·rhs_c + U|lhs_c − rhs_c| + Dη₀. Then on B, |l_i − l_j| ≥ lhs* − rhs* > 0, i.e. |S_i| ≠ |S_j|.

**L (uniform weights).** On B ∩ E all log|S_i| are equal, so (1/n)Σ_i log|S_i| = log F. Let
s* = Σ_i (a_i + e_i) + Σ_k |Σ_i G_ik| H_k. Then Σ_i log|S_i| ≤ s* on B. The inner sums g_k = Σ_i G_ik carry error
≤ nU·Σ_i|G_ik|, and the outer sum has ≤ 2n + D terms. The products |g_k|·H_k and ga·H_k can underflow (≤ 2Dη₀ in
total). So |s_c − s*| ≤ (3n + D + 2)U·mag_c(1 + U) + 2Dη₀ ≤ 31U·mag_c + 2Dη₀ for d ≤ 7 (3n + D + 2 = 30 at d = 7). The right-hand side is
Y_c = fl(fl(n·L_c) − MARG) ≤ nL_c + 2.1U·n|L_c| − MARG(1 − U). The test fl(s_c + MARG(1 + mag_c)) < Y_c gives, by (F3),
s_c + MARG(1 + mag_c)(1 − 2U) < Y_c, hence
s* < nL_c + 2.1U·n|L_c| − MARG(1 − U) − MARG(1 + mag_c)(1 − 2U) + 31U·mag_c + 2Dη₀ < nL_c ≤ n log c. So log F < log c on B ∩ E. More precisely, (1/n)Σ log|S_i| < log c on
all of B, which is the uniform strict form the Prop R corollary needs.

**Uniformity.** Each implication above gives a strict inequality whose gap is bounded below on B. Indeed B is
compact, and each bound is a fixed real computed from (a, G, e, hi, lo) with slack ≥ MARG. This is the "uniform
strict test" of the Prop R corollary.

## A.8 Geometry and exclusion

**Outside / symmetry (`modrange`, `geom_discard`).** In the domain, dx = max(0, |C_x| − H_x) and the analogous dy,
ex, ey are exact (A.2), and min_B |u_v| = |(dx, dy)|, max_B |u_v| = |(ex, ey)|. By L3, mn_c ≤ min |u_v| and
mx_c ≥ max |u_v|. Multiplying by (1 ∓ 4U) and rounding moves each further in the safe direction (F3).
* `mn > 1 + 10^−12` ⇒ |u_v| > 1 on B.
* `fl(mx(1 + 10^−12)) < fl(mn′(1 − 10^−12))` ⇒ mx < mn′ (F3) ⇒ |u_v| < |u_{v+1}| on B.
(Here and below, 1 ± 10^−12, 10^−13 denote the doubles the literals round to; only 1 + 10^−12 > 1 > 1 − 10^−12 and
10^−13 > 0.99·10^−13 are used.)
* `C[1] + H[1] < 0` is exact ⇒ Im u_2 < 0 on B.

**Exclusion, mode 2 (`excluded`).**
Let p̂ be the computed equality point. `c/regress/eqpoints_v2.py` **certifies** ‖p̂ − p‖² < (1.1·10^−16)² in arb
(python-flint ball arithmetic, 300 bits; exact points via arb cos_pi/sin_pi; the doubles enter exactly) for every
output point, d = 4..7 (runs/regress_v2/eqpoints_v2.txt: ALL PASS). It also compares with mpmath in both directions,
which is provenance only, since soundness needs only the computed → exact direction. Each computed point is within 1.1·10^−16 (Euclidean) of an exact
equality point, each exact point has a computed point that close, and the counts are (n−1)!. The maximum Euclidean
errors are d = 4: 7.1·10^−17, d = 5: 0, d = 6: 8.0·10^−17, d = 7: 1.0·10^−16. We use ‖p̂ − p‖ ≤ 1.6·10^−16. (Even the
weaker per-coordinate bound 10^−15 asserted by c/regress/regress_v2.py, i.e. ≤ 3.2·10^−15 in norm, would leave
9.6·10^−14 of slack below.)

For u ∈ B, |u − p̂|₂ ≤ (Σ_k (|C_k − p̂_k| + H_k)²)^{1/2} =: Φ. Each fx_c = fl(fl(|fl(C − p̂)|) + H) ≥ (|C − p̂| + H)
(1−U)² − η₀. far2 is a sum of 2nv ≤ 10 squares, so far2_c ≥ Φ²(1−U)^{15} − 20η₀, and
fl(√far2_c) ≥ Φ(1 − 8.5U) − 2^−535. The test is fl(Y) ≤ r_excl with Y = fl(√far2_c·(1 + 10^−12)) + 10^−13, the last
sum taken exactly. If Y > r_excl + ½ulp(r_excl), then fl(Y) > r_excl. So the test implies Y ≤ r_excl(1 + U), i.e.
Φ(1 − 8.5U)(1 + 10^−12)(1 − U) − 2^−534 + 10^−13 ≤ r_excl(1 + U). Since (1 − 8.5U)(1 + 10^−12)(1 − U) > 1, we get
Φ ≤ r_excl(1 + U) − 10^−13 + 2^−534. Therefore sup_B |u − p| ≤ Φ + 1.6·10^−16 < r_excl − 9·10^−14 < r_excl, and B lies
in the closed ball. Here r_excl is the double parsed from "0.05", which is 0.05 + 2.8·10^−18. Since
sup_B |u − p| < r_excl − 9·10^−14 < 0.05, B also lies in the closed ball of exact radius 0.05 used by
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

