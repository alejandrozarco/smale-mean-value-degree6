# Local certificate at the equality configuration (d = 6; also d = 5)

The code is `local/local_cert.py`, which uses python-flint arb/acb balls at 128 bits, together with
`local/exact_facts.py`, which does exact arithmetic in Q(ω). The run logs are
`runs/local_cert_v3_d6_T0.0527.log` and `runs/local_cert_v3_d5_T0.0527.log`; each takes about 1 minute. The logs
`runs/local_cert_v3_d6_T{0.01,0.02,0.03,0.04}.log` are the same computation at smaller radii (used in the figure).

## Statement

Let d = 6, n = 5, c = 5/6. Chart: b_1 = 1 and b_j = ω^{j−1} e^{ε_j} (j = 2..5, ω = e^{2πi/5}), with ε ∈ C^4,
x = (Re ε_2, Im ε_2, …) ∈ R^8 and t = |x|_2. Put ℓ_i = log|S_i|, L = (1/n) Σ ℓ_i, and E = {|S_1| = … = |S_5|}.

**Claim.** For every x ∈ E with 0 < |x|_2 ≤ T := 0.0527, we have log F(x) = L(x) ≤ log(5/6) − 0.00757 |x|_2²
< log(5/6). At x = 0, F = 5/6. T is the exact decimal, and the code uses an arb ball that contains it.

**Hand-over to the branch and bound.** `c/smale_bb_v2` (excl_mode 2) sets a box aside only when it has verified that
the box lies in a closed Euclidean u-ball Σ_j |u_j − p_j|² ≤ r² (r = 0.05) around one of the 24 equality points p.
Relabelling maps each such ball to the ball around p = (ω̄, ω̄², ω̄³, ω̄⁴), where u_j = ω^{−(j−1)} e^{−ε_j}. With
w_j = u_j/p_j − 1 we have |ε_j| = |log(1 + w_j)| ≤ |w_j|/(1 − |w_j|), so |ε|_2 ≤ r/(1 − r) = 0.0526316 < T.
`c/check_done.py` checks r/(1 − r) ≤ 0.0527 exactly for the IEEE value of r.

## Proof structure (all constants computed in arb)

1. **Taylor models.** h_i = log S_i − log c is holomorphic on the polydisk |ε|_∞ ≤ T. Taylor-model arithmetic (degree
   N = 6 in 4 complex variables) gives h_i = P_i + R_i. Here P_i is the degree-≤ 6 Taylor polynomial with acb ball
   coefficients, and |R_i| ≤ r_i on the polydisk (r_i ≤ 2.9e−6). The code uses |R_i| ≤ r_i (t/T)² for t ≤ T.
2. **Exact low-order facts** (`local/exact_facts.py`, asserted by `local_cert.py`). In Q(ω): S_i(0) = (d−1)/d for all
   i, and Σ_i ∂S_i/∂ε_j(0) = 0 for j = 2..n. Hence h_i(0) = 0, L(0) = log c and ∇L(0) = 0 exactly (checked for
   d = 4..7). With J = Dℓ(0) (a 5×8 matrix whose rows sum to 0), G_i := ℓ_i − L = (Jx)_i + G_i^{≥2}(x).
3. **Tube around K = ker J.** On E all G_i = 0, so y := |Jx| = |G^{≥2}(x)|. Write x = Πx + w with Π = VVᵀ (formed in
   arb) and |w| ≤ β y + γ t, where β = |J⁺| ≤ 1.7929 and γ ≤ 2.04e−15. With:
   - q = sup_{|v|≤1} |G²(Πv)| ≤ 0.23673 (a rigorous branch-and-bound over the unit sphere in K);
   - b₂ = (Σ_i |A_i|²)^{1/2} ≤ 1.2650;
   - the degree 3..6 parts and remainders of G,

   the fixed-point inequality y ≤ q t² + 2b₂|Π| t(βy + γt) + b₂(βy + γt)² + g₃(t), started from y ≤ |J| t, gives
   |Jx| ≤ κ t² with κ = 0.4787 for all x ∈ E with |x| ≤ T.
4. **Quadratic part.** For σ = 5, H_σ := Hess L(0) − σ JᵀJ ≤ −μ I with μ = 0.03192459 (arb LDLᵀ of −H_σ − μI).
   Hence ½xᵀ Hess L(0) x ≤ −(μ/2)t² + (σ/2)κ² t⁴.
5. **Cubic part.**
   - sup_{|z|=1} |L³(Vz)| ≤ 0.0325473, by a rigorous sphere branch-and-bound with exact rational sphere-intersection
     tests.
   - |L³(a+w) − L³(a)| ≤ C₃((|a|+|w|)³ − |a|³) with C₃ ≤ 0.5192.
   - Degrees 4..6 are bounded by Σ_α |C_α| m_α t^k, with m_α = Π_j (α_j/k)^{α_j/2} computed in arb.
   - The remainder contributes (1/n) Σ r_i (t/T)².
6. **Conclusion.** L(x) − log c ≤ t² g(t), where every term of g is increasing in t:
   - quadratic: −0.015962;
   - penalty: 0.001591;
   - cubic on K: 0.001715;
   - cubic off K: 0.003883;
   - degrees 4–6: 0.000782;
   - remainder: 0.000417.

   So g(T) ≤ −0.0075740986 (an arb upper bound; the displayed terms are rounded).

For d = 5 the same code gives g(T) ≤ −0.02284 (μ = 0.05039). At T = 0.08 the d = 6 bound is not negative.

## Numerical comparison (not rigorous)

`local/e_profile.py` projects points onto E by Newton's method, in float64, from 400 random K-directions per radius.
The sampled minimum of (log c − log F)/|x|² on E is 0.0198–0.0215 for t = 0.005 … 0.2 (`runs/e_profile_d6.tsv`).
`local/check_tm.py` samples |h_i − P_i| at random points of the polydisk and compares it with r_i.
