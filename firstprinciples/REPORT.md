# Independent first-principles check, $`d = 5`$ and $`d = 6`$

Status: computational checks, not peer reviewed; produced by AI models (see [`AI_DISCLOSURE.md`](../AI_DISCLOSURE.md)).

This check was written by a separate agent instance. That instance did not read the evaluator, the arb checker or its
specification, the error analysis, `REDUCTION.md`, `local/` or `local2/`. Its only inputs were:
- the problem statement;
- the coordinate conventions;
- the claim table;
- the formats of the tree and tasks files.

Its critical values $`V_i = P(b_i)/b_i`$ come from multiplying out $`P' = \prod_j (1 - z/b_j)`$ and integrating the
coefficients (`fpcore.py`), not from any integral formula.

Run:
- `python3 derive.py` (algebra, about 20 s; log `derive.log`);
- `cc -O2 -o walk walk.c`, then `python3 sample_tree.py D all` with the tree in `FP_TREE_D5` / `FP_TREE_D6`.
  The $`d = 6`$ tree is the release asset `d6.canonical.tree.zst`. The $`d = 5`$ tree is written by `c/export_tree.c`.
  The second $`d = 6`$ pass set `FP_DEEP=66` (every leaf at depth $`\ge 66`$); the selection it printed is in
  `d6_pass2_deep.log`. Intermediate files go to `FP_WORK` (about 0.6 GB for both degrees).

## Results per claim
**(1) Coordinates** x = (Re u₂, Im u₂, …) are consistent with the trees: every outside and symmetry leaf satisfies its
claim under this convention (d = 5: 37 156 and 205 334; d = 6: 18.66 M and 173.6 M).

**(2) Formulas** S_i = V_i: PASS.
- Hand derivation: V_i = ∫₀¹ ∏_j (1 − t b_i/b_j) dt. For i ≥ 2 the factors j = 1 and j ∉ {1, i} give the
  denominator u_i^{n−1}, so V_i = T_i/u_i^{n−1}.
- sympy, d = 4–7: V_i built from P minus the claimed S_i is exactly 0 for every i.
- 30-digit quadrature agrees to 3e-29.

**(3) Domain reduction:** PASS. Each symmetry was checked symbolically for d = 4–7.
- Scaling: P(λz)/λ keeps V_i, index by index.
- Any permutation σ of all n critical points: V_i(b∘σ) = V_σ(i)(b) (generators (1 2) and the n-cycle).
- Conjugation: V_i ∈ Q(b), so V_i(conj b) = conj V_i(b).
- Choosing b₁ of minimal modulus needs permutations of all n indices, which holds.
- The reduction sends 2000 random configurations per degree into R with the sorted |V| preserved.
- Closedness and u_j = 0 only add points to be covered. The top boxes are the side-1/2 grid cells of [−1,1]^D with
  Im u₂ ≥ −1/2, and every omitted cell has Im u₂ ≤ −1/2 < 0.

**(4) Equality points:** PASS, exactly in Q[w]/Φ_n for all 2, 6, 24 and 120 points (d = 4–7). b_j = ω^{a_j} gives
P' = 1 − zⁿ, P = z − z^d/d and V_i = 1 − 1/d.

**(5) Hand-over to the local chart:** PASS.
- u_j = p_j e^{−ε_j} and w_j = u_j/p_j − 1 give |w_j| = |u_j − p_j|, and
  |ε_j| ≤ −log(1 − |w_j|) = |w_j| g(|w_j|) with g(s) = −log(1−s)/s increasing. Hence ‖ε‖ ≤ −log(1 − r).
- The bound is sharp, attained at u₂ = (1 − r)p₂. At r = 1/20 it equals 0.0512933.
- Numerically the largest ratio was 0.99992.
- Relabelling b₂..b_n is an isometry: it sends each equality point to the identity point and permutes V₂..V_n.

**Discrepancy (no consequence).** The task headers' r_excl is 0x3fa999999999999a = 1/20 + 2.8e-18. Every excluded
leaf was re-checked against the exact ball of radius 1/20 and lies inside it: smallest margin 1/400 − max‖x−p‖² is
5.456e-10 for d = 6 (60-digit re-check) and 1.37e-6 for d = 5.

## Trees (all records; `walk.c`)
| | d = 5 | d = 6 |
|---|---|---|
| records | 3072 | 49152 |
| each parses as exactly one complete tree; ids unique and complete; size = Σ(14 + nbytes); #internal = #leaves − 1 | yes | yes |
| invalid codes / code 6 | 0 / 0 | 0 / 0 |
| max depth; max splits of one coordinate | 47; 8 | 71; 9 |
| F / E / L | 757 590 / 143 562 / 100 008 | 561.2 M / 145.4 M / 37.9 M |
| outside / symmetry / excluded | 37 156 / 205 334 / 2 916 | 18.66 M / 173.6 M / 591 073 |

**Geometric leaves:** all were checked exactly, with 0 violations.

| smallest margin | d = 5 | d = 6 |
|---|---|---|
| outside | 6.1e-5 | 3.8e-6 |
| symmetry | 4.6e-5 | 9.5e-7 |
| excluded | 1.37e-6 | 5.456e-10 |

**F/E/L leaves** (`sample_tree.py`): the true V_i were evaluated at sample points, with 0 violations.

What was sampled:
- d = 5: all 1 001 160 leaves, each at 64 corners, the centre and 12 random points (77 M points).
- d = 6, pass 1: 1.33 M hash-sampled leaves over all 48 138 task ids that have F/E/L leaves.
- d = 6, pass 2: 0.85 M leaves, including every leaf at depth ≥ 66 (410 633 leaves), each with 256 corners.
- d = 6 total: 423 M points.
- In float64, anything within 1e-9 plus an error estimate of a threshold was recomputed in 50-digit mpmath. That
  happened once, and it passed.

How each claim was tested:
- **F:** one fixed i with |V_i| < c at all points of the leaf.
- **E:** one fixed pair with log|V_i| − log|V_j| of the same strict sign at all points.
- **L:** Σ log|V_i| < n log c at every point.

| | d = 5 | d = 6 |
|---|---|---|
| smallest F margin | 3.13e-6 | 1.287e-6 (task 18097, depth 70) |
| smallest E margin (log ratio) | 2.50e-5 | 1.122e-5 |
| smallest L margin (log) | 9.07e-5 | 6.49e-5 |
| largest min_i \|V_i\| sampled | c − 1.27e-4 | c − 7.6e-5 |

Further searches:
- An adversarial L-BFGS-B search inside the closest 950 (d = 5) and 2 500 (d = 6) leaves found nothing below the
  corner margins.
- A Nelder–Mead search for min_i |V_i| > c peaked at 0.7999927 (d = 5) and 0.8333076 (d = 6), both below c. This is
  weak evidence: the method converges poorly on nonsmooth objectives.

## Points that need an extra argument
1. E leaves do not by themselves give min_i |V_i| ≤ c. Two things are needed:
   - (a) A maximiser, or the supremum, must be handled although R is not compact (u_j → 0).
   - (b) The maximal value must be attained with all |V_i| equal.

   1 642 (d = 5) and 1 070 847 (d = 6) E leaves contain points with some u_j = 0.
2. Excluded leaves rely on the local certificate covering ε-radius ≥ 0.0512933. That certificate was not checked here.

Note, not part of this check: item 1 is what `REDUCTION.md` supplies. Proposition R and its corollary work with the
closure of the image, so no attainment is assumed, and use only nondegenerate points. Item 2 is supplied by `local/`
($`T = 0.0527`$) and `local2/` ($`T = 1/19`$).

## Limits
- This is falsification by sampling, not a proof: about 0.3 % of the d = 6 F/E/L leaves, at finitely many points
  each, without interval arithmetic.
- The evaluator, its error analysis, the local certificates and the global argument were not examined.
- d = 7 was covered only by the algebraic checks.
