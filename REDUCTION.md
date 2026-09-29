# Equal-modulus reduction

## Notation

- d ≥ 3 and n = d − 1.
- Ω = (C*)^n is the set of critical-point tuples B = (b_1, …, b_n). All b_i ≠ 0, which is equivalent to P′(0) ≠ 0,
  and repetitions are allowed.
- S = (S_1, …, S_n): Ω → C^n, with S_i(B) = P(b_i)/(b_i P′(0)) = ∫₀¹ Π_j (1 − t b_i/b_j) dt. Each S_i is a Laurent
  polynomial in B and is invariant under B ↦ λB.
- M_d := sup_Ω min_i |S_i|. Smale's theorem gives M_d ≤ 4, and z − z^d/d gives M_d ≥ (d−1)/d.

## Proposition R

There is a point z* in the Euclidean closure of S(Ω) with |z*_1| = … = |z*_n| = M_d.

*Proof.*

**Step 1: the image is dense in a hypersurface.** Let Ω′ = {B ∈ Ω : b_1 = 1} ≅ (C*)^{n−1}. It is irreducible, and
S(Ω′) = S(Ω). At B = (1, …, 1) we have

- ∂S_i/∂b_j = −(n−1)/(n(n+1)) for i = j;
- ∂S_i/∂b_j = 1/(n(n+1)) for i ≠ j.

With b_1 fixed, the minor on rows and columns 2..n has determinant (−1)^{n−1}/(n(n+1)^{n−1}) ≠ 0. For d = 6 this is
1/6480. `local/rank_check.py` also checks rank n − 1 exactly at a Gaussian-rational point.

So the Zariski closure of S(Ω) is an irreducible hypersurface Z_f = {f = 0}. By Chevalley's theorem S(Ω) is
constructible, so it contains a dense Zariski-open subset of Z_f. That subset is dense in Z_f in the Euclidean
topology (Mumford, *The Red Book*, I.10, Thm 1). Hence Z_f is the Euclidean closure of S(Ω).

**Step 2.** Let Z* = Z_f ∩ (C*)^n and C := sup_{z ∈ Z*} min_i |z_i|. The function min_i |z_i| is continuous, and
S(Ω) ∩ (C*)^n is dense in Z*. So C = M_d, and C is finite by Smale's bound.

**Step 3: the amoeba argument.** Let A = Log(Z*) ⊂ R^n, where Log(z) = (log|z_1|, …, log|z_n|). A is closed. Every
connected component of R^n \ A is convex (Gelfand–Kapranov–Zelevinsky 1994, Ch. 6 Cor. 1.6; Forsberg–Passare–Tsikh
2000).

Let γ = log C and p = (γ, …, γ). The open orthant Q = {x : x_i > γ for all i} misses A. Suppose p ∉ A. Then some ball
B(p, r) misses A. The set Q ∪ B(p, r) is connected, so it lies in one component U of R^n \ A, and by convexity
conv(Q ∪ B(p, r)) ⊂ U.

Take 0 < δ < r/√n. For any x with min_i x_i > γ − δ/2, write x = ½(p − δ·1) + ½q with q = 2x − p + δ·1. Then q ∈ Q,
so x ∈ U. Hence no point of A has min_i x_i > γ − δ/2, which contradicts γ = sup_A min_i x_i.

Therefore p ∈ A, and z* can be taken in Z* with |z*_i| = M_d for all i. ∎

## Corollary (the form used by the computation)

**Chart.** Let D be the closed chart: b_1 = 1 is a critical point of minimal modulus, u_j = 1/b_j and |u_j| ≤ 1
(j = 2..n), cut down by the symmetries |u_2| ≥ … ≥ |u_n| and Im u_2 ≥ 0. D includes the degenerate boundary u_j = 0.
A point is *nondegenerate* if all u_j ≠ 0.

**Hypothesis.** Suppose D is covered by finitely many closed boxes β, and each β satisfies at least one of the
following. Here c = (d−1)/d, each U < c or U < log c, each g > 0, and every inequality holds at all nondegenerate
points of β.

- β is contained in a closed Euclidean u-ball of radius ρ around an equality point, or a symmetric image of one.
  On that ball every nondegenerate point with |S_1| = … = |S_n| has min_i |S_i| ≤ c.
- (F) sup |S_i| ≤ U < c for some i.
- (E) |log|S_i| − log|S_j|| ≥ g > 0 for some i, j.
- (L) Σ_i w_i log|S_i| ≤ U < log c for some real weights w with Σ w_i = 1.

**Conclusion.** Then M_d = c.

*Proof.* Take z* = lim S(B_k) from Proposition R.

1. Normalise each B_k into D. Normalisation permutes the S_i and possibly conjugates them, so the moduli are only
   permuted. Pass to a subsequence with a fixed relabelling and a fixed conjugation choice. Then u_k → u* ∈ D.
2. Every coordinate of the limit has modulus M_d > 0, so the logarithms converge. Infinitely many u_k lie in one box β.
3. That β cannot satisfy (F), (E) or (L), because each would give a contradiction in the limit. So β lies in one of
   the balls.
4. The balls lie at positive distance from {u_j = 0}, and S is continuous there. So z* is an equal-modulus value
   S(B*), and M_d ≤ c. ∎

Nothing here uses continuity of S at u_j = 0, and none is available there.

## Literature

Crane (PhD thesis, Cambridge 2003, Thm 5.4; *Comput. Methods Funct. Theory* 6 (2006) 145–163) gives equal-modulus
statements for maximisers, under the hypothesis that a maximiser exists. T. W. Ng (survey talk, 2018) notes a proof
of Crane's max–min theorem via convexity of amoeba complements. The statement above removes the attainment
hypothesis: it only uses a point of the closure of the image.
