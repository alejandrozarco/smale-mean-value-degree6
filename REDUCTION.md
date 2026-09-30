# Equal-modulus reduction

## Notation

- $`d \ge 3`$ and $`n = d - 1`$.
- $`\Omega = (\mathbb{C}^*)^n`$ is the set of critical-point tuples $`B = (b_1, \dots, b_n)`$. All $`b_i \ne 0`$, which is
  equivalent to $`P'(0) \ne 0`$, and repetitions are allowed.
- $`S = (S_1, \dots, S_n) \colon \Omega \to \mathbb{C}^n`$, with
  $`S_i(B) = \frac{P(b_i)}{b_i\,P'(0)} = \int_0^1 \prod_j \bigl(1 - t\,b_i/b_j\bigr)\,dt`$.
  Each $`S_i`$ is a Laurent polynomial in $`B`$ and is invariant under $`B \mapsto \lambda B`$.
- $`M_d := \sup_\Omega \min_i \lvert S_i\rvert`$. Smale's theorem gives $`M_d \le 4`$, and $`z - z^d/d`$ gives
  $`M_d \ge (d-1)/d`$.

## Proposition R

There is a point $`z^*`$ in the Euclidean closure of $`S(\Omega)`$ with $`\lvert z^*_1\rvert = \dots = \lvert z^*_n\rvert = M_d`$.

*Proof.*

**Step 1: the image is dense in a hypersurface.** Let $`\Omega' = \{B \in \Omega : b_1 = 1\} \cong (\mathbb{C}^*)^{n-1}`$.
It is irreducible, and $`S(\Omega') = S(\Omega)`$. At $`B = (1, \dots, 1)`$ we have

- $`\partial S_i/\partial b_j = -(n-1)/(n(n+1))`$ for $`i = j`$;
- $`\partial S_i/\partial b_j = 1/(n(n+1))`$ for $`i \ne j`$.

With $`b_1`$ fixed, the minor on rows and columns $`2, \dots, n`$ has determinant
$`(-1)^{n-1}/\bigl(n(n+1)^{n-1}\bigr) \ne 0`$. For $`d = 6`$ this is $`1/6480`$. `local/rank_check.py` also checks rank
$`n - 1`$ exactly at a Gaussian-rational point.

So the Zariski closure of $`S(\Omega)`$ is an irreducible hypersurface $`Z_f = \{f = 0\}`$. By Chevalley's theorem
$`S(\Omega)`$ is constructible, so it contains a dense Zariski-open subset of $`Z_f`$. That subset is dense in $`Z_f`$ in the
Euclidean topology (Mumford, *The Red Book*, I.10, Thm 1). Hence $`Z_f`$ is the Euclidean closure of $`S(\Omega)`$.

**Step 2.** Let $`Z^* = Z_f \cap (\mathbb{C}^*)^n`$ and $`C := \sup_{z \in Z^*} \min_i \lvert z_i\rvert`$. The function
$`\min_i \lvert z_i\rvert`$ is continuous, and $`S(\Omega) \cap (\mathbb{C}^*)^n`$ is dense in $`Z^*`$. So $`C = M_d`$, and $`C`$ is
finite by Smale's bound.

**Step 3: the amoeba argument.** Let $`A = \mathrm{Log}(Z^*) \subset \mathbb{R}^n`$, where
$`\mathrm{Log}(z) = (\log\lvert z_1\rvert, \dots, \log\lvert z_n\rvert)`$. $`A`$ is closed. Every connected component
of $`\mathbb{R}^n \setminus A`$ is convex (Gelfand–Kapranov–Zelevinsky 1994, Ch. 6 Cor. 1.6; Forsberg–Passare–Tsikh 2000).

Let $`\gamma = \log C`$ and $`p = (\gamma, \dots, \gamma)`$. The open orthant $`Q = \{x : x_i \gt  \gamma \text{ for all } i\}`$
misses $`A`$. Suppose $`p \notin A`$. Then some ball $`B(p, r)`$ misses $`A`$. The set $`Q \cup B(p, r)`$ is connected, so it lies
in one component $`U`$ of $`\mathbb{R}^n \setminus A`$, and by convexity $`\mathrm{conv}(Q \cup B(p, r)) \subset U`$.

Take $`0 \lt  \delta \lt  r/\sqrt{n}`$. For any $`x`$ with $`\min_i x_i \gt  \gamma - \delta/2`$, write
$`x = \tfrac12(p - \delta\mathbf{1}) + \tfrac12 q`$ with $`q = 2x - p + \delta\mathbf{1}`$. Then $`q \in Q`$, so $`x \in U`$.
Hence no point of $`A`$ has $`\min_i x_i \gt  \gamma - \delta/2`$, which contradicts $`\gamma = \sup_A \min_i x_i`$.

Therefore $`p \in A`$, and $`z^*`$ can be taken in $`Z^*`$ with $`\lvert z^*_i\rvert = M_d`$ for all $`i`$. $`\square`$

## Corollary (the form used by the computation)

**Chart.** Let $`D`$ be the closed chart: $`b_1 = 1`$ is a critical point of minimal modulus, $`u_j = 1/b_j`$ and
$`\lvert u_j\rvert \le 1`$ ($`j = 2, \dots, n`$), cut down by the symmetries $`\lvert u_2\rvert \ge \dots \ge \lvert u_n\rvert`$
and $`\mathrm{Im} u_2 \ge 0`$. $`D`$ includes the degenerate boundary $`u_j = 0`$. A point is *nondegenerate* if all
$`u_j \ne 0`$.

**Hypothesis.** Suppose $`D`$ is covered by finitely many closed boxes $`\beta`$, and each $`\beta`$ satisfies at least one
of the following. Here $`c = (d-1)/d`$, each $`U \lt  c`$ or $`U \lt  \log c`$, each $`g \gt  0`$, and every inequality holds at all
nondegenerate points with all $`S_i \ne 0`$ of $`\beta`$.

- $`\beta`$ is contained in a closed Euclidean $`u`$-ball of radius $`\rho`$ around an equality point, or a symmetric image
  of one. On that ball every nondegenerate point with $`\lvert S_1\rvert = \dots = \lvert S_n\rvert`$ has
  $`\min_i \lvert S_i\rvert \le c`$.
- (F) $`\sup \lvert S_i\rvert \le U \lt  c`$ for some $`i`$.
- (E) $`\bigl\lvert \log\lvert S_i\rvert - \log\lvert S_j\rvert \bigr\rvert \ge g \gt  0`$ for some $`i \ne j`$.
- (L) $`\sum_i w_i \log\lvert S_i\rvert \le U \lt  \log c`$ for some real weights $`w`$ with $`\sum_i w_i = 1`$.

**Conclusion.** Then $`M_d = c`$.

*Proof.* Take $`z^* = \lim S(B_k)`$ from Proposition R.

1. Normalise each $`B_k`$ into $`D`$. Normalisation permutes the $`S_i`$ and possibly conjugates them, so the moduli are only
   permuted. Pass to a subsequence with a fixed relabelling and a fixed conjugation choice. Then $`u_k \to u^* \in D`$.
2. Every coordinate of the limit has modulus $`M_d \gt  0`$, so the logarithms converge. Infinitely many $`u_k`$ lie in one
   box $`\beta`$ of the cover. (In the computation, leaves discarded as outside the chart or by symmetry contain no point
   of $`D`$, so the remaining leaves — F, E, L and excluded — cover $`D`$ and form the cover in the hypothesis.)
3. That $`\beta`$ cannot satisfy (F), (E) or (L), because each would give a contradiction in the limit. So $`\beta`$ is
   contained in one of the closed balls, and since the ball is closed, $`u^*`$ lies in it too.
4. The balls lie at positive distance from $`\{u_j = 0\}`$, and $`S`$ is continuous there. So $`z^*`$ is an equal-modulus
   value $`S(B^*)`$, and $`M_d \le c`$. $`\square`$

Nothing here uses continuity of $`S`$ at $`u_j = 0`$, and none is available there.

## Literature

Crane (PhD thesis, Cambridge 2003, Thm 5.4; *Comput. Methods Funct. Theory* 6 (2006) 145–163) gives equal-modulus
statements for maximisers, under the hypothesis that a maximiser exists. T. W. Ng (survey talk, 2018) notes a proof
of Crane's max–min theorem via convexity of amoeba complements. The statement above removes the attainment
hypothesis: it only uses a point of the closure of the image.
