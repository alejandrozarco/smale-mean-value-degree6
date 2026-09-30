# Independent arb check of a subdivision tree — specification

The checker (`arbcheck.c`) was written from this specification and the mathematics only, independently of
`c/smale_bb_v2.c` and its error analysis, in a different arithmetic (FLINT/arb ball arithmetic).

## 1. The functions

Fix $`d \in \{5, 6\}`$ (also 4 and 7 should work), $`n = d - 1`$, and $`nv = d - 2`$ complex variables $`u_2, \dots, u_n`$. The real
coordinates are $`x \in \mathbb{R}^D`$ with $`D = 2 \cdot nv`$, ordered $`x = (\mathrm{Re}\, u_2, \mathrm{Im}\, u_2, \mathrm{Re}\, u_3, \mathrm{Im}\, u_3, \dots)`$. So variable $`u_{v+2}`$
($`v = 0..nv-1`$) has real part $`x_{2v}`$ and imaginary part $`x_{2v+1}`$.

- $`S_1(u) = \int_0^1 (1 - t) \prod_{j=2}^{n} (1 - t u_j)\, dt`$.
- For $`i = 2..n`$: $`S_i(u) = T_i(u) / u_i^{n-1}`$, where
  $`T_i(u) = \int_0^1 (1 - t)(u_i - t) \prod_{j=2, j\ne i}^{n} (u_i - t u_j)\, dt`$.

  $`T_i`$ is a polynomial, and $`S_i`$ is undefined at $`u_i = 0`$.
- $`c = (d-1)/d`$. These $`S_i`$ are the normalised critical values $`P(b_i)/b_i`$ of a polynomial with $`P(0) = 0`$, $`P'(0) = 1`$ and
  critical points $`b_1 = 1`$, $`b_j = 1/u_j`$. The integrals of monomials are exact: $`\int_0^1 (1-t) t^k\, dt = 1/((k+1)(k+2))`$.

A point $`u`$ is **nondegenerate** if all $`u_j \ne 0`$.

## 2. Boxes, tasks and the tree file

A box is a closed product of real intervals $`[C_k - H_k, C_k + H_k]`$, $`k = 0..D-1`$. $`C_k`$ and $`H_k`$ are binary64 values, and
the box is the exact real set they describe.

**Tasks file** (e.g. `runs/d6_v2/d6.tasks`, text):
- The first line is a header containing `half_width=<hex float>` (for example `0x1p-2`).
- Every other line is `id c_0 … c_{D−1}` with centres as C99 hex floats (parse with `strtod`).
- Every top box has all half-widths equal to the header's half_width.

**Tree file** (`*.tree`, binary, little endian): a sequence of records `"TK"` (2 bytes) | uint32 id | uint64 nbytes
| nbytes node codes. Each record is the subdivision tree of the top box with that id, one byte per node, in
depth-first pre-order.

Node codes:

| code | meaning | claim to verify for the node's box $`B`$ |
|---|---|---|
| 16 + k | bisect coordinate $`k`$ | none: the children are $`B`$ with coordinate $`k`$ replaced by $`[C_k - H_k, C_k]`$ (centre $`C_k - H_k/2`$, half-width $`H_k/2`$) and $`[C_k, C_k + H_k]`$; the LOWER child's subtree comes first, then the upper child's |
| 0 (F) | leaf | there is an $`i`$ with $`\lvert S_i\rvert \lt  c`$ at every nondegenerate point of $`B`$, uniformly: $`\sup \lt  c`$ |
| 1 (E) | leaf | there are $`i \ne j`$ with $`\log\lvert S_i\rvert - \log\lvert S_j\rvert`$ bounded away from 0, with constant sign, at the nondegenerate points of $`B`$ where both are defined; hence $`B`$ contains no nondegenerate point with $`\lvert S_1\rvert = \dots = \lvert S_n\rvert`$ |
| 2 (L) | leaf | $`\sum_{i=1}^{n} \log\lvert S_i\rvert \le n\cdot\log c - g`$ for some $`g \gt  0`$ at every nondegenerate point of $`B`$ |
| 3 (outside) | leaf | for some variable $`v`$, $`\lvert u_v\rvert \gt  1`$ at every point of $`B`$ |
| 4 (symmetry) | leaf | for some $`v`$, $`\lvert u_v\rvert \lt  \lvert u_{v+1}\rvert`$ at every point of $`B`$ (consecutive variables in the order $`u_2, u_3, \dots`$); or $`\mathrm{Im}\, u_2 \lt  0`$ at every point of $`B`$, i.e. $`C_1 + H_1 \lt  0`$ |
| 5 (excluded) | leaf | $`B \subset`$ the closed Euclidean ball in $`\mathbb{R}^D`$ of radius exactly 1/20 around an equality point $`p`$ |
| 6 (unresolved) | leaf | must not occur; report it as a failure |

**Equality points.** The points $`p`$ are, with $`\omega = e^{2\pi i/n}`$, all tuples $`(u_2, \dots, u_n) = (\mathrm{conj}(\omega^{a_2}), \dots, \mathrm{conj}(\omega^{a_n}))`$
where $`(a_2, \dots, a_n)`$ is a permutation of $`(1, \dots, n-1)`$. There are $`(n-1)!`$ of them, and their coordinates are to be
enclosed in arb (for example with `arb_cos_pi` / `arb_sin_pi` of rationals).

A leaf's claim may be verified by ANY sound method: the checker does not have to use the test the original program
used. It only has to prove the stated claim for that leaf. A different-but-sound test is welcome.

**Refinement.** If the checker cannot prove a leaf's claim directly, it may bisect that leaf itself, up to a fixed
extra depth (e.g. 12), and prove, for every sub-box, EITHER the same claim OR any claim in the table (F/E/L/outside/
symmetry/excluded). This is sound because the leaves' claims only need to hold piecewise. It must count such
refinements, and report a failure only if refinement also fails.

**Coverage.** The checker must verify that each record's node stream parses as exactly one complete tree (no bytes
left over, none missing). It must also check that every task id in the tasks file has exactly one record. The union of
the tasks' top boxes covering the chart is checked separately by `c/check_done.py` and is not the checker's job.

## 3. Suggested method (the checker may do better)

Plain interval evaluation of $`S_i`$ over a whole box overestimates too much. The following centred forms work.

1. **Exact expansion at the centre.** Let $`m`$ be the box centre (exact doubles; $`u_v = m_v + \delta_v`$). Expand $`T_i(m + \delta)`$
   exactly as a polynomial in $`\delta`$ (and $`S_1`$ likewise), with coefficients as acb balls. Multiply out the linear factors as
   polynomials in $`(t, \delta)`$, then integrate each power of $`t`$ with the exact moments above. $`T_i`$ has degree $`\le n-1`$ in $`\delta_i`$ and
   degree $`\le 1`$ in each other $`\delta_j`$.
2. **Disc bounds.** On the box, $`\lvert\delta_v\rvert \le \rho_v := \sqrt{H_{2v}^2 + H_{2v+1}^2}`$. Then
   $`\lvert T_i(m+\delta) - c_0 - \sum_{\deg 1} c_\alpha \delta^\alpha\rvert \le R := \sum_{\deg \ge 2} \lvert c_\alpha\rvert \rho^\alpha`$.
3. **Interval bounds for $`\log\lvert S_i\rvert`$.**
   - $`\lvert T_i\rvert \in [\lvert c_0\rvert - L - R, \lvert c_0\rvert + L + R]`$ with $`L = \sum_{\deg 1} \lvert c_\alpha\rvert \rho^\alpha`$.
   - $`\lvert u_i\rvert \in [\lvert m_i\rvert - \rho_i, \lvert m_i\rvert + \rho_i]`$.
   - Then $`\log\lvert S_i\rvert = \log\lvert T_i\rvert - (n-1) \log\lvert u_i\rvert`$.
4. **Affine models.** If $`q := (L + R)/\lvert c_0\rvert \lt  1/2`$, write $`T = c_0(1 + w)`$ with $`\lvert w\rvert \le q`$. Then
   $`\log\lvert T\rvert = \log\lvert c_0\rvert + \mathrm{Re}\left(\sum \gamma_\alpha \delta^\alpha\right) \pm \left(R/\lvert c_0\rvert + q^2/(2(1-q))\right)`$, with $`\gamma_\alpha = c_\alpha/c_0`$ on the degree-1 terms. Likewise, if
   $`p := \rho_i/\lvert m_i\rvert \lt  1/2`$, then $`\log\lvert u_i\rvert = \log\lvert m_i\rvert + \mathrm{Re}(\delta_i/m_i) \pm p^2/(2(1-p))`$. This gives
   $`\log\lvert S_i\rvert = a_i + G_i\cdot(x - x_m) \pm e_i`$, where $`G_i`$ is a real gradient vector in $`\mathbb{R}^D`$. Over the box,
   $`\lvert G_i\cdot(x - x_m)\rvert \le \sum_k \lvert G_{ik}\rvert H_k`$.
5. **The tests.**
   - F: some $`a_i + \sum_k \lvert G_{ik}\rvert H_k + e_i \lt  \log c`$, or an interval upper bound of $`\log\lvert S_i\rvert \lt  \log c`$.
   - E: some $`\lvert a_i - a_j\rvert \gt  e_i + e_j + \sum_k \lvert G_{ik} - G_{jk}\rvert H_k`$, or an interval lower bound of one $`\log\lvert S_i\rvert`$ above an
     interval upper bound of another $`\log\lvert S_j\rvert`$, with $`i \ne j`$.
   - L: $`\sum_i a_i + \sum_k \left\lvert\sum_i G_{ik}\right\rvert H_k + \sum_i e_i \lt  n \log c`$.
   All comparisons are done in arb (`arb_lt` etc.; a comparison that is not certain counts as failed).
6. **Degenerate points.** Near $`u_i = 0`$, $`S_i`$ is undefined, and $`\lvert S_i\rvert \to \infty`$ as $`u_i \to 0`$ with $`T_i \ne 0`$. A lower bound
   $`\log\lvert S_i\rvert \ge \log\inf\lvert T_i\rvert - (n-1) \log\sup\lvert u_i\rvert`$ is valid at the nondegenerate points even if the box meets $`u_i = 0`$.
   An upper bound for $`\log\lvert S_i\rvert`$ needs $`\lvert u_i\rvert`$ bounded away from 0 on the box.

## 4. Engineering requirements

- C with FLINT $`\ge 3`$ (arb/acb); Homebrew FLINT 3.6.0 is installed (`/opt/homebrew/include/flint`). Use about 64–128 bit
  precision.
- Preallocate every arb/acb variable per thread. No allocation in the per-leaf inner loop.
- Use pthreads over task records, with a thread count from the command line. Support sharding: process only records
  with `id % NSHARD == SHARD` (and optionally an id range), so that parts can run on another machine.
- Stream the tree file (it can be ~2 GB); do not load it whole.
- Output one line per task: id, leaves by type, leaves verified directly, leaves verified after refinement,
  sub-boxes created by refinement, failures. Then a final summary line. The exit status is non-zero if any failure
  occurred or any id is missing.
- Also provide a `--node-limit` / `--sample` mode for benchmarking, e.g. check only the first N records, or a
  random 1/K of the records.
