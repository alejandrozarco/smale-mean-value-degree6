# Smale's mean value conjecture, degree 6: computation and certificates

Status: **computational results, not peer reviewed.** First published 2026-09-30.

For a complex polynomial P of degree d and a point z with P′(z) ≠ 0, Smale's mean value conjecture concerns the
smallest constant K such that some critical point c of P satisfies |P(z) − P(c)| ≤ K |z − c| |P′(z)|. The
conjectured sharp value is K = (d−1)/d, attained by z − z^d/d at z = 0.

This repository contains a computation for d = 6 (constant 5/6), intended to show that the sharp value holds in
degree 6. It has three parts:
1. an equal-modulus reduction (`REDUCTION.md`);
2. a branch-and-bound cover of a compact 8-dimensional chart in rigorous floating-point ball arithmetic (`c/`,
   `runs/d6_v2/`);
3. an interval-arithmetic second-order certificate near the extremal configuration (`local/`).

The same code applied to d = 5 is included (`runs/d5_v2/`).

## Normalisation

Take z = 0, P(0) = 0 and P′(0) = 1, with critical points b_1, …, b_n (n = d − 1). The quantity to bound is

F(B) = min_i |S_i|, where S_i = P(b_i)/b_i = ∫₀¹ Π_j (1 − t b_i/b_j) dt.

The chart: b_1 = 1 is a critical point of minimal modulus, u_j = 1/b_j and |u_j| ≤ 1, reduced by symmetry to
|u_2| ≥ … ≥ |u_n| and Im u_2 ≥ 0. The chart has real dimension 2(d − 2) = 8. For i ≥ 2,
S_i = T_i/u_i^{n−1} with T_i = ∫₀¹ (1−t)(u_i − t) Π_{j≠1,i} (u_i − t u_j) dt, a polynomial. The extremal
configurations are u = (ω̄^{a_2}, …, ω̄^{a_n}) with ω = e^{2πi/n} and (a_j) a permutation of 1..n−1; there are 24
for d = 6.

## Method (outline)

- **Reduction** (`REDUCTION.md`). Some point of the closure of the image of (S_1, …, S_n) attains sup F with all
  |S_i| equal. The proof follows Crane and Ng and uses convexity of amoeba complements; it does not assume that a
  maximiser exists. A finite cover of the closed chart by closed boxes then suffices, provided every box either lies
  in a small ball around an equality point, where the local certificate applies, or passes one of the following
  tests uniformly on the box. Here c = (d−1)/d.
  - **F:** some |S_i| < c.
  - **E:** some |S_i| ≠ |S_j|, so the box misses E = {|S_1| = … = |S_n|}.
  - **L:** (1/n) Σ log|S_i| < log c.
- **Branch and bound** (`c/smale_bb_v2.c`).
  - Each T_i is expanded exactly at the box centre, with coefficients as complex balls. This gives interval bounds
    and first-order affine models of log|S_i|, from which the F, E and L tests are evaluated.
  - The start grid is 4^8 boxes of half-width 1/4. Boxes that pass no test are bisected.
  - A box is set aside only if it lies inside the closed Euclidean u-ball of radius 0.05 around an equality point.
  - The arithmetic is binary64 round-to-nearest ball arithmetic with an explicit error term for every operation, and
    the correctly rounded logarithm of CORE-MATH (`third_party/core-math-log/`). The error analysis is
    `c/ERROR_ANALYSIS_v2.md`.
  - Each finished top-level task is written as one checksummed record, bound to a hash of the source and the
    configuration. `c/check_done.py` checks the records and the grid coverage independently.
- **Local certificate** (`local/LOCAL_CERT.md`). On E, within |ε|₂ ≤ 0.0527 of the equality point in the chart
  b_j = ω^{j−1} e^{ε_j}, it shows log F ≤ log(5/6) − 0.00757 |ε|₂². It uses Taylor models of degree 6 in arb ball
  arithmetic, and exact arithmetic in Q(ω) for the zeroth- and first-order terms. The u-balls of radius 0.05 map into
  |ε|₂ ≤ 0.05/0.95 = 0.05263.

## Results

| | d = 5 | d = 6 |
|---|---|---|
| top-level tasks | 3 072 | 49 152 |
| boxes processed | 2 490 060 | 1 874 817 938 |
| leaves closed by F / E / L | 757 590 / 143 562 / 100 008 | 561 202 241 / 145 383 090 / 37 946 415 |
| leaves outside chart / symmetry / excluded ball | 37 156 / 205 334 / 2 916 | 18 661 249 / 173 649 477 / 591 073 |
| unresolved boxes | 0 | 0 |
| maximum subdivision depth | 47 | 71 |
| task CPU time | 7 s | 23 735 s |
| `c/check_done.py` | CERTIFICATE COMPLETE (`runs/d5_v2/CHECK.txt`) | CERTIFICATE COMPLETE (`runs/d6_v2/CHECK.txt`) |
| local certificate bound g(T), T = 0.0527 | ≤ −0.02284 | ≤ −0.0075741 |

Source sha256 `8ddcec1767e0b64b404843cc0e0d99684b6526ef9c6d34288129a89fe0b9487f` (`c/smale_bb_v2.c`). Binary
sha256 `e99fdb885004db53961c411065372595540137ca3f78b654b0da96019aa08723` (Apple clang 21.0.0, arm64; see
`runs/d6_v2/BUILD.txt`). `runs/d6_v2/d6.done` sha256
`0f8199fb3b360623e1ee7e493bfb530da7868483340f85222caaa96a2db0506b`.

## Figures

Regenerate with `python3 figures/make_figures.py` (matplotlib).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/outcomes_dark.svg">
  <img alt="d = 6: leaf boxes by outcome (F, E, L tests, symmetry, outside chart, excluded ball, unresolved)" src="figures/outcomes_light.svg">
</picture>

Leaf boxes of the d = 6 subdivision by the test or rule that closed them (`runs/d6_v2/d6.done`).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/tasks_dark.svg">
  <img alt="Histograms of the number of boxes processed per top-level task, d = 5 and d = 6" src="figures/tasks_light.svg">
</picture>

Boxes processed per top-level task, for d = 5 and d = 6.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/distance_dark.svg">
  <img alt="d = 6: boxes processed per task against the distance of the task centre to the nearest equality point" src="figures/distance_light.svg">
</picture>

d = 6: boxes processed per task, against the distance from the task centre to the nearest equality point
(`runs/d6_v2/d6.done`, `runs/d6_v2/d6.tasks`).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/local_dark.svg">
  <img alt="Deficit (log c − log F)/|x|^2 on E near the equality point: numerical samples and certified lower bounds" src="figures/local_light.svg">
</picture>

(log c − log F)/|x|² on E near the equality point (d = 6). The dots are float64 samples (`runs/e_profile_d6.tsv`),
and the line is their minimum at each radius.
Each segment is a certified lower bound, valid for |x| ≤ T (`runs/local_cert_v3_d6_T*.log`). The dashed line marks
the radius that the excluded balls map into.

## Contents

| path | content |
|---|---|
| `REDUCTION.md` | equal-modulus reduction (Proposition R) and the form in which the computation uses it |
| `c/smale_bb_v2.c` | branch-and-bound evaluator (its sha256 is bound into every run record) |
| `c/build_v2.sh` | build command; writes `BUILD.txt` (compiler, hashes, fused-multiply-add scan) |
| `c/ERROR_ANALYSIS_v2.md` | per-operation floating-point error analysis of the evaluator |
| `c/check_done.py` | independent completion checker (exact-rational grid, records, configuration) |
| `c/progress.py` | progress of a running job |
| `c/regress/` | task replay, equality-point certification (arb), primitive checks, targeted regression |
| `third_party/core-math-log/` | CORE-MATH correctly rounded log (MIT licence, upstream commit in `UPSTREAM_COMMIT`) |
| `runs/d6_v2/`, `runs/d5_v2/` | build records, task lists, per-task records (`*.done`), unresolved lists (empty), checker output |
| `local/` | local certificate (`local_cert.py`, `exact_facts.py`, `LOCAL_CERT.md`), Jacobian rank check, numerical checks |
| `runs/local_cert_v3_*.log` | local certificate outputs |
| `runs/e_profile_d6.tsv` | numerical samples on E (figure data) |
| `runs/regress_v2/` | outputs of the regression scripts |
| `src/smale.py` | float64 evaluation of S_i (used by the numerical checks) |
| `figures/` | figures and the script that makes them |
| `MANIFEST.sha256` | sha256 of every file |

## Reproduce

Requirements: a C compiler; Python 3 with `mpmath`, `numpy`, `sympy`, `python-flint` (arb) and `matplotlib`.
Versions used: Apple clang 21.0.0 on macOS arm64 (Apple M5), Python 3.9.6, python-flint 0.6.0, mpmath 1.3.0,
numpy 1.26.4, sympy 1.14.0, matplotlib 3.9.4.

```sh
# build (writes OUT/smale_bb_v2.bin and OUT/BUILD.txt)
c/build_v2.sh OUT

# check the published run: records, configuration, source hash, exact grid coverage
python3 c/check_done.py runs/d6_v2/d6 --build runs/d6_v2/BUILD.txt
python3 c/check_done.py runs/d5_v2/d5 --build runs/d5_v2/BUILD.txt
# with an identical toolchain the binary hash also matches: add --bin OUT/smale_bb_v2.bin

# recompute a random sample of d = 6 tasks and compare with the records (add task ids to include specific ones)
python3 c/regress/replay_tasks.py OUT/smale_bb_v2.bin runs/d6_v2/d6 50

# full rerun (about 24 000 CPU-seconds for d = 6; resumable), then check
OUT/smale_bb_v2.bin run 6 0.05 4 3 OUT/d6 2 1
python3 c/check_done.py OUT/d6 --build OUT/BUILD.txt --bin OUT/smale_bb_v2.bin

# equality points (certified in arb), primitives, targeted regression
python3 c/regress/eqpoints_v2.py OUT/smale_bb_v2.bin
cc -O2 -ffp-contract=off -DSRC_SHA256_RAW=0 -Ithird_party/core-math-log -o OUT/prim_test c/regress/prim_test_v2.c third_party/core-math-log/log.c -lpthread
OUT/prim_test | python3 c/regress/prim_check_v2.py
python3 c/regress/regress_v2.py OUT/smale_bb_v2.bin runs/d5_v2/d5 OUT/regress

# local certificate (about 1 minute each), exact facts, Jacobian rank
python3 local/local_cert.py 6 0.0527
python3 local/local_cert.py 5 0.0527
python3 local/exact_facts.py
python3 local/rank_check.py 6

# figures
python3 figures/make_figures.py
```

## Assumptions

The arithmetic argument (`c/ERROR_ANALYSIS_v2.md` §A.1) assumes:
- IEEE-754 binary64 operations, correctly rounded to nearest with gradual underflow;
- a compiler that neither contracts nor reassociates floating-point expressions;
- correct rounding of the compiled CORE-MATH `cr_log`.

The local certificate assumes the correctness of arb (python-flint).

## References

- S. Smale, The fundamental theorem of algebra and complexity theory, Bull. Amer. Math. Soc. 4 (1981) 1–36.
- D. Tischler, Critical points and values of complex polynomials, J. Complexity 5 (1989) 438–456.
- E. Crane, Topics in conformal geometry and dynamics, PhD thesis, Cambridge, 2003.
- E. Crane, Extremal polynomials in Smale's mean value conjecture, Comput. Methods Funct. Theory 6 (2006) 145–163.
- E. Crane, A bound for Smale's mean value conjecture for complex polynomials, Bull. London Math. Soc. 39 (2007)
  781–791.
- A. Conte, E. Fujikawa, N. Lakic, Smale's mean value conjecture and the coefficients of univalent functions,
  Proc. Amer. Math. Soc. 135 (2007) 3295–3301.
- T. W. Ng, Smale's mean value conjecture and related problems, survey talk, 2018.
- Z. Liu, Smale's mean value conjecture in degree six: certified enclosures and a cancellation obstruction,
  Zenodo record 22390113.
- A. Jatar, T. W. Ng, Smale's mean value conjecture and its dual conjecture for complex polynomials,
  arXiv:2608.27047.
- https://github.com/tadamcz/mean-value-problem (Lean formalisation concerning the conjecture with K = 1 in large
  degree; not examined here).
- I. M. Gelfand, M. M. Kapranov, A. V. Zelevinsky, Discriminants, resultants, and multidimensional determinants,
  Birkhäuser, 1994.
- M. Forsberg, M. Passare, A. Tsikh, Laurent determinants and arrangements of hyperplane amoebas, Adv. Math. 151
  (2000) 45–70.
- A. Sibidanov, P. Zimmermann, S. Glondu, The CORE-MATH project, ARITH 2022.
