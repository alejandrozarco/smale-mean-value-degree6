# arbcheck: an independent FLINT/arb check of the subdivision tree

`arbcheck` checks a subdivision tree (`*.tree`) against its tasks file (`*.tasks`) for Smale's mean
value conjecture in degree $`d`$ ($`4 \le d \le 7`$). It was written from `SPEC.md` and from the mathematics
only. Every claim is decided in FLINT/arb ball arithmetic.

## Build

```
make                       # or:
cc -O2 -Wall -ffp-contract=off -fcx-limited-range -I/opt/homebrew/include -L/opt/homebrew/lib \
   -o arbcheck arbcheck.c -lflint -lgmp -lmpfr -lpthread
make selftest              # compares the expansions with direct evaluation for d = 4..7
```

Requirements: FLINT $`\ge 3`$ (tested with 3.6.0), GMP, MPFR and pthreads. The compiler flag
`-fcx-limited-range` only affects the double-precision predictor, which never accepts a claim (see
below). On Linux, set `FLINT_PREFIX=/usr` (or wherever FLINT is installed).

## Usage

```
arbcheck [options] d tasks_file tree_file

  -t N          threads (default 1); records are distributed over threads, largest first
  --shard S/N   only records with id % N == S
  --ids A-B     only records with A <= id <= B
  --sample K    only a pseudo-random 1/K of the records (hash of id; --seed X changes it)
  --limit N     only the first N selected records in file order
  --largest N   also check the N largest available records
  --prec P      arb precision in bits (default 64)
  --maxref R    extra refinement depth per leaf (default 12; 0 = strict, no refinement)
  --timing      time every leaf and report ns/leaf and leaves/s/thread by leaf type
  --quiet       print only failing task lines and the summary
  --selftest    run the self test (formulas, bounds, equality points) for degree d and exit
```

Examples:

```
./arbcheck -t 8 5 runs/d5_v2/d5.tasks OUT/d5.tree
./arbcheck -t 8 --shard 0/4 6 runs/d6_v2/d6.tasks OUT/d6.tree   # 1 of 4 machines
./arbcheck -t 2 --timing --sample 74 --largest 5 6 runs/d6_v2/d6.tasks OUT/d6.tree
```

The tree file is streamed. A first pass reads only the 14-byte record headers and seeks over the
node data. Each worker thread then reads its records through its own `FILE*` and a 1 MB buffer, so
memory use is small and does not depend on the size of the file. Only complete records are used.
A record that runs past the end of the file (for example while the file is still being written)
counts as an incomplete tail and is not checked.

Exit status: 0 only if every checked leaf is verified and every coverage check passes. The status is
1 if any leaf or record fails, or if in full mode an id is missing, duplicated, unknown, the header
magic is wrong, or the tail is incomplete. It is 2 for usage or I/O errors. Sample, limit and largest
modes skip the check for missing ids and say so in the summary.

## Output

One line per checked record, printed as records finish:

```
task <id> nodes=<n> maxdepth=<d> F=.. E=.. L=.. O=.. S=.. X=.. U=.. direct=.. refined=.. subboxes=.. fail=.. parse=<status>
```

- F, E, L, O, S, X, U: leaves by code (0 F, 1 E, 2 L, 3 outside, 4 symmetry, 5 excluded, 6 unresolved).
- `direct`: leaves whose own claim was proved on the leaf box.
- `refined`: leaves proved only after the checker bisected them itself.
- `subboxes`: sub-boxes created by refinement.
- `fail`: leaves not proved, plus 1 if the record does not parse.
- `parse`: `ok`, `missing_bytes` (the tree is incomplete), `leftover_bytes` (bytes remain after one
  complete tree), `invalid_code` (a code from 7 to 15, or 16+k with $`k \ge D`$), `inexact_split` or
  `too_deep`.

Then comes one `SUMMARY` line with the totals, `failed_records`, `missing`, `dup`, `unknown`,
`badheader`, `tail_incomplete`, wall time, CPU time and `status=PASS|FAIL`. With `--timing`, one
`TIMING` line per leaf type follows. stderr gets the up to 20 first `FAIL` lines and up to 20
`REFINED` notes, each with the id, the code and the box as hex centre $`\pm`$ half-width.

## Mathematics used

Set $`u_1 := 1`$, $`n = d - 1`$ and $`\mu_k = \int_0^1 (1-t) t^k\, dt = 1/((k+1)(k+2))`$. Expanding the products in $`t`$ and
integrating term by term gives

- $`S_1 = \sum_{k=0}^{n-1} (-1)^k \mu_k e_k(u_2,\dots,u_n)`$
- $`T_i = \sum_{k=0}^{n-1} (-1)^k \mu_k u_i^{n-1-k} e_k(\{u_j : j \ne i,\ 1 \le j \le n\})`$, for $`i \ge 2`$

Here $`e_k`$ is the elementary symmetric polynomial, and $`u_1 = 1`$ is one of the variables in the second
formula, because $`(u_i - t) = (u_i - t\cdot u_1)`$.

**Taylor expansion at the box centre.** Let $`m`$ be the box centre. The coordinates $`m_v`$ are exact
doubles and are loaded exactly into acb. The checker computes $`e_r(m_W)`$ as balls for every subset $`W`$ of
the variables, by the recurrence $`e_r(W \cup \{v\}) = e_r(W) + m_v e_{r-1}(W)`$. The subsets are computed
lazily and shared by all $`S_i`$. $`S_1`$ is multi-affine, so the coefficient of $`\delta^J`$ is
$`\sum_r (-1)^{r+|J|} \mu_{r+|J|} e_r(m_{\mathrm{all}\setminus J})`$.

$`T_i`$ is affine in each $`\delta_j`$ with $`j \ne i`$. The coefficient of $`\delta^J`$ is a polynomial $`A_J(u_i)`$ whose
coefficients are $`\pm\mu\cdot E1`$, where $`E1_r(W) = e_r(W) + e_{r-1}(W) = e_r(\{1\} \cup W)`$ and $`W = O_i \setminus J`$. Its Taylor
coefficients at $`m_i`$ are

```math
c_{J,a} = \sum_{p \ge a} C(p,a) (-1)^{n-1-p} \mu_{n-1-p} m_i^{p-a} E1_{n-1-p-|J|}(W)
```

Each coefficient is evaluated as a real `arb_dot` over interleaved real and imaginary parts, from a
per-variable vector and a per-subset vector that are both precomputed per box. There is no acb_poly
code and no allocation in the per-leaf loop: all arb, acb and mag variables live in a per-thread
workspace.

**Bounds.** With $`\rho_v = \sqrt{H_{2v}^2 + H_{2v+1}^2}`$, rounded up, we have $`\lvert\delta_v\rvert \le \rho_v`$ on the box. The
checker computes, all as upper bounds in mag arithmetic:

- $`L = \sum_{\mathrm{order}\ 1} \lvert c_\alpha\rvert \rho^\alpha`$
- $`R = \sum_{\mathrm{order} \ge 2} \lvert c_\alpha\rvert \rho^\alpha`$
- $`\lvert T\rvert \in [\lvert c_0\rvert - L - R, \lvert c_0\rvert + L + R]`$
- $`\lvert u_i\rvert \in [\lvert m_i\rvert - \rho_i, \lvert m_i\rvert + \rho_i]`$

These give $`S_{\mathrm{lo}} \le \lvert S_f\rvert \le S_{\mathrm{hi}}`$ at every nondegenerate point. $`S_{\mathrm{hi}} = \infty`$ when $`\lvert u_i\rvert`$ is not bounded away
from 0. The lower bound remains valid even when the box meets $`u_i = 0`$.

**Affine model.** The model exists when $`q = (L+R)/\lvert c_0\rvert \lt  1`$ and, for $`i \ge 2`$, $`p = \rho_i/\lvert m_i\rvert \lt  1`$.
Then, for every $`x`$ in the box,

```math
\log\lvert S_f(x)\rvert \in a + G\cdot(x - m) \pm e
```

with:

- $`a = \log\lvert c_0\rvert - (n-1) \log\lvert m_i\rvert`$
- $`G`$ from $`\mathrm{Re}(c_{\mathrm{lin}} \delta / c_0) - (n-1) \mathrm{Re}(\delta_i/m_i)`$
- $`e = R/\lvert c_0\rvert + q^2/(2(1-q)) + (n-1) p^2/(2(1-p))`$

The last two terms use $`\lvert \log(1+w) - w\rvert \le q^2/(2(1-q))`$ for $`\lvert w\rvert \le q \lt  1`$. Over the box,
$`\lvert G\cdot(x-m)\rvert \le \sum_k \lvert G_k\rvert H_k`$.

## Tests per leaf type

A leaf is verified by any one of the listed tests. All comparisons are certain arb or mag comparisons.

| code | claim | test(s) |
|---|---|---|
| F (0) | some $`S_i`$ with $`\sup\lvert S_i\rvert \lt  c`$ | $`S_{\mathrm{hi}} \lt  c`$ (mag); or $`a + e + \sum\lvert G_k\rvert H_k \lt  \log c`$ (arb_lt) |
| E (1) | some $`i \ne j`$ with $`\log\lvert S_i\rvert - \log\lvert S_j\rvert`$ of constant sign, bounded away from 0 | $`S_{\mathrm{lo},i} \gt  S_{\mathrm{hi},j}`$ (mag); or the ball $`(a_i - a_j) \pm (e_i + e_j + \sum\lvert G_{ik} - G_{jk}\rvert H_k)`$ excludes 0 |
| L (2) | $`\sum \log\lvert S_i\rvert \lt  n \log c - g`$ | $`\prod S_{\mathrm{hi}} \lt  c^n`$ (mag); or $`\sum`$ over $`i`$ of the affine model where it exists (the gradients are summed before taking $`\lvert\cdot\rvert`$), or $`\log S_{\mathrm{hi},i}`$ where it does not, is $`\lt  n \log c`$ (arb_lt) |
| outside (3) | some $`\lvert u_v\rvert \gt  1`$ on $`B`$ | $`\inf \lvert u_v\rvert^2 = \mathrm{lo}(x)^2 + \mathrm{lo}(y)^2 \gt  1`$, with $`\mathrm{lo} = \max(0, \lvert C\rvert - H)`$ |
| symmetry (4) | $`C_1 + H_1 \lt  0`$, or some $`\lvert u_v\rvert \lt  \lvert u_{v+1}\rvert`$ on $`B`$ | the sum is negative in arb; or $`\sup \lvert u_v\rvert^2 \lt  \inf \lvert u_{v+1}\rvert^2`$ |
| excluded (5) | $`B \subset`$ closed ball of radius 1/20 about an equality point | $`\sum_k (\lvert C_k - p_k\rvert + H_k)^2 \le 1/400`$, with $`p_k = \cos/\sin(2\pi a/n)`$ as balls (`arb_sin_cos_pi_fmpq`); all $`(n-1)!`$ points are tried |
| unresolved (6) | none | always a failure |

**Which functions and pairs are tried.** A double-precision copy of the same model predicts which $`S_i`$
(for F) or which pair (for E) will succeed, and those are tried first. The prediction orders the
attempts, and in `test_F`, `test_L` and the exclusion test it can also skip an attempt that is predicted to fail.
Every acceptance is decided in arb, and the others are tried if a candidate fails. Skipping can only cause a
rejection (followed by refinement), never an acceptance. On the $`d=6`$
sample, F needs exactly 1 arb model per leaf and E exactly 2.

**Refinement (§2).** If a leaf's claim is not proved directly, the leaf is bisected along its widest
coordinate, up to `--maxref` levels. Each sub-box must satisfy the leaf's own claim or any other claim
of the table (outside, symmetry, excluded, F, E, L). Bisection is exact, which is checked with TwoSum.

**Coverage (§2).**
- Each record must parse as exactly one complete pre-order tree: no bytes missing, none left over,
  and valid codes only.
- Every task id must have exactly one record: missing, duplicate and unknown ids are counted.
- Every top box has half-width `half_width` from the tasks header.

**Self test.**
- Compares the centre values with a direct double-precision evaluation of the defining integrals.
- Checks $`S_{\mathrm{lo}} \le \lvert S_f\rvert \le S_{\mathrm{hi}}`$ and the affine enclosure at random points of random boxes.
- Checks that all $`(n-1)!`$ equality points satisfy $`\lvert S_f(p)\rvert = c`$ for every $`f`$ and are distinct.

## Results

The recorded runs are `runs/d5_arbcheck/arbcheck.log` and `runs/d6_arbcheck/arbcheck.log`, with versions and hashes in
`runs/d6_arbcheck/BUILD.txt`.

- $`d=5`$: status PASS; 1,246,566 leaves, all proved directly.
- $`d=6`$: status PASS; 937,433,545 leaves; 118 proved after further bisection (240 sub-boxes in total); 0 failures and
  0 missing ids.

During development the checker was also run on corrupted copies of the $`d=5`$ tree. The corruptions were relabelled
leaves, truncated records, extra bytes, a duplicated record and a missing record, and every one was reported with
`--maxref 0`. With the default `--maxref 12`, relabelled leaves can be re-proved by other claims of the table, which
§2 of the specification allows. These tests are not recorded in this repository.

**Reading the refinement counts.** A clean tree verifies almost entirely directly, so a nonzero
`refined` count deserves a look, and `--maxref 0` checks each leaf's own claim strictly.

**Input checks.** Non-finite centres or half-widths in the tasks file are rejected. The sampling modes (`--sample`,
`--limit`, `--largest`) check only part of a tree and exit 0 even though ids are missing, so they are for
benchmarking only. A certificate run checks all records and requires `status=PASS` with `missing=0`.
