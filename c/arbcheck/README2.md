# arbcheck2: arbcheck with a second model of $`S_i`$ (the w-form)

`arbcheck2` is `arbcheck` (see `README.md`) plus one more way to prove leaf claims: the w-form of
`SPEC_W.md`. It was written from `SPEC.md`, `SPEC_W.md` and the mathematics only, starting from a copy of
`arbcheck.c`. `arbcheck.c` is unchanged (sha256 `68b18fda…`). `arbcheck2.c` and this file were written by an
AI model; for reviews see the section Review below.

Nothing else changes: the claims, the file formats, the coverage checks, refinement, sharding and sampling
are those of `arbcheck`. Every claim is still decided by certain comparisons in FLINT/arb ball arithmetic.

## What is added

For $`i = 2..n`$, with $`w = 1/u_i`$ and $`O_i = \{2..n\} \setminus \{i\}`$,

```math
S_i = \int_0^1 (1-t)(1-tw) \prod_{j \in O_i} (1 - t u_j w)\, dt = \sum_{k=0}^{n-1} (-1)^k \mu_k\, w^k\, e_k(\{1\} \cup \{u_j : j \in O_i\}),
```

with $`\mu_k = 1/((k+1)(k+2))`$. This is $`T_i / u_i^{n-1}`$: in `arbcheck`,
$`T_i = \sum_k (-1)^k \mu_k u_i^{n-1-k} e_k(\{1\} \cup \{u_j\})`$, so the coefficient of $`u_i^{n-1-k}`$ in $`T_i`$ is the
coefficient of $`w^k`$ here. The self-test also checks the identity numerically (see below).

The w-form models $`\log\lvert S_i\rvert`$ directly. It has no $`(n-1)\log\lvert u_i\rvert`$ term, whose variation made the
T-form too weak on large leaves.

**Two passes per claim.** Each F, E or L claim is first tried with the T-form only. This pass runs exactly
the tests of `arbcheck`. Only if it fails is the w-form pass run, which tries every test that uses at least
one w-form quantity:

- **F:** $`\sup\lvert S_i\rvert \lt  c`$ from the w-form bound, or the w-form affine upper bound $`\lt  \log c`$.
- **E:** magnitude bounds that use, for each function, the better of the T-form and w-form bounds. Then the
  affine pairs (T, w), (w, T) and (w, w).
- **L:** for each $`f`$, one of: the T-form affine model, the w-form affine model, or $`\log`$ of the best $`\sup\lvert S_f\rvert`$.
  A double-precision estimate chooses the combination (at most $`3^n`$ of them). The chosen combination is then
  decided in arb. Any combination is valid, because each term is a valid upper bound on its own.

During refinement, a sub-box first gets the T-form pass for all claims (O, S, X, F, E, L). Only after that
does it get the w-form pass for F, E and L.

So anything `arbcheck` proves, `arbcheck2` proves the same way. `--no-wform` turns the second pass off; the
output is then the same as `arbcheck`'s, apart from the version line and the added counters (task lines and
FAIL/REFINED lines are identical).

**Counting.**
- `wdirect` (per task) and `wform_direct` (summary, also split by F/E/L) count the leaves whose own claim
  failed in the T-form pass and was proved in the w-form pass.
- `wsub` / `wform_subboxes` count refinement sub-boxes proved the same way.

## How the w-form is made rigorous

Let $`m`$ be the box centre, loaded exactly into acb, and $`\rho_v`$ the upper bound of
$`\sqrt{H_{2v}^2 + H_{2v+1}^2}`$ from `arbcheck` (`mag_hypot`, rounded up), so $`\lvert\delta_v\rvert \le \rho_v`$ on the box.
`mlo` is a certified lower bound of $`\lvert m_i\rvert`$ (`acb_get_mag_lower`). Each item below gives the quantity, how
it is computed, and why the result is valid.

1. **Domain.** The w-form for $`S_i`$ is used only if `mag_cmp(mlo, rho) > 0`, that is $`\lvert m_i\rvert \ge`$ `mlo` $`\gt `$ `rho` $`\ge \lvert\delta_i\rvert`$.
   Then $`\lvert u_i\rvert \ge \lvert m_i\rvert - \lvert\delta_i\rvert \gt  0`$ on the whole box, so $`S_i`$ is defined everywhere on it. Boxes that
   meet $`\lvert u_i\rvert \le \rho_i`$ never use the w-form for $`S_i`$.
2. **$`1/m_i`$.** `acb_inv(w0, m)` gives a ball that contains the exact $`1/m_i`$ ($`m_i`$ is exact and nonzero).
   `acb_sqr(w0)` gives a ball that contains $`1/m_i^2`$. The expansion point is the exact $`1/m_i`$, and the balls
   only enclose it.
3. **$`\rho_w`$.** $`\lvert\delta_w\rvert = \lvert 1/u_i - 1/m_i\rvert = \lvert\delta_i\rvert / (\lvert m_i\rvert\lvert u_i\rvert) \le \rho_i / (\lvert m_i\rvert(\lvert m_i\rvert - \rho_i))`$.
   The right side decreases in $`\lvert m_i\rvert`$ for $`\lvert m_i\rvert \gt  \rho_i`$, so replacing $`\lvert m_i\rvert`$ by `mlo` keeps it an upper
   bound. In mag arithmetic: `mag_sub_lower(mlo, rho)` and `mag_mul_lower` give a lower bound of the
   denominator, and `mag_div(rho, ·)` an upper bound of the quotient. Zero or infinite intermediates make the
   w-form unavailable on that box.
4. **The η bound.** $`\delta_w = -\delta_i/m_i^2 + \eta`$ with $`\eta = \delta_i^2 / (m_i^2 u_i)`$. So
   $`\lvert\eta\rvert \le \rho_i^2 / (\lvert m_i\rvert^2(\lvert m_i\rvert - \rho_i))`$. This is computed like $`\rho_w`$, from lower bounds of the
   denominator and the upper bound $`\rho_i^2`$ of the numerator.
5. **Coefficients.** The coefficient of $`\delta^J`$ ($`J \subseteq O_i`$) is $`B_J(w) = \sum_k (-1)^k \mu_k E1(W, k-\lvert J\rvert)\, w^k`$.
   Here $`W = O_i \setminus J`$ and $`E1(W, r) = e_r(\{1\} \cup m_W)`$. These are the same per-subset symmetric functions
   that `arbcheck` already uses, held as balls. The Taylor coefficients at $`w_0`$ are
   $`c_{J,a} = \sum_{k \ge \max(a, \lvert J\rvert)} \binom{k}{a} (-1)^k \mu_k\, w_0^{k-a} E1(W, k-\lvert J\rvert)`$, computed with `acb_dot`.
   The constants $`\binom{k}{a}(-1)^k \mu_k`$ are exact rationals evaluated as arb balls. Each ball contains the exact
   coefficient at the exact $`1/m_i`$, because it is computed from balls that contain the exact inputs.
6. **Disc bounds.**
   - $`L = \sum_{\text{order } 1} \lvert c\rvert \rho^\alpha`$ and $`R = \sum_{\text{order} \ge 2} \lvert c\rvert \rho^\alpha`$, using $`\rho_w`$ for $`\delta_w`$. Both use
     `acb_get_mag` (an upper bound over the ball) and mag products, which round up.
   - $`\lvert S_i\rvert \le \lvert c_0\rvert_{\text{hi}} + L + R`$, and $`\lvert S_i\rvert \ge \lvert c_0\rvert_{\text{lo}} - L - R`$ (`mag_sub_lower`), at every point of the box.
7. **Affine model.** If $`c_{0,\text{lo}} \gt  0`$ and $`q = (L+R)/c_{0,\text{lo}} \lt  1`$ (an upper bound in mag):
   - $`\log\lvert S_i\rvert = \log\lvert c_0\rvert + \mathrm{Re}(z) \pm q^2/(2(1-q))`$ with $`S_i = c_0(1+z)`$ and $`\lvert z\rvert \le q`$. This uses
     $`\lvert\log(1+z) - z\rvert \le \sum_{k \ge 2} q^k/k \le q^2/(2(1-q))`$.
   - $`\mathrm{Re}(z) = \mathrm{Re}(\gamma_w \delta_w + \sum_j \gamma_j \delta_j) \pm R/\lvert c_0\rvert`$, with $`\gamma = c_{\text{lin}}/c_0`$ (`acb_div`, which encloses the
     exact quotient).
   - $`\mathrm{Re}(\gamma_w \delta_w) = \mathrm{Re}(-\gamma_w w_0^2\, \delta_i) \pm \lvert\gamma_w\rvert_{\text{hi}} \cdot \lvert\eta\rvert`$-bound.
   - So $`a = \log\lvert c_0\rvert`$ (`arb_log_hypot`), the gradient is $`(\mathrm{Re}\,g, -\mathrm{Im}\,g)`$ per variable with $`g = \gamma_j`$ or
     $`g = -\gamma_w w_0^2`$, and $`e = R/c_{0,\text{lo}} + q^2/(2(1-q)) + \lvert\gamma_w\rvert_{\text{hi}}\, \eta_{\text{bound}}`$, all rounded up.
   - `SPEC_W.md` suggests $`q \lt  1/2`$. The bound holds for every $`q \lt  1`$, and `arbcheck` already uses $`q \lt  1`$.
     A larger $`q`$ only gives a larger $`e`$.
8. **Tests.** These are the same certain comparisons as in `arbcheck`: `mag_cmp`, `arb_lt`, `arb_is_positive` /
   `arb_is_negative` on $`a \pm (e + \sum_k \lvert G_k\rvert H_k)`$, with $`\lvert G_k\rvert`$ from `arb_get_mag`. An uncertain comparison
   counts as failed.

Double precision is used only to choose the order of the attempts and, in the L test, the combination of
models. It never accepts a claim.

## Self-test

`./arbcheck2 --selftest d` runs the `arbcheck` self-test and adds the following:

- **Identity.** At up to 6 points of each box (the centre, 4 corners and 1 random point), for every $`i \ge 2`$,
  it evaluates $`T_i/u_i^{n-1}`$ and the w-form integral from their definitions in acb at 128 bits. It requires
  the two balls to overlap and their relative difference to be $`\le 10^{-20}`$. It also compares them in double
  precision.
- **Centre value.** The w-form constant term $`c_0`$ must equal $`S_i(m)`$ (direct evaluation), up to $`10^{-9}`$
  relative.
- **Bounds.** At 24 points per box (the centre, 4 corners and 19 random points), $`\lvert S_i\rvert`$ and $`\log\lvert S_i\rvert`$ come
  from direct evaluation of the defining integral. They must lie within the w-form bounds
  $`[S_{\text{lo}}, S_{\text{hi}}]`$ and within $`a + G\cdot(x-m) \pm e`$. The tolerances are those of `arbcheck`'s self-test: $`10^{-12}`$
  relative for magnitudes, and $`10^{-9}`$ relative plus $`10^{-12}`$ for logs.
- **Boxes.** 4000 random boxes: 2000 as in `arbcheck`, 1000 with $`\lvert m_i\rvert = r\rho_i`$ and $`r \in (1, 3)`$, and
  1000 with $`r \in (1, 1.062)`$. The last two groups put the boxes near the edge of the w-form domain.

Results (FLINT 3.6.0, 64-bit precision):

| d | identity points | max rel. diff (arb) | w-form models | not applicable | affine | point checks | max $`\lvert\log\lvert S\rvert-\text{affine}\rvert/e`$ | problems |
|---|---|---|---|---|---|---|---|---|
| 4 | 48000 | 4.3e-36 | 7939 | 61 | 5653 | 190536 | 0.970 | 0 |
| 5 | 72000 | 2.4e-35 | 11891 | 109 | 9017 | 285384 | 0.875 | 0 |
| 6 | 96000 | 2.8e-35 | 15853 | 147 | 12139 | 380472 | 0.887 | 0 |
| 7 | 120000 | 2.0e-34 | 19815 | 185 | 15313 | 475560 | 0.841 | 0 |

The T-form part of the self-test is unchanged and also reports 0 problems.

**Mutation check of the self-test** (scratch copies, not kept). Five deliberate errors were injected into
the w-form code, one at a time:

- dropping the $`q^2/(2(1-q))`$ term;
- dropping the $`\eta`$ term;
- using $`\rho_i`$ instead of $`\rho_w`$;
- dropping the minus sign of $`-\gamma_w/m_i^2`$;
- conjugating $`1/m_i`$.

The self-test reported problems for each of them at d = 5. The $`\eta`$ error, which is small, was detected at
d = 5 but not at d = 6.

**Adversarial check** (scratch files, not kept). Boxes containing an equality point were labelled F, E and L.
All three claims are false on such boxes. The check used d = 5 and 6, half-widths $`2^{-2}`$ to $`2^{-20}`$, 40
boxes per equality point, and `--maxref 0`. All 7200 claims were rejected.

## Runs

Settings: `-t 2` at `nice -n 19`, on a machine also running other jobs, with arbcheck2.c sha256
`01afc9afcbac359db9b57bf69bf33c209cc2a54f2c2714ffe4940f7fc487b3ad`. The logs are in `runs/d5_arbcheck2/`
(`regress_v2.log`, `v3.log`, run with `--quiet`).

**Regression, d = 5 tree of smale_bb_v2** (`arbcheck` reports PASS with 0 refined):

```
SUMMARY d=5 records=3072 nodes=2490060 leaves=1246566 F=757590 E=143562 L=100008 O=37156 S=205334 X=2916 U=0 direct=1246566 refined=0 subboxes=0 fail=0 failed_records=0 missing=0 dup=0 unknown=0 badheader=0 tail_incomplete=0 wform_direct=0 (F=0 E=0 L=0) wform_subboxes=0 maxdepth=47 wall=18.1s cpu=10.6s status=PASS
```

**New d = 5 tree of smale_bb_v3:**

```
SUMMARY d=5 records=3072 nodes=489000 leaves=246036 F=117338 E=39459 L=42686 O=3927 S=40953 X=1673 U=0 direct=246036 refined=0 subboxes=0 fail=0 failed_records=0 missing=0 dup=0 unknown=0 badheader=0 tail_incomplete=0 wform_direct=119738 (F=55375 E=22374 L=41989) wform_subboxes=0 maxdepth=46 wall=8.6s cpu=5.6s status=PASS
```

- Every leaf is proved directly, with no refinement.
- 119,738 of the 246,036 leaves (48.7 %) need the w-form: 55,375 F, 22,374 E and 41,989 L.
- The same result is obtained with `--maxref 0` (strict) and with `--prec 128`.

For comparison, on the same tree:

- `arbcheck`, and `arbcheck2 --no-wform`, both give `direct=126298 refined=119691 subboxes=1439397 fail=47 failed_records=44 status=FAIL`
  (wall 43.6 s and 58.4 s).
- 126,298 + 119,738 = 246,036, so every leaf the T-form could not prove directly is proved directly by the
  w-form.

## Build and usage

```
make arbcheck2             # same flags as arbcheck, with flint_compat.h force-included, and -lm
make selftest2             # ./arbcheck2 --selftest 4 .. 7
./arbcheck2 -t 2 5 runs/d5_v3/d5.tasks /path/to/d5_v3.tree
```

The options are those of `arbcheck` (`README.md`), plus `--no-wform` (T-form only). The first line on stdout
is a version line:

```
VERSION arbcheck2 2.0 (SPEC.md T-form + SPEC_W.md w-form), built with FLINT 3.6.0, w-form on
```

The task lines gain `wdirect=` and `wsub=` before `parse=`. The `SUMMARY` line gains
`wform_direct=N (F= E= L=) wform_subboxes=` before `maxdepth=`. A certificate run checks all records and needs
`status=PASS` with `missing=0`, as for `arbcheck`.

## Review
See `AI_DISCLOSURE.md`.
