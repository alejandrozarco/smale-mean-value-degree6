# Independent local certificate: specification

The certificate in this directory was written from this specification and the mathematics only, independently of
`local/local_cert.py` and its description.

## 1. Setting

$`d = 6`$, $`n = 5`$, $`c = 5/6`$, $`\omega = e^{2\pi i/5}`$.

For critical points $`b = (b_1, \dots, b_5) \in (\mathbb{C}^*)^5`$,

```math
S_i(b) = \int_0^1 \prod_{j=1}^{5} \Bigl(1 - t\,\frac{b_i}{b_j}\Bigr)\,dt, \qquad i = 1, \dots, 5 .
```

The factor $`j = i`$ is $`(1 - t)`$. $`S_i`$ is the normalised critical value $`P(b_i)/b_i`$ of the polynomial with
$`P(0) = 0`$, $`P'(0) = 1`$ and critical points $`b`$.

**Chart.** $`b_1 = 1`$ and $`b_j = \omega^{j-1} e^{\varepsilon_j}`$ for $`j = 2, \dots, 5`$. Here $`\varepsilon \in \mathbb{C}^4`$ and
$`x = (\mathrm{Re}\,\varepsilon_2, \mathrm{Im}\,\varepsilon_2, \dots, \mathrm{Re}\,\varepsilon_5, \mathrm{Im}\,\varepsilon_5) \in \mathbb{R}^8`$.
At $`x = 0`$ the critical points are the 5th roots of unity, and $`S_i = 5/6`$ for all $`i`$.

**Equal-modulus set.** $`E = \{x : \lvert S_1\rvert = \lvert S_2\rvert = \dots = \lvert S_5\rvert\}`$. On $`E`$,
$`\min_i \lvert S_i\rvert = \lvert S_1\rvert`$ and $`\log \min_i\lvert S_i\rvert = L(x) := \tfrac15 \sum_i \log\lvert S_i(x)\rvert`$.

## 2. Claim to certify (rigorously, in ball arithmetic)

**Claim.** There is an explicit $`T \ge 1/19`$ such that for every $`x \in E`$ with $`\lVert x\rVert_2 \le T`$ we have
$`L(x) \le \log c`$, i.e. $`\min_i \lvert S_i\rvert \le 5/6`$.

A quantitative form is preferred: $`L(x) \le \log c - \kappa \lVert x\rVert_2^2`$ for some explicit $`\kappa \gt  0`$.

Every constant and inequality in the final argument must be certified: exact arithmetic, or arb/acb ball arithmetic
where every comparison is certain. Floating point may be used only to choose things (bases, directions, targets) that
are verified afterwards.

Near $`x = 0`$ the deficit $`\log c - L`$ tends to 0 quadratically, so a plain box subdivision cannot close the claim.
Some second-order argument at $`x = 0`$ is needed. The method is otherwise the implementer's choice. Some possible
ingredients, none mandatory:
- Taylor models of $`\log S_i`$ in $`\varepsilon`$ with rigorous remainders;
- the structure of $`E`$ near 0 (the linear parts of $`\log\lvert S_i\rvert`$, their kernel, and bounds on how far $`E`$
  deviates from that kernel);
- a Hessian or Lagrangian bound;
- a branch and bound over spheres or shells;
- exact symbolic facts at 0.

If you use exact facts at $`x = 0`$ (such as $`S_i(0) = 5/6`$ or vanishing first derivatives of $`L`$), prove them yourself.
Examples are exact arithmetic in $`\mathbb{Q}(\omega)`$ or an analytic argument.

## 3. Hand-over (verify as part of the task)

The global computation sets aside only boxes contained in a closed Euclidean ball
$`\sum_{j=2}^{5} \lvert u_j - p_j\rvert^2 \le r^2`$ with $`r = 1/20`$. The coordinates are $`u_j = 1/b_j`$ (with $`b_1 = 1`$),
and $`p`$ is one of the 24 points $`p = (\bar\omega^{a_2}, \dots, \bar\omega^{a_5})`$, where $`(a_2, \dots, a_5)`$ is a
permutation of $`(1, 2, 3, 4)`$.

1. Show that relabelling (a permutation of $`b_2, \dots, b_5`$, which permutes the $`S_i`$ and leaves $`E`$ and
   $`\min_i \lvert S_i\rvert`$ invariant) maps each of these balls onto the ball around $`p = (\bar\omega, \bar\omega^2, \bar\omega^3, \bar\omega^4)`$.
2. Show that this ball, written in the $`\varepsilon`$-chart above, lies in $`\{\lVert x\rVert_2 \le 1/19\}`$.
   Alternatively, certify the claim directly on each of the 24 balls.
3. State exactly which set your certificate covers, and why it contains the image of every set-aside box.

## 4. Deliverables (in `local2/`)

- The code (python-flint 0.6 arb/acb, or C with FLINT 3.6 at `/opt/homebrew`). Note: in python-flint 0.6,
  `arb ** k` returns NaN when the ball contains 0; use repeated multiplication. Pure-Python exact arithmetic is fine
  for exact facts.
- `README.md`: the statement proved, the method step by step, every certified constant, how to run it, and the
  run time.
- The log of a full run, e.g. `local2/run_T....log`.

The run must be reproducible with one command and should take at most about 30 minutes on one core.
