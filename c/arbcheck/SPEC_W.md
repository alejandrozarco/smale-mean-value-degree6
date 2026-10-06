# Addendum to SPEC.md: a second model of $`S_i`$ in the variable $`w = 1/u_i`$

This addendum extends `SPEC.md` for a second checker version, `arbcheck2`. Everything in `SPEC.md` still applies:
the functions, the claims, the file formats, refinement and coverage. Only the suggested method (§3) gains one
ingredient. As before, the checker must be written from the specifications and the mathematics only.

## Why
Trees produced by a newer subdivision program have larger leaves. On many of them the centred form of §3 for
$`T_i = u_i^{\,n-1} S_i`$ is too weak. The relative variation of the factor $`u_i^{\,n-1}`$ is about $`(n-1)\rho_i/\lvert m_i\rvert`$, and it
pushes $`q`$ above $`1/2`$, even where $`\log\lvert S_i\rvert`$ itself varies little.

## The w-form ($`i = 2..n`$)
With $`w = 1/u_i`$:

```math
S_i(u) = \int_0^1 (1-t)\,(1 - t w) \prod_{j=2,\ j\ne i}^{n} (1 - t\,u_j\,w)\,dt .
```

This is $`T_i/u_i^{\,n-1}`$ written out. The checker should verify the identity, e.g. symbolically or in a self-test. It is
a polynomial of degree $`\le n-1`$ in $`w`$ and of degree $`\le 1`$ in each other $`u_j`$.

Let $`m`$ be the box centre, $`\delta_v = u_v - m_v`$ and $`\rho_v \ge \sqrt{H_{2v}^2 + H_{2v+1}^2}`$, as in §3.

1. **Domain.** Use the w-form for $`S_i`$ only if $`\lvert m_i\rvert \gt  \rho_i`$, certified. Then on the box
   $`\lvert u_i\rvert \ge \lvert m_i\rvert - \rho_i \gt  0`$, and $`w`$ lies in the closed disc
   $`\lvert w - 1/m_i\rvert \le \rho_w := \rho_i/(\lvert m_i\rvert(\lvert m_i\rvert - \rho_i))`$.
2. **Expansion.** Expand $`S_i`$ exactly at $`(w, u_j) = (1/m_i, m_j)`$ in the variables $`\delta_w = w - 1/m_i`$ and $`\delta_j`$. Use acb
   balls for the coefficients, with $`1/m_i`$ enclosed in a ball. Integrate the powers of $`t`$ with the exact moments
   $`1/((k+1)(k+2))`$.
3. **Disc bounds.** As in §3, with $`\rho_w`$ for the variable $`\delta_w`$. This gives bounds for $`\lvert S_i\rvert`$, and so for
   $`\log\lvert S_i\rvert`$, directly. No $`(n-1)\log\lvert u_i\rvert`$ term is needed.
4. **Affine model.** If $`q := (L + R)/\lvert c_0\rvert \lt  1/2`$, then
   $`\log\lvert S_i\rvert = \log\lvert c_0\rvert + \mathrm{Re}\left(\gamma_w \delta_w + \sum_j \gamma_j \delta_j\right) \pm \left(R/\lvert c_0\rvert + q^2/(2(1-q))\right)`$.
   To make it affine in the box coordinates, use

   ```math
   \delta_w = \frac{1}{m_i+\delta_i} - \frac{1}{m_i} = -\frac{\delta_i}{m_i^2} + \eta, \qquad
   \lvert\eta\rvert = \frac{\lvert\delta_i\rvert^2}{\lvert m_i\rvert^2\,\lvert m_i+\delta_i\rvert} \le \frac{\rho_i^2}{\lvert m_i\rvert^2(\lvert m_i\rvert-\rho_i)} .
   ```

   So $`\mathrm{Re}(\gamma_w\delta_w) = \mathrm{Re}(-\gamma_w\delta_i/m_i^2) \pm \lvert\gamma_w\rvert\,\rho_i^2/(\lvert m_i\rvert^2(\lvert m_i\rvert - \rho_i))`$. This gives a gradient vector and an error
   $`e_i`$ as in §3.4.
5. **Use.** For each $`i \ge 2`$ the checker may use the T-form, the w-form or both. Any bound from either form is valid
   on its own, and it may take the better one. The tests F/E/L are those of §3.5. Degenerate points (§3.6): the
   w-form needs $`\lvert u_i\rvert`$ bounded away from 0, so it is not used on boxes that meet $`u_i = 0`$.

The checker remains free to use any other sound method.
