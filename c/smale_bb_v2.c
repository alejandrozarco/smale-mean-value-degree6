/* smale_bb_v2.c — rigorous E-restricted branch and bound for Smale's mean value conjecture, degree d <= 7.
 * Version 2 (2026-09-29).  c/smale_bb.c (v1, sha256 1d943ca1…) is kept unchanged for provenance; v2 answers the
 * two arithmetic referee reports audit_C/astra_report.md and audit_C/AUDIT_CLAUDE.md.  Changes vs v1:
 *   A1  cb_addto / cb_subto add TINY (underflow allowance), like cb_mul.
 *   A2  The affine model of log|S_i| is declared INVALID unless |c0|^2 (constant Taylor coefficient) and |m_i|^2
 *       are >= 2^-900 and every model quantity is finite.  v1 silently dropped the numerator gradient when
 *       |c0|^2 underflowed to 0 (Astra's 2^-200 counterexample) — impossible now.
 *   A3  Complex moduli: no libm hypot.  cmod_hi_fast (upper bound, unscaled, absolute slack 2^-530) and cmod_bounds
 *       (two-sided, exact power-of-two rescaling, so no underflow of squares); proofs below.  modrange (symmetry /
 *       outside discards) uses cmod_bounds: Astra's s = 2^-537 box is no longer symmetry-discarded.
 *   A4  log = CORE-MATH cr_log (correctly rounded in round-to-nearest, third_party/core-math-log/, upstream commit in
 *       UPSTREAM_COMMIT); no libm pow (integer powers by repeated multiplication, powu below); no libm sin/cos
 *       (equality points from closed forms with IEEE sqrt only, root_of_unity below).  The only floating-point
 *       operations used are IEEE-754 +,-,*,/,sqrt (correctly rounded), fabs/frexp/ldexp (exact) and cr_log.
 *   A5  Domain hypotheses are ENFORCED, not assumed: every box handed to the evaluator must satisfy
 *         H_k a power of two, 2^-40 <= H_k <= 1,  C_k = m_k H_k with m_k an integer,  |C_k| + H_k <= 1.
 *       (All B&B boxes of a run satisfy this by construction; it is checked for every top box and asserted on every
 *       bisection.)  The grid parameter s must be a power of two <= 2^20 (s = 3 has gaps: rejected).
 *       Out-of-domain boxes are REFUSED (box/eval/models modes); the env var SMALE_V2_DIAG_NO_DOMAIN=1 lifts the
 *       check in the diagnostic modes only (never in `run`), to exhibit the arithmetic fixes at extreme scales.
 *   A6  Exclusion: only excl_mode 2 (box inside the closed Euclidean u-ball of radius r_excl around an equality
 *       point, = local/LOCAL_CERT.md) is accepted by `run`, and it is the default; r_excl <= 0.05 enforced.
 *       Absolute margin 1e-13 (>= sqrt(8) * max coordinate error of the closed-form equality points, < 1e-15).
 *   A7  Crash-safe, configuration-bound checkpoints: see "CHECKPOINT FORMAT" below.  Independent completion
 *       checker: c/check_done.py.
 *   A8  Build: c/build_v2.sh (records command, compiler version, hashes in runs/<run>/BUILD.txt).
 *       -O2 -ffp-contract=off, no -ffast-math.  The binary contains fmadd instructions ONLY from CORE-MATH's
 *       explicit __builtin_fma calls (part of cr_log's proved algorithm); none from contraction.
 *
 * Problem (see JOURNAL.md, refs/crane_equal_modulus.md):
 *   n = d-1 critical points, b_1 = 1 (minimal modulus), u_j = 1/b_j, |u_j| <= 1 (j = 2..n),
 *   real coordinates x = (Re u_2, Im u_2, Re u_3, ...), D = 2(n-1).  Symmetry: |u_2| >= ... >= |u_n|,
 *   Im u_2 >= 0.   S_1 = int_0^1 (1-t) prod_{j>=2} (1 - t u_j) dt,
 *   S_i = T_i / u_i^{n-1},  T_i = int_0^1 (1-t)(u_i - t) prod_{j != 1,i} (u_i - t u_j) dt.
 *   Goal: every box is (F) F < c, (E) disjoint from E = {|S_1| = ... = |S_n|}, (L) mean_i log|S_i| < log c,
 *   or inside the closed Euclidean u-ball of radius r_excl around an equality point (excl_mode 2).
 *
 * ARITHMETIC.  Complex balls (midpoint + radius) in IEEE-754 binary64, round-to-nearest (never changed).
 *   Model: fl(a op b) = (a op b)(1+δ) + η, |δ| <= U = 2^-53, |η| <= 2^-1075 (η = 0 for + and -), δη = 0.
 *   L1 (cb_mul/cb_addto/cb_subto).  Radius = propagated input radii + U * (sum of |computed partial results|)
 *      [covers δ-terms, since |fl(x) - x| <= U |fl(x)| in round-to-nearest for normal results] then * INFL
 *      (INFL - 1 = 2^-40 >> the <= 12 relative roundings committed while summing the nonnegative radius terms)
 *      + TINY (1e-300 >> the <= 8 underflow errors η of 2^-1075 per operation).
 *   L2 (cmod_hi_fast).  For |x|,|y| <= 2^500:  sqrt(x^2+y^2) <= fl(fl(sqrt(fl(x*x + y*y))) * (1+8U)) + 2^-530.
 *      Proof: s = fl(fl(x^2)+fl(y^2)) >= (x^2+y^2)(1-U)^2 - 2^-1073 ... see JOURNAL §13; the 2^-530 term dominates
 *      sqrt(2^-1073) = 2^-536.5 and the relative factor dominates (1-U)^-2 after its own rounding.
 *   L3 (cmod_bounds).  With 2^e the binade of max(|x|,|y|), xs = x 2^-e, ys = y 2^-e (exact unless ys underflows,
 *      error <= 2^-1075 << U relative to s >= 1/4), s = xs^2 + ys^2 in [1/4, 2], r = fl(sqrt(fl(s))) satisfies
 *      r = |z| 2^-e (1+θ), |θ| <= 2.5U.  hi = ldexp(fl(r(1+8U)), e) + 2^-1074, lo = max(0, ldexp(fl(r(1-8U)), e)
 *      - 2^-1074) (ldexp exact unless the result is subnormal, then error <= 2^-1075).  lo <= |z| <= hi.
 *   L4 (powu).  x >= 0, a <= 8, no underflow (x >= 2^-60 in all calls):  x^a <= fl(fl(x...x) (1 + 2(a+1)U)) + TINY.
 *   L5 (cr_log).  |cr_log(x) - log x| <= 0.5 ulp(log x) <= U |log x|; every log enters a test with an explicit
 *      allowance REL (1 + |log|) , REL = 1e-13 (>= 400 x the combined error of the scalar formulas).
 *   L6 (scalar model formulas, affine models, comparisons): as audited in v1 (audit_C/AUDIT_CLAUDE.md §1), with
 *      REL = 1e-13 per scalar quantity and MARG = 1e-12 relative margins on every final comparison.
 *
 * CHECKPOINT FORMAT (outprefix.done).  Every line is  "<payload> <chk>\n", chk = first 16 hex digits of
 *   sha256(payload).  First line: "H smale_bb_v2 cfg=<cfg16> <cfgstring>", cfgstring =
 *   "smale_bb_v2;src=<sha256 of this source>;d=<d>;r_excl=<16 hex digits of the IEEE bits>;excl_mode=<m>;s=<s>;seed=<seed>",
 *   cfg16 = first 16 hex digits of sha256(cfgstring).  Records: "R id processed F E L outside sym excl unresolved
 *   maxdepth sec vol cfg=<cfg16>".  Each record is written with ONE write(2) on an O_APPEND descriptor after the
 *   task has completely finished, then fsync'd (F_FULLFSYNC on macOS).  On resume every complete line is validated
 *   (checksum, config hash, id range, uniqueness, tree count processed = 2*leaves - 1); any invalid complete line
 *   aborts the run (manual inspection needed); a torn trailing line (no '\n') is truncated away (ftruncate + fsync).
 *   The file is flock'ed: two processes cannot write the same checkpoint.  outprefix.tasks lists the task boxes
 *   (hex floats) and is verified byte-for-byte on resume.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/stat.h>

#ifndef SRC_SHA256_RAW
#error "build with c/build_v2.sh (defines SRC_SHA256_RAW = sha256 of this file)"
#endif
#define STR2_(x) #x
#define STR_(x) STR2_(x)
#define SRC_SHA256 STR_(SRC_SHA256_RAW)

extern double cr_log(double);   /* CORE-MATH, third_party/core-math-log/log.c */
#define LOG cr_log

#define NMAX 6            /* max number of critical points (d <= 7) */
#define VMAX (NMAX - 1)   /* max number of complex variables */
#define DMAX (2 * VMAX)
#define TDEG (NMAX + 1)

static const double U = 1.1102230246251565e-16;   /* 2^-53 */
static const double INFL = 1.0 + 9.094947017729282e-13; /* 1 + 2^-40 */
static const double TINY = 1e-300;
static const double REL = 1e-13;   /* error allowance for scalar model quantities */
static const double MARG = 1e-12;  /* comparison margin */
static const double HMIN = 0x1p-40; /* domain: smallest admissible half-width */

typedef struct { double re, im, rad; } cb;

static inline double nrm1(double re, double im) { return fabs(re) + fabs(im); }

static inline cb cb_mul(cb a, cb b) {
    double p1 = a.re * b.re, p2 = a.im * b.im, p3 = a.re * b.im, p4 = a.im * b.re;
    cb r; r.re = p1 - p2; r.im = p3 + p4;
    double na = nrm1(a.re, a.im), nb = nrm1(b.re, b.im);
    double rad = na * b.rad + a.rad * nb + a.rad * b.rad
               + U * (fabs(p1) + fabs(p2) + fabs(r.re) + fabs(p3) + fabs(p4) + fabs(r.im));
    r.rad = rad * INFL + TINY;
    return r;
}
static inline void cb_addto(cb *q, cb a) {   /* q += a */
    double re = q->re + a.re, im = q->im + a.im;
    double rad = q->rad + a.rad + U * (fabs(re) + fabs(im));
    q->re = re; q->im = im; q->rad = rad * INFL + TINY;      /* v2: + TINY */
}
static inline void cb_subto(cb *q, cb a) {   /* q -= a */
    double re = q->re - a.re, im = q->im - a.im;
    double rad = q->rad + a.rad + U * (fabs(re) + fabs(im));
    q->re = re; q->im = im; q->rad = rad * INFL + TINY;      /* v2: + TINY */
}

/* L2: upper bound for |x + iy|, valid for |x|, |y| <= 2^500 */
static inline double cmod_hi_fast(double x, double y) {
    return sqrt(x * x + y * y) * (1 + 8 * U) + 0x1p-530;
}
/* L3: lo <= |x + iy| <= hi, any finite x, y */
static void cmod_bounds(double x, double y, double *lo, double *hi) {
    x = fabs(x); y = fabs(y);
    double m = x > y ? x : y;
    if (m == 0) { *lo = 0; *hi = 0; return; }
    int e; (void)frexp(m, &e);
    double xs = ldexp(x, -e), ys = ldexp(y, -e);
    double r = sqrt(xs * xs + ys * ys);
    *hi = ldexp(r * (1 + 8 * U), e) + 0x1p-1074;
    double l = ldexp(r * (1 - 8 * U), e) - 0x1p-1074;
    *lo = l > 0 ? l : 0;
}
/* L4: upper bound for x^a, x >= 0, 0 <= a <= 8 */
static inline double powu(double x, int a) {
    double r = 1.0;
    for (int k = 0; k < a; k++) r *= x;
    return r * (1 + 2 * (a + 1) * U) + TINY;
}

/* ---------------- problem data ---------------- */
typedef struct {
    int d, n, nv, D;
    double logc_lo;        /* lower bound for log((d-1)/d) */
    cb mom[TDEG];          /* int_0^1 (1-t) t^k dt = 1/((k+1)(k+2)) as balls */
    int next;              /* number of equality points */
    double *ext;           /* next * nv * 2 (re, im) */
    double r_excl;
    int excl_mode;         /* 2: box inside closed Euclidean ball sum_j |u_j-p_j|^2 <= r^2 (matches local/LOCAL_CERT.md).
                              0/1 (v1 rules) are available ONLY in the diagnostic `box` mode. */
} prob_t;

/* per-box models */
typedef struct {
    double a[NMAX], G[NMAX][DMAX], e[NMAX], lo[NMAX], hi[NMAX];
    int valid[NMAX];
} model_t;

/* domain hypotheses (A5) */
static int in_domain(int D, const double *C, const double *H) {
    for (int k = 0; k < D; k++) {
        double h = H[k], c = C[k];
        if (!(h >= HMIN && h <= 1.0)) return 0;
        int e; if (frexp(h, &e) != 0.5) return 0;          /* power of two */
        double m = c / h;                                   /* exact: h is a power of two in the normal range */
        if (!(m == floor(m))) return 0;                     /* integer (also rejects NaN/inf) */
        if (!(fabs(c) + h <= 1.0)) return 0;               /* exact: both are multiples of h, |.| <= 2 */
    }
    return 1;
}
static int diag_no_domain(void) {
    const char *s = getenv("SMALE_V2_DIAG_NO_DOMAIN");
    return s && !strcmp(s, "1");
}

/* Taylor coefficients of T_i at centre m (point values, exact doubles), as balls.  (Unchanged from v1.)
 * Output layout: coef[a][mask]; for i >= 1: a = power of delta_{vi}, mask over "others" (bit k <-> others[k]);
 * for i = 0: a = 0, mask over all nv variables. */
static void taylor_T(const prob_t *P, int i, const double *mre, const double *mim,
                     cb coef[NMAX][1 << VMAX], int *others, int *nbits) {
    int nv = P->nv;
    static __thread cb A[TDEG][NMAX][1 << VMAX], B[TDEG][NMAX][1 << VMAX];
    int tmax = 0, amax = 0, nb = 0;
    A[0][0][0] = (cb){1.0, 0.0, 0.0};
    if (i == 0) {
        for (int v = 0; v < nv; v++) {
            others[v] = v;
            cb c1 = {-mre[v], -mim[v], 0.0};
            int nm = 1 << nb, bit = 1 << nb;
            for (int t = 0; t <= tmax + 1; t++) for (int mk = 0; mk < 2 * nm; mk++) B[t][0][mk] = (cb){0, 0, 0};
            for (int t = 0; t <= tmax; t++) for (int mk = 0; mk < nm; mk++) {
                cb p = A[t][0][mk];
                cb_addto(&B[t][0][mk], p);                       /* * 1 */
                cb_addto(&B[t + 1][0][mk], cb_mul(p, c1));       /* * (-t m_v) */
                cb_subto(&B[t + 1][0][mk | bit], p);             /* * (-t delta_v) */
            }
            tmax++; nb++;
            for (int t = 0; t <= tmax; t++) for (int mk = 0; mk < (1 << nb); mk++) A[t][0][mk] = B[t][0][mk];
        }
        *nbits = nb;
        for (int mk = 0; mk < (1 << nb); mk++) {
            cb s = {0, 0, 0};
            for (int t = 0; t <= tmax; t++) cb_addto(&s, cb_mul(A[t][0][mk], P->mom[t]));
            coef[0][mk] = s;
        }
        return;
    }
    int vi = i - 1, k = 0;
    for (int v = 0; v < nv; v++) if (v != vi) others[k++] = v;
    cb c0 = {mre[vi], mim[vi], 0.0};
    /* factor (m_i - t) + delta_i */
    {
        for (int t = 0; t <= 1; t++) for (int a = 0; a <= 1; a++) B[t][a][0] = (cb){0, 0, 0};
        cb p = A[0][0][0];
        cb_addto(&B[0][0][0], cb_mul(p, c0));
        cb_subto(&B[1][0][0], p);
        cb_addto(&B[0][1][0], p);
        tmax = 1; amax = 1;
        for (int t = 0; t <= 1; t++) for (int a = 0; a <= 1; a++) A[t][a][0] = B[t][a][0];
    }
    for (int q = 0; q < nv - 1; q++) {
        int v = others[q];
        cb c1 = {-mre[v], -mim[v], 0.0};
        int nm = 1 << nb, bit = 1 << nb;
        for (int t = 0; t <= tmax + 1; t++) for (int a = 0; a <= amax + 1; a++)
            for (int mk = 0; mk < 2 * nm; mk++) B[t][a][mk] = (cb){0, 0, 0};
        for (int t = 0; t <= tmax; t++) for (int a = 0; a <= amax; a++) for (int mk = 0; mk < nm; mk++) {
            cb p = A[t][a][mk];
            if (p.re == 0.0 && p.im == 0.0 && p.rad == 0.0) continue;
            cb_addto(&B[t][a][mk], cb_mul(p, c0));          /* * m_i */
            cb_addto(&B[t + 1][a][mk], cb_mul(p, c1));      /* * (-t m_v) */
            cb_addto(&B[t][a + 1][mk], p);                  /* * delta_i */
            cb_subto(&B[t + 1][a][mk | bit], p);            /* * (-t delta_v) */
        }
        tmax++; amax++; nb++;
        for (int t = 0; t <= tmax; t++) for (int a = 0; a <= amax; a++)
            for (int mk = 0; mk < (1 << nb); mk++) A[t][a][mk] = B[t][a][mk];
    }
    *nbits = nb;
    for (int a = 0; a <= amax; a++) for (int mk = 0; mk < (1 << nb); mk++) {
        cb s = {0, 0, 0};
        for (int t = 0; t <= tmax; t++) cb_addto(&s, cb_mul(A[t][a][mk], P->mom[t]));
        coef[a][mk] = s;
    }
}

static inline double up(double x) { return x * INFL + TINY; }          /* for x >= 0 computed with few roundings */
static inline int isfin(double x) { return x - x == 0; }

static void build_models(const prob_t *P, const double *C, const double *H, model_t *M) {
    int n = P->n, nv = P->nv;
    double mre[VMAX], mim[VMAX], rho[VMAX], absm[VMAX], absm_lo[VMAX], absm_hi[VMAX];
    for (int v = 0; v < nv; v++) {
        mre[v] = C[2 * v]; mim[v] = C[2 * v + 1];
        double rl;
        cmod_bounds(H[2 * v], H[2 * v + 1], &rl, &rho[v]);          /* rho >= |(h_x, h_y)| */
        cmod_bounds(mre[v], mim[v], &absm_lo[v], &absm_hi[v]);
        absm[v] = sqrt(mre[v] * mre[v] + mim[v] * mim[v]);          /* midpoint estimate only (model centre) */
    }
    cb coef[NMAX][1 << VMAX];
    int others[VMAX], nb;
    for (int i = 0; i < n; i++) {
        taylor_T(P, i, mre, mim, coef, others, &nb);
        int amax = (i == 0) ? 0 : n - 1;
        double A0lo, A0hi;
        cmod_bounds(coef[0][0].re, coef[0][0].im, &A0lo, &A0hi);
        double rT = coef[0][0].rad, linmaj = 0.0;
        double gre[VMAX] = {0}, gim[VMAX] = {0}, gerr = 0.0;
        double c0re = coef[0][0].re, c0im = coef[0][0].im, c0n2 = c0re * c0re + c0im * c0im;
        /* A2: the affine model needs gamma = c/c0 for EVERY linear coefficient; if |c0|^2 is not safely in the
           normal range the model is declared invalid (never: silently drop the gradient term). */
        int grad_ok = (c0n2 >= 0x1p-900) && isfin(c0n2);
        for (int a = 0; a <= amax; a++) {
            double ra = (i == 0) ? 1.0 : powu(rho[i - 1], a);       /* L4: no libm pow */
            for (int mk = 0; mk < (1 << nb); mk++) {
                cb c = coef[a][mk];
                int deg = a + __builtin_popcount(mk);
                if (deg == 0) continue;
                double rp = ra;
                for (int q = 0; q < nb; q++) if (mk & (1 << q)) rp *= rho[others[q]];
                double ac = cmod_hi_fast(c.re, c.im);
                if (deg == 1) {
                    linmaj += ac * rp;
                    rT += c.rad * rp;
                    int var = (a == 1) ? i - 1 : others[__builtin_ctz(mk)];
                    if (grad_ok) {
                        /* gamma = c / c0 = c * conj(c0) / |c0|^2 ; numerator underflow <= 2^-1074 abs, / c0n2 >= 2^-900
                           => <= 2^-174 absolute, covered by the 0x1p-160 term */
                        double gr = (c.re * c0re + c.im * c0im) / c0n2;
                        double gi = (c.im * c0re - c.re * c0im) / c0n2;
                        gre[var] += gr; gim[var] += gi;
                        gerr += (REL * (fabs(gr) + fabs(gi)) + 0x1p-160) * rho[var];
                    }
                } else {
                    rT += (ac + c.rad) * rp;
                }
            }
        }
        rT = up(up(rT));          /* <= 40 nonnegative terms, each with <= 8 roundings: INFL^2 - 1 = 2^-39 >> 400 U */
        linmaj = up(up(linmaj));
        double supT = up(A0hi + linmaj + rT);
        double infT = (A0lo - linmaj - rT);
        infT = infT - 4 * U * (A0lo + linmaj + rT) - TINY;
        double dlo_log = 0, dhi_log = 0, dmag = 0;
        int dlo_ok = 1;
        if (i > 0) {
            int v = i - 1;
            double plo = (absm_lo[v] - rho[v]) * (1 - 4 * U);
            double phi = (absm_hi[v] + rho[v]) * (1 + 4 * U);
            if (plo > 0) { dlo_log = (n - 1) * LOG(plo); } else { dlo_ok = 0; }
            dhi_log = (n - 1) * LOG(phi);
            dmag = (n - 1) * (fabs(LOG(phi)) + (plo > 0 ? fabs(LOG(plo)) : 0));
        }
        /* interval bounds for log|S_i| */
        double ls = LOG(supT);
        double hi = (dlo_ok && supT < INFINITY) ? ls - dlo_log + REL * (1 + fabs(ls) + dmag) : INFINITY;
        double lo;
        if (infT > 0) { double li = LOG(infT); lo = li - dhi_log - REL * (1 + fabs(li) + dmag); }
        else lo = -INFINITY;
        if (!(hi == hi)) hi = INFINITY;
        if (!(lo == lo)) lo = -INFINITY;
        M->hi[i] = hi; M->lo[i] = lo;
        /* affine model */
        int ok = 0;
        double ai = 0, ei = INFINITY, amag = 0;
        if (grad_ok && A0lo > 0) {
            double q = up((linmaj + rT) / A0lo);
            if (q < 0.5) {
                ok = 1;
                double la = LOG(sqrt(c0n2));
                ai = la; amag = fabs(la);
                ei = up(rT / A0lo) + up(q * q / (2 * (1 - q)) * (1 + 4 * U));
                ei += gerr;
                if (i > 0) {
                    int v = i - 1;
                    double m2 = mre[v] * mre[v] + mim[v] * mim[v];
                    double p = (absm_lo[v] > 0) ? up(rho[v] / (absm_lo[v] * (1 - 4 * U))) : INFINITY;
                    if (!(p < 0.5) || !(m2 >= 0x1p-900)) ok = 0;
                    else {
                        double lm = LOG(absm[v]);
                        ai -= (n - 1) * lm; amag += (n - 1) * fabs(lm);
                        /* - (n-1)/m = -(n-1) conj(m)/|m|^2 */
                        double gr = -(n - 1) * mre[v] / m2, gi = (n - 1) * mim[v] / m2;
                        gre[v] += gr; gim[v] += gi;
                        ei += (REL * (fabs(gr) + fabs(gi)) + 0x1p-160) * rho[v];
                        ei += up((n - 1) * p * p / (2 * (1 - p)) * (1 + 4 * U));
                    }
                }
                ei = up(ei * (1 + 1e-12) + REL * (1 + amag + fabs(ai)));
            }
        }
        /* A2: explicit finiteness invariant */
        if (ok) {
            if (!isfin(ai) || !isfin(ei)) ok = 0;
            for (int v = 0; v < nv && ok; v++) if (!isfin(gre[v]) || !isfin(gim[v])) ok = 0;
        }
        M->valid[i] = ok;
        M->a[i] = ok ? ai : 0.0;
        M->e[i] = ok ? ei : INFINITY;
        for (int v = 0; v < nv; v++) {
            M->G[i][2 * v] = ok ? gre[v] : 0.0;
            M->G[i][2 * v + 1] = ok ? -gim[v] : 0.0;
        }
    }
}

/* returns 0 = not pruned, 1 = F, 2 = E, 4 = L ; also returns split coordinate.  (Unchanged from v1.) */
static int eval_box(const prob_t *P, const double *C, const double *H, int *split) {
    int n = P->n, D = P->D;
    model_t M;
    build_models(P, C, H, &M);
    double logc = P->logc_lo;
    int res = 0;
    /* F test */
    for (int i = 0; i < n && !res; i++) {
        if (M.hi[i] + MARG * (1 + fabs(M.hi[i])) < logc) res = 1;
        if (M.valid[i]) {
            double s = M.a[i], mag = fabs(M.a[i]);
            for (int k = 0; k < D; k++) { double t = fabs(M.G[i][k]) * H[k]; s += t; mag += t; }
            s += M.e[i]; mag += M.e[i];
            if (s + MARG * (1 + mag) < logc) res = 1;
        }
    }
    /* E test */
    if (!res) {
        double mlo = -INFINITY, mhi = INFINITY;
        for (int i = 0; i < n; i++) { if (M.lo[i] > mlo) mlo = M.lo[i]; if (M.hi[i] < mhi) mhi = M.hi[i]; }
        if (mlo - MARG * (1 + fabs(mlo)) > mhi + MARG * (1 + fabs(mhi))) res = 2;
        for (int i = 0; i < n && !res; i++) for (int j = i + 1; j < n && !res; j++) {
            if (!(M.valid[i] && M.valid[j])) continue;
            double lhs = fabs(M.a[i] - M.a[j]);
            double rhs = M.e[i] + M.e[j];
            for (int k = 0; k < D; k++) rhs += fabs(M.G[i][k] - M.G[j][k]) * H[k];
            if (lhs - rhs > MARG * (1 + fabs(M.a[i]) + fabs(M.a[j]) + rhs)) res = 2;
        }
    }
    /* L test (uniform weights): sum a + sum_k |sum_i G_ik| h_k + sum e < n log c */
    if (!res) {
        int allv = 1;
        for (int i = 0; i < n; i++) allv &= M.valid[i];
        if (allv) {
            double s = 0, mag = 0;
            for (int i = 0; i < n; i++) { s += M.a[i] + M.e[i]; mag += fabs(M.a[i]) + M.e[i]; }
            for (int k = 0; k < D; k++) {
                double g = 0, ga = 0;
                for (int i = 0; i < n; i++) { g += M.G[i][k]; ga += fabs(M.G[i][k]); }
                s += fabs(g) * H[k]; mag += ga * H[k];
            }
            if (s + MARG * (1 + mag) < n * logc - MARG) res = 4;
        }
    }
    if (!res && split) {
        double best = -1; int bk = 0;
        for (int k = 0; k < D; k++) {
            double gm = 0;
            for (int i = 0; i < n; i++) if (fabs(M.G[i][k]) > gm) gm = fabs(M.G[i][k]);
            double sc = H[k] * (1 + gm);
            if (sc > best) { best = sc; bk = k; }
        }
        *split = bk;
    }
    return res;
}

/* ---------------- geometry: outside / symmetry / exclusion ---------------- */
/* A3: rigorous two-sided bounds mn <= min_{box} |u_v|, mx >= max_{box} |u_v| (no underflow of squares). */
static void modrange(const double *C, const double *H, int v, double *mn, double *mx) {
    double dx = fabs(C[2 * v]) - H[2 * v], dy = fabs(C[2 * v + 1]) - H[2 * v + 1];   /* exact in the domain */
    if (dx < 0) dx = 0;
    if (dy < 0) dy = 0;
    double l, h;
    cmod_bounds(dx, dy, &l, &h);
    *mn = l;
    double ex = fabs(C[2 * v]) + H[2 * v], ey = fabs(C[2 * v + 1]) + H[2 * v + 1];
    cmod_bounds(ex, ey, &l, &h);
    *mx = h;
    /* In the domain dx, dy, ex, ey are exact.  Outside it (diagnostic modes) a normal-range result of |C| -+ H has
       relative error <= U (subnormal results of +,- are exact); the factors below absorb it. */
    *mn *= (1 - 4 * U); *mx *= (1 + 4 * U);
}
/* 1 = outside (|u_v| > 1 on whole box), 2 = symmetry-discarded, 0 = keep.  Conservative margins. */
static int geom_discard(const prob_t *P, const double *C, const double *H) {
    int nv = P->nv;
    double mn[VMAX], mx[VMAX];
    for (int v = 0; v < nv; v++) modrange(C, H, v, &mn[v], &mx[v]);
    for (int v = 0; v < nv; v++) if (mn[v] > 1.0 + 1e-12) return 1;
    for (int v = 0; v + 1 < nv; v++) if (mx[v] * (1 + 1e-12) < mn[v + 1] * (1 - 1e-12)) return 2;
    if (C[1] + H[1] < 0) return 2;   /* exact in the domain (and conservative otherwise: fl is monotone) */
    return 0;
}
static int excluded(const prob_t *P, const double *C, const double *H) {
    int nv = P->nv;
    if (P->excl_mode != 2) {       /* v1 rules, diagnostic `box` mode only */
        double hmax = 0;
        for (int k = 0; k < P->D; k++) if (H[k] > hmax) hmax = H[k];
        for (int e = 0; e < P->next; e++) {
            const double *p = P->ext + (size_t)e * nv * 2;
            int in = 1; double dmax = 0;
            for (int v = 0; v < nv && in; v++) {
                double ax = fabs(C[2 * v] - p[2 * v]), ay = fabs(C[2 * v + 1] - p[2 * v + 1]);
                if (P->excl_mode == 1) {
                    double fx = ax + H[2 * v], fy = ay + H[2 * v + 1];
                    if (!(sqrt(fx * fx + fy * fy) * (1 + 1e-12) + 1e-13 <= P->r_excl)) in = 0;
                } else {
                    double nx = ax - H[2 * v], ny = ay - H[2 * v + 1];
                    if (nx < 0) nx = 0;
                    if (ny < 0) ny = 0;
                    double nd = sqrt(nx * nx + ny * ny);
                    if (nd > dmax) dmax = nd;
                }
            }
            if (P->excl_mode == 1 && in) return 1;
            if (P->excl_mode == 0 && dmax < P->r_excl && hmax < P->r_excl) return 1;
        }
        return 0;
    }
    /* mode 2: sup_{u in box} |u - p|_2 <= sqrt(sum_k (|C_k - p_k| + H_k)^2); the computed value has relative error
       <= (D+4)U << 1e-12, underflow of squares <= D 2^-1075 << 1e-13^2, and the computed p differs from the true
       equality point by <= 1e-15 per coordinate (closed forms, verified in c/regress_v2.py), i.e. <= 3e-15 in the
       Euclidean norm << 1e-13.  Hence the test below implies box ⊂ closed ball(p_true, r_excl). */
    for (int e = 0; e < P->next; e++) {
        const double *p = P->ext + (size_t)e * nv * 2;
        double far2 = 0;
        for (int v = 0; v < nv; v++) {
            double fx = fabs(C[2 * v] - p[2 * v]) + H[2 * v], fy = fabs(C[2 * v + 1] - p[2 * v + 1]) + H[2 * v + 1];
            far2 += fx * fx + fy * fy;
        }
        if (sqrt(far2) * (1 + 1e-12) + 1e-13 <= P->r_excl) return 1;
    }
    return 0;
}

/* ---------------- branch and bound on one top box ---------------- */
typedef struct {
    uint64_t processed, F, E, L, outside, sym, excl, unresolved;
    int maxdepth;
    double vol_done;
} stats_t;

typedef struct { double c[DMAX], h[DMAX]; int depth; } box_t;

typedef struct {
    FILE *samp; double samp_p; uint64_t rng; pthread_mutex_t *mtx;
    int unres_fd; FILE *unres_fp; const char *cfg16;
} aux_t;

static inline double urand(uint64_t *s) {
    *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17;
    return (*s >> 11) * (1.0 / 9007199254740992.0);
}

static void full_sync(int fd) {
#ifdef F_FULLFSYNC
    if (fcntl(fd, F_FULLFSYNC) == 0) return;
#endif
    if (fsync(fd) != 0) { perror("fsync"); exit(3); }
}
static void write_all(int fd, const char *buf, size_t len) {
    while (len) {
        ssize_t w = write(fd, buf, len);
        if (w < 0) { if (errno == EINTR) continue; perror("write"); exit(3); }
        buf += w; len -= (size_t)w;
    }
}

static void bb_box(const prob_t *P, const double *c0, const double *h0, stats_t *st, aux_t *ax) {
    int D = P->D;
    if (!in_domain(D, c0, h0)) { fprintf(stderr, "FATAL: top box outside the domain hypotheses (A5)\n"); exit(4); }
    size_t cap = 4096, sp = 0;
    box_t *stk = malloc(cap * sizeof(box_t));
    if (!stk) { perror("malloc"); exit(3); }
    memcpy(stk[0].c, c0, D * sizeof(double)); memcpy(stk[0].h, h0, D * sizeof(double)); stk[0].depth = 0;
    sp = 1;
    while (sp) {
        box_t b = stk[--sp];
        st->processed++;
        double vol = 1; for (int k = 0; k < D; k++) vol *= 2 * b.h[k];
        int g = geom_discard(P, b.c, b.h);
        if (g) { if (g == 1) st->outside++; else st->sym++; st->vol_done += vol; continue; }
        if (excluded(P, b.c, b.h)) { st->excl++; st->vol_done += vol; continue; }
        int sk = 0;
        int r = eval_box(P, b.c, b.h, &sk);
        if (ax && ax->samp && urand(&ax->rng) < ax->samp_p) {
            pthread_mutex_lock(ax->mtx);
            fprintf(ax->samp, "%d", r);
            for (int k = 0; k < D; k++) fprintf(ax->samp, " %.17g", b.c[k]);
            for (int k = 0; k < D; k++) fprintf(ax->samp, " %.17g", b.h[k]);
            fprintf(ax->samp, "\n");
            pthread_mutex_unlock(ax->mtx);
        }
        if (r) {
            if (r == 1) st->F++; else if (r == 2) st->E++; else st->L++;
            st->vol_done += vol; continue;
        }
        if (b.h[sk] < 1e-9) {   /* give up: record as unresolved (must not happen in a proof run) */
            st->unresolved++; st->vol_done += vol;
            if (ax && (ax->unres_fd >= 0 || ax->unres_fp)) {
                char line[1024]; int L = snprintf(line, sizeof line, "UNRESOLVED cfg=%s", ax->cfg16 ? ax->cfg16 : "-");
                for (int k = 0; k < D; k++) L += snprintf(line + L, sizeof line - L, " %a", b.c[k]);
                for (int k = 0; k < D; k++) L += snprintf(line + L, sizeof line - L, " %a", b.h[k]);
                L += snprintf(line + L, sizeof line - L, "\n");
                pthread_mutex_lock(ax->mtx);
                if (ax->unres_fd >= 0) { write_all(ax->unres_fd, line, (size_t)L); full_sync(ax->unres_fd); }
                else { fputs(line, ax->unres_fp); fflush(ax->unres_fp); }
                pthread_mutex_unlock(ax->mtx);
            }
            continue;
        }
        if (sp + 2 > cap) { cap *= 2; stk = realloc(stk, cap * sizeof(box_t)); if (!stk) { perror("realloc"); exit(3); } }
        box_t b1 = b, b2 = b;
        double hh = b.h[sk] * 0.5;           /* exact (dyadic) */
        if (!(hh >= HMIN)) { fprintf(stderr, "FATAL: bisection below 2^-40 (A5)\n"); exit(4); }
        b1.h[sk] = hh; b2.h[sk] = hh;
        b1.c[sk] = b.c[sk] - hh; b2.c[sk] = b.c[sk] + hh;
        b1.depth = b2.depth = b.depth + 1;
        if (b1.depth > st->maxdepth) st->maxdepth = b1.depth;
        stk[sp++] = b2; stk[sp++] = b1;
    }
    free(stk);
}

/* ---------------- setup ---------------- */
/* equality points: u_j = conj(omega^{a_j}), omega = e^{2 pi i/n}; closed forms using IEEE sqrt only (A4).
   Coordinate error <= 1e-15 (checked against mpmath in c/regress_v2.py). */
static void root_of_unity(int n, int k, double *re, double *im) {
    k %= n; if (k < 0) k += n;
    double c = 0, s = 0;
    if (n == 4) { static const double C4[4] = {1, 0, -1, 0}, S4[4] = {0, 1, 0, -1}; c = C4[k]; s = S4[k]; }
    else if (n == 3 || n == 6) {
        int kk = (n == 3) ? 2 * k : k;         /* angle = kk * pi/3 */
        double h3 = sqrt(3.0) / 2;
        static const double C6[6] = {1, 0.5, -0.5, -1, -0.5, 0.5};
        static const int S6[6] = {0, 1, 1, 0, -1, -1};
        c = C6[kk]; s = S6[kk] * h3;
    } else if (n == 5) {
        double r5 = sqrt(5.0);
        double c1 = (r5 - 1) / 4, c2 = -(r5 + 1) / 4, s1 = sqrt((5 + r5) / 8), s2 = sqrt((5 - r5) / 8);
        double Cs[5] = {1, c1, c2, c2, c1}, Ss[5] = {0, s1, s2, -s2, -s1};
        c = Cs[k]; s = Ss[k];
    } else { fprintf(stderr, "root_of_unity: n=%d unsupported\n", n); exit(1); }
    *re = c; *im = s;
}
static void perm_rec(int *a, int k, int m, double *out, int *cnt, int n) {
    if (k == m) {
        for (int j = 0; j < m; j++) {
            double re, im;
            root_of_unity(n, a[j], &re, &im);
            out[(size_t)(*cnt) * m * 2 + 2 * j] = re;
            out[(size_t)(*cnt) * m * 2 + 2 * j + 1] = -im;       /* conj */
        }
        (*cnt)++;
        return;
    }
    for (int i = k; i < m; i++) {
        int t = a[k]; a[k] = a[i]; a[i] = t;
        perm_rec(a, k + 1, m, out, cnt, n);
        t = a[k]; a[k] = a[i]; a[i] = t;
    }
}
static void setup(prob_t *P, int d, double r_excl, int excl_mode) {
    P->d = d; P->n = d - 1; P->nv = d - 2; P->D = 2 * (d - 2);
    P->logc_lo = LOG((double)(d - 1) / d) - 1e-14;
    for (int k = 0; k < TDEG; k++) {
        double mid = 1.0 / ((double)(k + 1) * (k + 2));
        P->mom[k] = (cb){mid, 0.0, 2 * U * mid};
    }
    int m = P->nv, fact = 1;
    for (int k = 2; k <= m; k++) fact *= k;
    P->ext = malloc(sizeof(double) * fact * m * 2);
    int a[VMAX], cnt = 0;
    for (int j = 0; j < m; j++) a[j] = j + 1;
    perm_rec(a, 0, m, P->ext, &cnt, P->n);
    P->next = cnt;
    P->r_excl = r_excl; P->excl_mode = excl_mode;
}

/* ---------------- SHA-256 (FIPS 180-4) ---------------- */
typedef struct { uint32_t h[8]; uint64_t len; uint8_t buf[64]; size_t bl; } sha_t;
static const uint32_t K256[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define ROR(x, r) (((x) >> (r)) | ((x) << (32 - (r))))
static void sha_block(sha_t *s, const uint8_t *p) {
    uint32_t w[64];
    for (int i = 0; i < 16; i++) w[i] = (uint32_t)p[4*i] << 24 | (uint32_t)p[4*i+1] << 16 | (uint32_t)p[4*i+2] << 8 | p[4*i+3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = ROR(w[i-15], 7) ^ ROR(w[i-15], 18) ^ (w[i-15] >> 3);
        uint32_t s1 = ROR(w[i-2], 17) ^ ROR(w[i-2], 19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    uint32_t a = s->h[0], b = s->h[1], c = s->h[2], d = s->h[3], e = s->h[4], f = s->h[5], g = s->h[6], h = s->h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t t1 = h + (ROR(e, 6) ^ ROR(e, 11) ^ ROR(e, 25)) + ((e & f) ^ (~e & g)) + K256[i] + w[i];
        uint32_t t2 = (ROR(a, 2) ^ ROR(a, 13) ^ ROR(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    s->h[0] += a; s->h[1] += b; s->h[2] += c; s->h[3] += d; s->h[4] += e; s->h[5] += f; s->h[6] += g; s->h[7] += h;
}
static void sha256_hex(const void *data, size_t len, char out[65]) {
    sha_t s = {{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19}, 0, {0}, 0};
    const uint8_t *p = data;
    uint64_t bits = (uint64_t)len * 8;
    while (len >= 64) { sha_block(&s, p); p += 64; len -= 64; }
    uint8_t buf[128] = {0};
    memcpy(buf, p, len); buf[len] = 0x80;
    size_t tot = (len + 9 <= 64) ? 64 : 128;
    for (int i = 0; i < 8; i++) buf[tot - 1 - i] = (uint8_t)(bits >> (8 * i));
    sha_block(&s, buf); if (tot == 128) sha_block(&s, buf + 64);
    for (int i = 0; i < 8; i++) snprintf(out + 8 * i, 9, "%08x", s.h[i]);
}
/* line = payload + " " + chk16 + "\n" */
static int make_line(char *dst, size_t cap, const char *payload) {
    char hx[65]; sha256_hex(payload, strlen(payload), hx);
    return snprintf(dst, cap, "%s %.16s\n", payload, hx);
}
static int check_line(const char *line, size_t len, char *payload_out, size_t cap) {
    /* line without the trailing '\n'; returns 1 if the checksum is valid */
    if (len < 18 || line[len - 17] != ' ' || len - 17 >= cap) return 0;
    memcpy(payload_out, line, len - 17); payload_out[len - 17] = 0;
    char hx[65]; sha256_hex(payload_out, len - 17, hx);
    return memcmp(hx, line + len - 16, 16) == 0;
}

/* ---------------- driver ---------------- */
typedef struct {
    const prob_t *P; int ntask; double *tc; double th; int *done;
    int next; pthread_mutex_t mtx; int done_fd; stats_t tot; double t0;
    FILE *samp; double samp_p; int unres_fd; uint64_t seed; const char *cfg16;
} work_t;

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static void *worker(void *arg) {
    work_t *W = arg;
    const prob_t *P = W->P;
    aux_t ax = {W->samp, W->samp_p, 0, &W->mtx, W->unres_fd, NULL, W->cfg16};
    pthread_mutex_lock(&W->mtx);
    ax.rng = W->seed = W->seed * 6364136223846793005ULL + 1442695040888963407ULL;
    if (!ax.rng) ax.rng = 1;
    pthread_mutex_unlock(&W->mtx);
    double h[DMAX];
    for (int k = 0; k < P->D; k++) h[k] = W->th;
    for (;;) {
        pthread_mutex_lock(&W->mtx);
        while (W->next < W->ntask && W->done[W->next]) W->next++;
        int t = W->next < W->ntask ? W->next++ : -1;
        pthread_mutex_unlock(&W->mtx);
        if (t < 0) break;
        stats_t st; memset(&st, 0, sizeof st);
        double t1 = now();
        bb_box(P, W->tc + (size_t)t * P->D, h, &st, &ax);
        double el = now() - t1;
        char payload[512], line[600];
        snprintf(payload, sizeof payload, "R %d %llu %llu %llu %llu %llu %llu %llu %llu %d %.3f %.17g cfg=%s", t,
                 (unsigned long long)st.processed, (unsigned long long)st.F, (unsigned long long)st.E,
                 (unsigned long long)st.L, (unsigned long long)st.outside, (unsigned long long)st.sym,
                 (unsigned long long)st.excl, (unsigned long long)st.unresolved, st.maxdepth, el, st.vol_done, W->cfg16);
        int L = make_line(line, sizeof line, payload);
        pthread_mutex_lock(&W->mtx);
        write_all(W->done_fd, line, (size_t)L);          /* one write(2) on O_APPEND */
        W->tot.processed += st.processed; W->tot.F += st.F; W->tot.E += st.E; W->tot.L += st.L;
        W->tot.outside += st.outside; W->tot.sym += st.sym; W->tot.excl += st.excl;
        W->tot.unresolved += st.unresolved; W->tot.vol_done += st.vol_done;
        if (st.maxdepth > W->tot.maxdepth) W->tot.maxdepth = st.maxdepth;
        pthread_mutex_unlock(&W->mtx);
        full_sync(W->done_fd);
    }
    return NULL;
}

static void usage(void) {
    fprintf(stderr,
        "usage:\n"
        "  smale_bb_v2 run  d r_excl s nthreads outprefix [excl_mode=2] [shuffle_seed=1]\n"
        "      certification run: grid of s^D top boxes (s a power of two), excl_mode must be 2, r_excl <= 0.05.\n"
        "      Resumable: validated records in outprefix.done are skipped (see CHECKPOINT FORMAT in the source).\n"
        "  smale_bb_v2 box  d r_excl excl_mode c1..cD h1..hD      (single box, diagnostic; prints stats)\n"
        "  smale_bb_v2 eval d        (stdin: lines c1..cD h1..hD; stdout: 'flag split' or 'REFUSED')\n"
        "  smale_bb_v2 models d      (stdin boxes; stdout per i: valid a e lo hi G_1..G_D, or 'REFUSED')\n"
        "  smale_bb_v2 geom d        (stdin boxes; stdout: geometric discard 0 keep / 1 outside / 2 symmetry)\n"
        "  smale_bb_v2 config d r_excl s [excl_mode=2] [seed=1]   (print the config string and hash)\n"
        "  smale_bb_v2 roots d       (print the equality points, hex)\n"
        "  smale_bb_v2 sha256 file   (self-test of the built-in SHA-256)\n"
        "  env SMALE_V2_DIAG_NO_DOMAIN=1 lifts the domain check in box/eval/models/geom (never in run).\n");
    exit(1);
}

static void cfg_string(char *out, size_t cap, int d, double r, int em, int s, unsigned long long seed) {
    uint64_t bits; memcpy(&bits, &r, 8);
    snprintf(out, cap, "smale_bb_v2;src=%s;d=%d;r_excl=%016llx;excl_mode=%d;s=%d;seed=%llu",
             SRC_SHA256, d, (unsigned long long)bits, em, s, seed);
}

static int read_box(int D, double *c, double *h) {
    int ok = 1;
    for (int k = 0; k < D && ok; k++) ok = scanf("%lf", &c[k]) == 1;
    for (int k = 0; k < D && ok; k++) ok = scanf("%lf", &h[k]) == 1;
    return ok;
}

int main(int argc, char **argv) {
    if (argc < 3) usage();
    if (!strcmp(argv[1], "sha256")) {
        FILE *f = fopen(argv[2], "rb"); if (!f) { perror(argv[2]); return 1; }
        fseek(f, 0, SEEK_END); long L = ftell(f); fseek(f, 0, SEEK_SET);
        char *b = malloc(L + 1); if (fread(b, 1, L, f) != (size_t)L) return 1; fclose(f);
        char hx[65]; sha256_hex(b, L, hx); printf("%s  %s\n", hx, argv[2]); return 0;
    }
    int d = atoi(argv[2]);
    if (d < 4 || d > NMAX + 1) { fprintf(stderr, "d must be in 4..%d\n", NMAX + 1); return 1; }
    prob_t P;
    int nodom = diag_no_domain();
    if (!strcmp(argv[1], "roots")) {
        setup(&P, d, 0.0, 2);
        for (int e = 0; e < P.next; e++) {
            for (int k = 0; k < P.D; k++) printf("%a ", P.ext[(size_t)e * P.D + k]);
            printf("\n");
        }
        return 0;
    }
    if (!strcmp(argv[1], "config")) {
        if (argc < 5) usage();
        char cs[512], hx[65];
        cfg_string(cs, sizeof cs, d, atof(argv[3]), argc > 5 ? atoi(argv[5]) : 2, atoi(argv[4]),
                   argc > 6 ? strtoull(argv[6], 0, 10) : 1ULL);
        sha256_hex(cs, strlen(cs), hx);
        printf("%s\n%.16s\n", cs, hx);
        return 0;
    }
    if (!strcmp(argv[1], "eval") || !strcmp(argv[1], "models") || !strcmp(argv[1], "geom")) {
        setup(&P, d, 0.0, 2);
        double c[DMAX], h[DMAX];
        model_t M;
        while (read_box(P.D, c, h)) {
            if (!nodom && !in_domain(P.D, c, h)) { printf("REFUSED\n"); continue; }
            if (!strcmp(argv[1], "eval")) {
                int sk = -1;
                int r = eval_box(&P, c, h, &sk);
                printf("%d %d\n", r, r ? -1 : sk);
            } else if (!strcmp(argv[1], "geom")) {
                printf("%d\n", geom_discard(&P, c, h));
            } else {
                build_models(&P, c, h, &M);
                for (int i = 0; i < P.n; i++) {
                    printf("%d %.17g %.17g %.17g %.17g", M.valid[i], M.a[i], M.e[i], M.lo[i], M.hi[i]);
                    for (int k = 0; k < P.D; k++) printf(" %.17g", M.G[i][k]);
                    printf("\n");
                }
            }
        }
        return 0;
    }
    if (!strcmp(argv[1], "box")) {
        if (argc < 5) usage();
        setup(&P, d, atof(argv[3]), atoi(argv[4]));
        if (argc < 5 + 2 * P.D) usage();
        if (P.excl_mode != 2) fprintf(stderr, "WARNING: excl_mode %d is NOT a certification mode (diagnostic only)\n", P.excl_mode);
        double c[DMAX], h[DMAX];
        for (int k = 0; k < P.D; k++) { c[k] = atof(argv[5 + k]); h[k] = atof(argv[5 + P.D + k]); }
        if (!in_domain(P.D, c, h)) {
            if (!nodom) { printf("{\"refused\": \"box outside the domain hypotheses (A5)\"}\n"); return 4; }
        }
        stats_t st; memset(&st, 0, sizeof st);
        pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
        aux_t ax = {NULL, 0, 1, &m, -1, stderr, NULL};
        if (getenv("SMALE_SAMPLE_FILE") && getenv("SMALE_SAMPLE_P")) {
            ax.samp = fopen(getenv("SMALE_SAMPLE_FILE"), "a");
            ax.samp_p = atof(getenv("SMALE_SAMPLE_P"));
            ax.rng = getenv("SMALE_SEED") ? strtoull(getenv("SMALE_SEED"), 0, 10) : 88172645463325252ULL;
        }
        double t0 = now();
        if (nodom && !in_domain(P.D, c, h)) {      /* single evaluation only (bisection needs the domain) */
            st.processed = 1;
            int g = geom_discard(&P, c, h);
            if (g == 1) st.outside++; else if (g == 2) st.sym++;
            else if (excluded(&P, c, h)) st.excl++;
            else { int r = eval_box(&P, c, h, NULL); if (r == 1) st.F++; else if (r == 2) st.E++; else if (r == 4) st.L++; }
        } else bb_box(&P, c, h, &st, &ax);
        double el = now() - t0;
        if (ax.samp) fclose(ax.samp);
        double vol = 1; for (int k = 0; k < P.D; k++) vol *= 2 * h[k];
        printf("{\"d\": %d, \"r_excl\": %g, \"excl_mode\": %d, \"processed\": %llu, \"F\": %llu, \"E\": %llu, \"L\": %llu, "
               "\"outside\": %llu, \"symmetry\": %llu, \"excluded\": %llu, \"unresolved\": %llu, \"maxdepth\": %d, "
               "\"vol_frac\": %.15g, \"time\": %.3f, \"boxes_per_s\": %.0f}\n",
               d, P.r_excl, P.excl_mode, (unsigned long long)st.processed, (unsigned long long)st.F,
               (unsigned long long)st.E, (unsigned long long)st.L, (unsigned long long)st.outside,
               (unsigned long long)st.sym, (unsigned long long)st.excl, (unsigned long long)st.unresolved,
               st.maxdepth, st.vol_done / vol, el, st.processed / (el > 0 ? el : 1e-9));
        return 0;
    }
    if (!strcmp(argv[1], "run")) {
        if (argc < 7) usage();
        double r = atof(argv[3]); int s = atoi(argv[4]), nth = atoi(argv[5]);
        const char *pre = argv[6];
        int em = argc > 7 ? atoi(argv[7]) : 2;
        unsigned long long seed = argc > 8 ? strtoull(argv[8], 0, 10) : 1;
        if (em != 2) { fprintf(stderr, "REFUSED: excl_mode %d is not a certification mode; only excl_mode 2 matches local/LOCAL_CERT.md\n", em); return 5; }
        if (!(r > 0 && r <= 0.05)) { fprintf(stderr, "REFUSED: r_excl must be in (0, 0.05] (local certificate radius 0.0527 >= r/(1-r))\n"); return 5; }
        if (s < 1 || s > (1 << 20) || (s & (s - 1))) { fprintf(stderr, "REFUSED: s must be a power of two <= 2^20 (s=%d: grid would not be exact)\n", s); return 5; }
        if (nth < 1 || nth > 64) { fprintf(stderr, "nthreads must be in 1..64\n"); return 1; }
        if (getenv("SMALE_V2_DIAG_NO_DOMAIN")) fprintf(stderr, "note: SMALE_V2_DIAG_NO_DOMAIN is ignored in run mode\n");
        setup(&P, d, r, em);
        int D = P.D;
        char cfgs[512], cfghex[65], cfg16[17];
        cfg_string(cfgs, sizeof cfgs, d, r, em, s, seed);
        sha256_hex(cfgs, strlen(cfgs), cfghex);
        memcpy(cfg16, cfghex, 16); cfg16[16] = 0;
        /* top boxes: grid of s per coordinate on [-1,1]^D (exact dyadics), top-level discards, deterministic shuffle */
        double th = 1.0 / s;
        long total = 1; for (int k = 0; k < D; k++) total *= s;
        double *tc = malloc(sizeof(double) * total * D);
        int ntask = 0; long ntriv = 0;
        double c[DMAX], h[DMAX];
        for (int k = 0; k < D; k++) h[k] = th;
        for (long idx = 0; idx < total; idx++) {
            long q = idx;
            for (int k = D - 1; k >= 0; k--) { c[k] = -1 + th * (2 * (q % s) + 1); q /= s; }
            if (!in_domain(D, c, h)) { fprintf(stderr, "FATAL: grid box outside domain\n"); return 4; }
            if (geom_discard(&P, c, h)) { ntriv++; continue; }
            memcpy(tc + (size_t)ntask * D, c, D * sizeof(double));
            ntask++;
        }
        uint64_t rs = seed * 2654435761ULL + 12345;
        for (int i = ntask - 1; i > 0; i--) {
            rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
            int j = rs % (i + 1);
            for (int k = 0; k < D; k++) { double t = tc[(size_t)i * D + k]; tc[(size_t)i * D + k] = tc[(size_t)j * D + k]; tc[(size_t)j * D + k] = t; }
        }
        char fn[1024];
        /* tasks file: written once, verified byte-for-byte on resume */
        {
            size_t cap = (size_t)ntask * (D * 26 + 16) + 1024, L = 0;
            char *buf = malloc(cap);
            L += snprintf(buf + L, cap - L, "# smale_bb_v2 tasks cfg=%s %s ntask=%d half_width=%a\n", cfg16, cfgs, ntask, th);
            for (int t = 0; t < ntask; t++) {
                L += snprintf(buf + L, cap - L, "%d", t);
                for (int k = 0; k < D; k++) L += snprintf(buf + L, cap - L, " %a", tc[(size_t)t * D + k]);
                L += snprintf(buf + L, cap - L, "\n");
            }
            snprintf(fn, sizeof fn, "%s.tasks", pre);
            FILE *f = fopen(fn, "rb");
            if (f) {
                char *old = malloc(L + 2); size_t got = fread(old, 1, L + 1, f); fclose(f);
                if (got != L || memcmp(old, buf, L)) { fprintf(stderr, "REFUSED: %s exists but differs from this configuration\n", fn); return 6; }
                free(old);
            } else {
                int fd = open(fn, O_WRONLY | O_CREAT | O_EXCL, 0644);
                if (fd < 0) { perror(fn); return 3; }
                write_all(fd, buf, L); full_sync(fd); close(fd);
            }
            free(buf);
        }
        /* checkpoint: open, lock, validate, truncate torn tail */
        int *done = calloc(ntask, sizeof(int));
        snprintf(fn, sizeof fn, "%s.done", pre);
        int fd = open(fn, O_RDWR | O_CREAT | O_APPEND, 0644);
        if (fd < 0) { perror(fn); return 3; }
        if (flock(fd, LOCK_EX | LOCK_NB) != 0) { fprintf(stderr, "REFUSED: %s is locked by another process\n", fn); return 6; }
        struct stat sb; fstat(fd, &sb);
        size_t fl = (size_t)sb.st_size;
        char *fb = malloc(fl + 1);
        if (fl && pread(fd, fb, fl, 0) != (ssize_t)fl) { perror("pread"); return 3; }
        size_t good = 0; int ndone = 0, have_header = 0, lineno = 0;
        work_t W; memset(&W, 0, sizeof W);
        char payload[1024];
        while (good < fl) {
            char *nl = memchr(fb + good, '\n', fl - good);
            if (!nl) break;                                  /* torn tail */
            size_t len = (size_t)(nl - (fb + good));
            lineno++;
            if (!check_line(fb + good, len, payload, sizeof payload)) {
                fprintf(stderr, "REFUSED: %s line %d: bad checksum (corrupt checkpoint; inspect manually)\n", fn, lineno); return 7; }
            if (lineno == 1) {
                char want[1024]; snprintf(want, sizeof want, "H smale_bb_v2 cfg=%s %s", cfg16, cfgs);
                if (strcmp(payload, want)) { fprintf(stderr, "REFUSED: %s header is for a different configuration:\n  %s\n", fn, payload); return 6; }
                have_header = 1;
            } else {
                int t; unsigned long long pr, F_, E_, L_, o, sy, ex, un; int md; double el, vd; char cf[64]; int pos = 0;
                if (sscanf(payload, "R %d %llu %llu %llu %llu %llu %llu %llu %llu %d %lf %lf cfg=%63s%n", &t, &pr, &F_, &E_, &L_,
                           &o, &sy, &ex, &un, &md, &el, &vd, cf, &pos) != 13 || payload[pos] != 0) {
                    fprintf(stderr, "REFUSED: %s line %d: malformed record\n", fn, lineno); return 7; }
                if (strcmp(cf, cfg16)) { fprintf(stderr, "REFUSED: %s line %d: record from another configuration\n", fn, lineno); return 7; }
                if (t < 0 || t >= ntask) { fprintf(stderr, "REFUSED: %s line %d: id %d out of range\n", fn, lineno, t); return 7; }
                if (done[t]) { fprintf(stderr, "REFUSED: %s line %d: duplicate id %d\n", fn, lineno, t); return 7; }
                unsigned long long leaves = F_ + E_ + L_ + o + sy + ex + un;
                if (pr != 2 * leaves - 1) { fprintf(stderr, "REFUSED: %s line %d: tree count mismatch\n", fn, lineno); return 7; }
                done[t] = 1; ndone++;
                W.tot.processed += pr; W.tot.F += F_; W.tot.E += E_; W.tot.L += L_; W.tot.outside += o;
                W.tot.sym += sy; W.tot.excl += ex; W.tot.unresolved += un; W.tot.vol_done += vd;
                if (md > W.tot.maxdepth) W.tot.maxdepth = md;
            }
            good += len + 1;
        }
        if (good < fl) {
            fprintf(stderr, "note: truncating torn trailing line (%zu bytes) of %s\n", fl - good, fn);
            if (ftruncate(fd, (off_t)good) != 0) { perror("ftruncate"); return 3; }
            full_sync(fd);
        }
        free(fb);
        if (!have_header) {
            char pl[1024], line[1100];
            snprintf(pl, sizeof pl, "H smale_bb_v2 cfg=%s %s", cfg16, cfgs);
            int L = make_line(line, sizeof line, pl);
            write_all(fd, line, (size_t)L); full_sync(fd);
        }
        fprintf(stderr, "smale_bb_v2 %s d=%d r_excl=%.17g excl_mode=%d s=%d seed=%llu cfg=%s: %d tasks (%ld trivial top boxes), %d already done\n",
                SRC_SHA256, d, r, em, s, seed, cfg16, ntask, ntriv, ndone);
        W.P = &P; W.ntask = ntask; W.tc = tc; W.th = th; W.done = done; W.next = 0; W.seed = seed; W.cfg16 = cfg16;
        pthread_mutex_init(&W.mtx, NULL);
        W.done_fd = fd;
        double sp_ = getenv("SMALE_SAMPLE_P") ? atof(getenv("SMALE_SAMPLE_P")) : 0.0;
        if (sp_ > 0) { snprintf(fn, sizeof fn, "%s.samples", pre); W.samp = fopen(fn, "a"); W.samp_p = sp_; }
        snprintf(fn, sizeof fn, "%s.unresolved", pre);
        W.unres_fd = open(fn, O_WRONLY | O_CREAT | O_APPEND, 0644);
        if (W.unres_fd < 0) { perror(fn); return 3; }
        W.t0 = now();
        pthread_t th_[64];
        for (int k = 0; k < nth; k++) pthread_create(&th_[k], NULL, worker, &W);
        for (int k = 0; k < nth; k++) pthread_join(th_[k], NULL);
        double el = now() - W.t0;
        double vtot = (double)ntask; for (int k = 0; k < D; k++) vtot *= 2 * th;
        printf("{\"version\": \"smale_bb_v2\", \"cfg\": \"%s\", \"d\": %d, \"r_excl\": %g, \"excl_mode\": %d, \"s\": %d, \"tasks\": %d, \"processed\": %llu, \"F\": %llu, "
               "\"E\": %llu, \"L\": %llu, \"outside\": %llu, \"symmetry\": %llu, \"excluded\": %llu, \"unresolved\": %llu, "
               "\"maxdepth\": %d, \"vol_frac\": %.15g, \"wall_this_session\": %.1f}\n",
               cfg16, d, r, em, s, ntask, (unsigned long long)W.tot.processed, (unsigned long long)W.tot.F,
               (unsigned long long)W.tot.E, (unsigned long long)W.tot.L, (unsigned long long)W.tot.outside,
               (unsigned long long)W.tot.sym, (unsigned long long)W.tot.excl, (unsigned long long)W.tot.unresolved,
               W.tot.maxdepth, W.tot.vol_done / vtot, el);
        fflush(stdout);
        full_sync(fd); close(fd);
        if (W.samp) fclose(W.samp);
        close(W.unres_fd);
        return W.tot.unresolved ? 2 : 0;
    }
    usage();
    return 1;
}
