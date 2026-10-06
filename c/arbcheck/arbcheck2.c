/*
 * arbcheck2.c -- independent FLINT/arb verifier of a subdivision tree for Smale's
 * mean value conjecture (see SPEC.md and SPEC_W.md in this directory).
 *
 * arbcheck2 = arbcheck (SPEC.md, T-form) + the w-form of SPEC_W.md as an additional,
 * independent way to bound S_i (i >= 2).  Written from SPEC.md, SPEC_W.md and mathematics only.
 *
 * Build (see Makefile, target arbcheck2):
 *   cc -O2 -Wall -ffp-contract=off -fcx-limited-range -include flint_compat.h -I/opt/homebrew/include \
 *      -L/opt/homebrew/lib -o arbcheck2 arbcheck2.c -lflint -lgmp -lmpfr -lpthread -lm
 *
 * Mathematics used (u_1 := 1, mu_k := 1/((k+1)(k+2)), e_k = elementary symmetric):
 *   S_1 = sum_{k=0}^{n-1} (-1)^k mu_k e_k(u_2..u_n)
 *   T_i = sum_{k=0}^{n-1} (-1)^k mu_k u_i^{n-1-k} e_k({u_j : j != i, j = 1..n})
 * (expand prod_{j != i} (u_i - t u_j) in t and integrate (1-t) t^k exactly).
 * w-form (w = 1/u_i, O_i = {2..n} \ {i}):
 *   S_i = int_0^1 (1-t)(1-tw) prod_{j in O_i} (1 - t u_j w) dt
 *       = sum_{k=0}^{n-1} (-1)^k mu_k w^k e_k({1} u {u_j : j in O_i})
 * which is T_i / u_i^{n-1} with w = 1/u_i (the coefficient of u_i^{n-1-k} in T_i is the
 * coefficient of w^k here).
 * The Taylor expansion at the box centre m is obtained exactly (up to ball
 * rounding) from the elementary symmetric functions of all subsets of the
 * centre coordinates, followed by a Taylor shift in delta_i (T-form) or in
 * delta_w = w - 1/m_i (w-form, 1/m_i enclosed in a ball).
 * See README.md and README2.md for the tests used for each leaf type.
 */
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/mag.h>
#include <flint/fmpq.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <complex.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/resource.h>

#define MAXNV 5
#define MAXDIM (2 * MAXNV)
#define MAXN (MAXNV + 1)
#define MAXSUB (1 << MAXNV)
#define NCODE 7
#define MAXDEPTH_TREE 400

enum { C_F = 0, C_E = 1, C_L = 2, C_O = 3, C_S = 4, C_X = 5, C_U = 6 };
static const char *CODE_NAME[NCODE] = {"F", "E", "L", "O", "S", "X", "U"};

/* ------------------------------------------------------------------ */
/* global configuration and constants (read-only once threads start)   */

static int g_d, g_n, g_nv, g_D, g_ALL;
static slong g_prec = 64;
static int g_maxref = 12;
static int g_timing = 0;
static int g_quiet = 0;
static int g_maxfailprint = 20;
static int g_debug = 0;
static int g_wform = 1;              /* use the w-form of SPEC_W.md (--no-wform disables it) */
#define ARBCHECK2_VERSION "arbcheck2 2.0 (SPEC.md T-form + SPEC_W.md w-form)"

static arb_t g_smu[MAXN + 1];       /* (-1)^k mu_k */
static double g_smu_d[MAXN + 1];
static arb_t g_logc, g_nlogc, g_r2excl, g_one;
static double g_logc_d;
static mag_t g_clo, g_cnlo;         /* lower bounds for c and c^n */
static int g_npts;                  /* number of equality points */
static arb_ptr g_peq;               /* g_npts * g_D */
static double *g_peq_d;

static int popc(unsigned x) { return __builtin_popcount(x); }

/* ------------------------------------------------------------------ */
typedef struct {
    double C[MAXDIM], H[MAXDIM];
} box_t;

/*
 * Per-thread workspace. Every arb/acb/mag variable is allocated once (ws_new) and reused;
 * everything that depends on the box is computed lazily and cached until the next box.
 */
typedef struct {
    slong prec;
    const box_t *b;
    int have_basic, have_est;
    uint32_t es_done, fw_done, v_done;
    acb_t m[MAXNV];                              /* centre coordinates */
    mag_t rho[MAXNV], mlo[MAXNV], mhi[MAXNV], Hm[MAXDIM];
    mag_t rhoJ[MAXSUB];                          /* prod_{v in J} rho_v */
    mag_t rpow[MAXNV][MAXN];                     /* rho_v^a */
    acb_t es[MAXSUB][MAXN + 1];                  /* e_r(m_W) */
    acb_t F[MAXSUB][MAXN + 1];                   /* F_W[p] = E1(W, |W|+1-p) */
    /* V[v][a][p] = C(p,a) (-1)^{n-1-p} mu_{n-1-p} m_v^{p-a}, stored interleaved as
       Vc = (Re V, -Im V) and Vs = (Im V, Re V), so that Re(sum V F) = Vc . F and Im(sum V F) = Vs . F */
    arb_t Vc[MAXNV][MAXN][2 * MAXN], Vs[MAXNV][MAXN][2 * MAXN];
    acb_t pw[MAXN], tv;
    /* non-rigorous double predictions (ordering only) */
    double est[MAXN], de[MAXN], lShi[MAXN], lSlo[MAXN], fup[MAXN], Gd[MAXN][MAXDIM];
    /* rigorous models */
    int mstate[MAXN], astate[MAXN], aff[MAXN];
    acb_t cs[MAXN][MAXNV + 1];                   /* c0 and the linear coefficients */
    mag_t Lf[MAXN], Rf[MAXN], c0lo[MAXN];
    arb_t a[MAXN];
    arb_t G[MAXN][MAXDIM];
    mag_t e[MAXN], Shi[MAXN], Slo[MAXN];
    int uc[MAXNV];
    arb_t lu[MAXNV], gu[MAXNV][2];
    acb_t sc;
    arb_t t1, t2, t3, t4, inv2, GS[MAXDIM];
    mag_t m1, m2, m3, m4, LR, c0hi;
    /* w-form of S_f, f = vi+1 >= 1 (SPEC_W.md); variables delta_w = w - 1/m_vi and delta_j, j != vi */
    int usew;                                    /* current pass: 0 = T-form only, 1 = w-form tests */
    uint32_t w0_done;
    int wdom[MAXNV];                             /* |m_v| > rho_v certified */
    acb_t w0[MAXNV], w0sq[MAXNV];                /* balls containing 1/m_v and 1/m_v^2 */
    mag_t rhow[MAXNV], etab[MAXNV];              /* upper bounds of |delta_w| and of |eta| on the box */
    mag_t rwpow[MAXNV][MAXN];                    /* rhow^a */
    /* VW[v][a][p] = C(n-1-p, a) (-1)^{n-1-p} mu_{n-1-p} w0^{n-1-p-a}, p = 0..n-1-a */
    acb_t VW[MAXNV][MAXN][MAXN];
    int wstate[MAXN], wok[MAXN], wastate[MAXN], waff[MAXN];
    acb_t cw[MAXN][MAXNV + 1];                   /* c0 and linear coefficients (slot 1+vi: delta_w) */
    mag_t LW[MAXN], RW[MAXN], c0loW[MAXN], ShiW[MAXN], SloW[MAXN], eW[MAXN];
    arb_t aW[MAXN], GW[MAXN][MAXDIM];
    acb_t wc, wg;
    mag_t wm1, wm2, wm3, wm4;
    /* statistics */
    uint64_t subboxes, nmodel, naff, nwmodel, nwaff, wsub;
} ws_t;

static ws_t *ws_new(slong prec)
{
    ws_t *w = calloc(1, sizeof(ws_t));
    int i, j, k;
    w->prec = prec;
    for (i = 0; i < MAXNV; i++) {
        acb_init(w->m[i]); mag_init(w->rho[i]); mag_init(w->mlo[i]); mag_init(w->mhi[i]);
        for (j = 0; j < MAXN; j++) mag_init(w->rpow[i][j]);
        for (j = 0; j < MAXN; j++) for (k = 0; k < 2 * MAXN; k++) { arb_init(w->Vc[i][j][k]); arb_init(w->Vs[i][j][k]); }
        arb_init(w->lu[i]); arb_init(w->gu[i][0]); arb_init(w->gu[i][1]);
    }
    for (i = 0; i < MAXN; i++) acb_init(w->pw[i]);
    acb_init(w->tv);
    for (i = 0; i < MAXDIM; i++) { mag_init(w->Hm[i]); arb_init(w->GS[i]); }
    for (i = 0; i < MAXSUB; i++) {
        mag_init(w->rhoJ[i]);
        for (j = 0; j <= MAXN; j++) acb_init(w->es[i][j]);
        for (j = 0; j <= MAXN; j++) acb_init(w->F[i][j]);
    }
    for (i = 0; i < MAXN; i++) {
        arb_init(w->a[i]);
        for (j = 0; j < MAXDIM; j++) arb_init(w->G[i][j]);
        for (j = 0; j <= MAXNV; j++) acb_init(w->cs[i][j]);
        mag_init(w->e[i]); mag_init(w->Shi[i]); mag_init(w->Slo[i]);
        mag_init(w->Lf[i]); mag_init(w->Rf[i]); mag_init(w->c0lo[i]);
    }
    acb_init(w->sc);
    for (i = 0; i < MAXNV; i++) {
        acb_init(w->w0[i]); acb_init(w->w0sq[i]); mag_init(w->rhow[i]); mag_init(w->etab[i]);
        for (j = 0; j < MAXN; j++) {
            mag_init(w->rwpow[i][j]);
            for (k = 0; k < MAXN; k++) acb_init(w->VW[i][j][k]);
        }
    }
    for (i = 0; i < MAXN; i++) {
        for (j = 0; j <= MAXNV; j++) acb_init(w->cw[i][j]);
        mag_init(w->LW[i]); mag_init(w->RW[i]); mag_init(w->c0loW[i]); mag_init(w->ShiW[i]); mag_init(w->SloW[i]); mag_init(w->eW[i]);
        arb_init(w->aW[i]);
        for (j = 0; j < MAXDIM; j++) arb_init(w->GW[i][j]);
    }
    acb_init(w->wc); acb_init(w->wg);
    mag_init(w->wm1); mag_init(w->wm2); mag_init(w->wm3); mag_init(w->wm4);
    arb_init(w->t1); arb_init(w->t2); arb_init(w->t3); arb_init(w->t4); arb_init(w->inv2);
    mag_init(w->m1); mag_init(w->m2); mag_init(w->m3); mag_init(w->m4); mag_init(w->LR); mag_init(w->c0hi);
    return w;
}

static void ws_free(ws_t *w)
{
    int i, j, k;
    for (i = 0; i < MAXNV; i++) {
        acb_clear(w->m[i]); mag_clear(w->rho[i]); mag_clear(w->mlo[i]); mag_clear(w->mhi[i]);
        for (j = 0; j < MAXN; j++) mag_clear(w->rpow[i][j]);
        for (j = 0; j < MAXN; j++) for (k = 0; k < 2 * MAXN; k++) { arb_clear(w->Vc[i][j][k]); arb_clear(w->Vs[i][j][k]); }
        arb_clear(w->lu[i]); arb_clear(w->gu[i][0]); arb_clear(w->gu[i][1]);
    }
    for (i = 0; i < MAXN; i++) acb_clear(w->pw[i]);
    acb_clear(w->tv);
    for (i = 0; i < MAXDIM; i++) { mag_clear(w->Hm[i]); arb_clear(w->GS[i]); }
    for (i = 0; i < MAXSUB; i++) {
        mag_clear(w->rhoJ[i]);
        for (j = 0; j <= MAXN; j++) acb_clear(w->es[i][j]);
        for (j = 0; j <= MAXN; j++) acb_clear(w->F[i][j]);
    }
    for (i = 0; i < MAXN; i++) {
        arb_clear(w->a[i]);
        for (j = 0; j < MAXDIM; j++) arb_clear(w->G[i][j]);
        for (j = 0; j <= MAXNV; j++) acb_clear(w->cs[i][j]);
        mag_clear(w->e[i]); mag_clear(w->Shi[i]); mag_clear(w->Slo[i]);
        mag_clear(w->Lf[i]); mag_clear(w->Rf[i]); mag_clear(w->c0lo[i]);
    }
    acb_clear(w->sc);
    for (i = 0; i < MAXNV; i++) {
        acb_clear(w->w0[i]); acb_clear(w->w0sq[i]); mag_clear(w->rhow[i]); mag_clear(w->etab[i]);
        for (j = 0; j < MAXN; j++) {
            mag_clear(w->rwpow[i][j]);
            for (k = 0; k < MAXN; k++) acb_clear(w->VW[i][j][k]);
        }
    }
    for (i = 0; i < MAXN; i++) {
        for (j = 0; j <= MAXNV; j++) acb_clear(w->cw[i][j]);
        mag_clear(w->LW[i]); mag_clear(w->RW[i]); mag_clear(w->c0loW[i]); mag_clear(w->ShiW[i]); mag_clear(w->SloW[i]); mag_clear(w->eW[i]);
        arb_clear(w->aW[i]);
        for (j = 0; j < MAXDIM; j++) arb_clear(w->GW[i][j]);
    }
    acb_clear(w->wc); acb_clear(w->wg);
    mag_clear(w->wm1); mag_clear(w->wm2); mag_clear(w->wm3); mag_clear(w->wm4);
    arb_clear(w->t1); arb_clear(w->t2); arb_clear(w->t3); arb_clear(w->t4); arb_clear(w->inv2);
    mag_clear(w->m1); mag_clear(w->m2); mag_clear(w->m3); mag_clear(w->m4); mag_clear(w->LR); mag_clear(w->c0hi);
    free(w);
}

/* attach a box to the workspace; everything else is computed lazily */
static void ws_set_box(ws_t *w, const box_t *b)
{
    int k;
    w->b = b;
    w->have_basic = w->have_est = 0;
    w->es_done = w->fw_done = w->v_done = 0;
    for (k = 0; k < MAXN; k++) w->mstate[k] = w->astate[k] = w->aff[k] = 0;   /* per-box model caches (review issue 2) */
    w->w0_done = 0;
    w->usew = 0;
    for (k = 0; k < MAXN; k++) w->wstate[k] = w->wok[k] = w->wastate[k] = w->waff[k] = 0;
}

static void ensure_basic(ws_t *w)
{
    int v, k, a, mask;
    if (w->have_basic) return;
    w->have_basic = 1;
    for (v = 0; v < g_nv; v++) {
        acb_set_d_d(w->m[v], w->b->C[2 * v], w->b->C[2 * v + 1]);   /* exact */
        mag_set_d(w->m1, w->b->H[2 * v]);
        mag_set_d(w->m2, w->b->H[2 * v + 1]);
        mag_hypot(w->rho[v], w->m1, w->m2);      /* upper bound of |delta_v| on the box */
        acb_get_mag_lower(w->mlo[v], w->m[v]);
        acb_get_mag(w->mhi[v], w->m[v]);
        mag_one(w->rpow[v][0]);
        for (a = 1; a < g_n; a++) mag_mul(w->rpow[v][a], w->rpow[v][a - 1], w->rho[v]);
        w->uc[v] = 0;
    }
    for (k = 0; k < g_D; k++) mag_set_d(w->Hm[k], w->b->H[k]);
    mag_one(w->rhoJ[0]);
    for (mask = 1; mask <= g_ALL; mask++) {
        v = __builtin_ctz(mask);
        mag_mul(w->rhoJ[mask], w->rhoJ[mask & (mask - 1)], w->rho[v]);
    }
    for (k = 0; k < g_n; k++) w->mstate[k] = w->astate[k] = w->aff[k] = 0;
    w->w0_done = 0;
    for (k = 0; k < g_n; k++) w->wstate[k] = w->wok[k] = w->wastate[k] = w->waff[k] = 0;
}

/* elementary symmetric functions e_r(m_W) of subsets W of the centre coordinates (lazy per subset) */
static void ensure_es_mask(ws_t *w, int mask)
{
    int v, prev, k, r;
    slong prec = w->prec;
    if (w->es_done & (1u << mask)) return;
    if (mask == 0) {
        acb_one(w->es[0][0]);
        w->es_done |= 1u;
        return;
    }
    v = __builtin_ctz(mask);
    prev = mask & (mask - 1);
    ensure_es_mask(w, prev);
    k = popc(prev);
    acb_mul(w->es[mask][k + 1], w->es[prev][k], w->m[v], prec);
    for (r = k; r >= 1; r--) {
        acb_mul(w->es[mask][r], w->es[prev][r - 1], w->m[v], prec);
        acb_add(w->es[mask][r], w->es[mask][r], w->es[prev][r], prec);
    }
    acb_one(w->es[mask][0]);
    w->es_done |= 1u << mask;
}

/* F_W[p] = E1(W, |W|+1-p), p = 0..|W|+1, where E1(W, r) = e_r(m_W) + e_{r-1}(m_W) = e_r({1} u m_W) */
static void ensure_fw(ws_t *w, int W)
{
    int wn, deg, p, r;
    slong prec = w->prec;
    if (w->fw_done & (1u << W)) return;
    ensure_es_mask(w, W);
    wn = popc(W);
    deg = wn + 1;
    for (p = 0; p <= deg; p++) {
        r = deg - p;
        if (r == 0) acb_one(w->F[W][p]);
        else if (r == wn + 1) acb_set(w->F[W][p], w->es[W][wn]);
        else acb_add(w->F[W][p], w->es[W][r], w->es[W][r - 1], prec);
    }
    w->fw_done |= 1u << W;
}

/* V[v][a][p] = K[a][p] m_v^{p-a},  K[a][p] = C(p,a) (-1)^{n-1-p} mu_{n-1-p} */
static arb_t g_K[MAXN][MAXN];
static void ensure_V(ws_t *w, int vi)
{
    int a, p, n = g_n;
    slong prec = w->prec;
    if (w->v_done & (1u << vi)) return;
    acb_one(w->pw[0]);
    for (p = 1; p < n; p++) acb_mul(w->pw[p], w->pw[p - 1], w->m[vi], prec);
    for (a = 0; a < n; a++)
        for (p = a; p < n; p++) {
            acb_mul_arb(w->tv, w->pw[p - a], g_K[a][p], prec);
            arb_set(w->Vc[vi][a][2 * p], acb_realref(w->tv));
            arb_neg(w->Vc[vi][a][2 * p + 1], acb_imagref(w->tv));
            arb_set(w->Vs[vi][a][2 * p], acb_imagref(w->tv));
            arb_set(w->Vs[vi][a][2 * p + 1], acb_realref(w->tv));
        }
    w->v_done |= 1u << vi;
}

/*
 * Non-rigorous double-precision mirror of model()/model_aff() below, used ONLY to choose
 * the order in which functions / pairs are tried (never to decide a claim):
 *   est[f] ~ log|S_f(m)|, Gd[f] ~ gradient, de[f] ~ affine error (inf if no affine model),
 *   lShi[f], lSlo[f] ~ logs of the magnitude bounds, fup[f] ~ predicted upper bound of log|S_f|.
 */
static inline double cabs2(double complex z) { return creal(z) * creal(z) + cimag(z) * cimag(z); }

static void ensure_est(ws_t *w)
{
    double complex mm[MAXNV], esd[MAXSUB][MAXN + 1], cf[MAXN];
    double rho[MAXNV], am[MAXNV], rJ[MAXSUB], rp[MAXNV][MAXN];
    int mask, v, prev, k, r, f, n = g_n, ALL = g_ALL;
    const double *H = w->b->H;
    if (w->have_est) return;
    w->have_est = 1;
    for (v = 0; v < g_nv; v++) {
        mm[v] = w->b->C[2 * v] + w->b->C[2 * v + 1] * _Complex_I;
        rho[v] = sqrt(H[2 * v] * H[2 * v] + H[2 * v + 1] * H[2 * v + 1]);
        am[v] = sqrt(cabs2(mm[v]));
        rp[v][0] = 1;
        for (k = 1; k < n; k++) rp[v][k] = rp[v][k - 1] * rho[v];
    }
    esd[0][0] = 1;
    rJ[0] = 1;
    for (mask = 1; mask <= ALL; mask++) {
        v = __builtin_ctz(mask);
        prev = mask & (mask - 1);
        k = popc(prev);
        esd[mask][k + 1] = esd[prev][k] * mm[v];
        for (r = k; r >= 1; r--) esd[mask][r] = esd[prev][r] + esd[prev][r - 1] * mm[v];
        esd[mask][0] = 1;
        rJ[mask] = rJ[prev] * rho[v];
    }
    for (f = 0; f < n; f++) {
        double L = 0, R = 0, ac0, Thi, Tlo, q, e, A, sp = 0;
        double complex c0 = 0, lin[MAXNV];
        int vi = f - 1, aff = 1;
        if (f == 0) {
            int J;
            for (J = 0; J <= ALL; J++) {
                int W = ALL ^ J, j = popc(J);
                double complex c = 0;
                for (r = 0; r <= popc(W); r++) c += g_smu_d[r + j] * esd[W][r];
                if (j == 0) c0 = c;
                else if (j == 1) { lin[__builtin_ctz(J)] = c; L += sqrt(cabs2(c)) * rJ[J]; }
                else R += sqrt(cabs2(c)) * rJ[J];
            }
        } else {
            int O = ALL ^ (1 << vi), J;
            for (J = O;; J = (J - 1) & O) {
                int W = O ^ J, j = popc(J), wn = popc(W), deg = n - 1 - j, p, s2, a;
                for (p = 0; p <= deg; p++) {
                    int kk = n - 1 - p, rr = kk - j;
                    double complex E1 = (rr <= wn ? esd[W][rr] : 0) + (rr >= 1 ? esd[W][rr - 1] : 0);
                    cf[p] = g_smu_d[kk] * E1;
                }
                for (s2 = 0; s2 < deg; s2++)
                    for (p = deg - 1; p >= s2; p--) cf[p] += mm[vi] * cf[p + 1];
                for (a = 0; a <= deg; a++) {
                    int tot = j + a;
                    double t;
                    if (!tot) { c0 = cf[0]; continue; }
                    if (tot == 1) lin[a ? vi : __builtin_ctz(J)] = cf[a];
                    t = sqrt(cabs2(cf[a])) * rJ[J] * rp[vi][a];
                    if (tot == 1) L += t; else R += t;
                }
                if (J == 0) break;
            }
        }
        ac0 = sqrt(cabs2(c0));
        Thi = ac0 + L + R;
        Tlo = ac0 - L - R;
        w->lShi[f] = log(Thi);
        w->lSlo[f] = Tlo > 0 ? log(Tlo) : -INFINITY;
        A = log(ac0);
        q = (L + R) / ac0;
        e = (q < 1) ? R / ac0 + q * q / (2 * (1 - q)) : INFINITY;
        if (!(q < 1)) aff = 0;
        for (v = 0; v < g_nv; v++) {
            double complex g = lin[v] / c0;
            w->Gd[f][2 * v] = creal(g);
            w->Gd[f][2 * v + 1] = -cimag(g);
        }
        if (f) {
            double pp = rho[vi] / am[vi], ulo = am[vi] - rho[vi];
            A -= (n - 1) * log(am[vi]);
            w->Gd[f][2 * vi] -= (n - 1) * creal(mm[vi]) / (am[vi] * am[vi]);
            w->Gd[f][2 * vi + 1] -= (n - 1) * cimag(mm[vi]) / (am[vi] * am[vi]);
            if (pp < 1) e += (n - 1) * pp * pp / (2 * (1 - pp)); else aff = 0;
            w->lShi[f] = ulo > 0 ? w->lShi[f] - (n - 1) * log(ulo) : INFINITY;
            w->lSlo[f] -= (n - 1) * log(am[vi] + rho[vi]);
        }
        for (k = 0; k < g_D; k++) sp += fabs(w->Gd[f][k]) * H[k];
        w->est[f] = isnan(A) ? INFINITY : A;
        w->de[f] = (aff && isfinite(e) && isfinite(sp)) ? e : INFINITY;
        w->fup[f] = w->lShi[f];
        if (isfinite(w->de[f]) && A + sp + e < w->fup[f]) w->fup[f] = A + sp + e;
        if (isnan(w->fup[f])) w->fup[f] = INFINITY;
        if (isnan(w->lShi[f])) w->lShi[f] = INFINITY;
        if (isnan(w->lSlo[f])) w->lSlo[f] = -INFINITY;
    }
}

/*
 * Rigorous model of S_f on the box (f = 0: S_1; f = v+1: S_{v+2}, variable v).
 * model():     Taylor coefficients at the centre (balls), c0 and linear coefficients kept;
 *              L = sum |c_lin| rho, R = sum_{order>=2} |c_alpha| rho^alpha (upper bounds);
 *              magnitude bounds Slo <= |S_f| <= Shi at all nondegenerate points of the box.
 * model_aff(): if possible, log|S_f(x)| in a + G.(x - m) +- e for all x in the box.
 */
static void model(ws_t *w, int f)
{
    slong prec = w->prec;
    int n = g_n, ALL = g_ALL;
    acb_ptr c0 = w->cs[f][0];
    mag_ptr L = w->Lf[f], R = w->Rf[f];
    if (w->mstate[f]) return;
    ensure_basic(w);
    w->nmodel++;
    w->mstate[f] = 1;
    mag_zero(L);
    mag_zero(R);

    if (f == 0) {
        /* S_1 = sum_J delta^J sum_r (-1)^{r+|J|} mu_{r+|J|} e_r(m_{ALL \ J}) */
        int J;
        for (J = 0; J <= ALL; J++) {
            int W = ALL ^ J, j = popc(J), wn = popc(W);
            acb_ptr c = (j == 0) ? c0 : (j == 1) ? w->cs[0][1 + __builtin_ctz(J)] : w->sc;
            ensure_es_mask(w, W);
            arb_dot(acb_realref(c), NULL, 0, g_smu[j], 1, acb_realref(w->es[W][0]), 2, wn + 1, prec);
            arb_dot(acb_imagref(c), NULL, 0, g_smu[j], 1, acb_imagref(w->es[W][0]), 2, wn + 1, prec);
            if (j == 0) continue;
            acb_get_mag(w->m1, c);
            mag_addmul(j == 1 ? L : R, w->m1, w->rhoJ[J]);
        }
    } else {
        /* T_i = sum_{J subset O} delta^J A_J(u_i),  A_J(u) = sum_p (-1)^{n-1-p} mu_{n-1-p} E1(W, n-1-p-|J|) u^p,
           Taylor coefficient of delta_i^a: c_{J,a} = sum_{p>=a} V[i][a][p] F_W[p]  (W = O \ J) */
        int vi = f - 1, O = ALL ^ (1 << vi), J;
        ensure_V(w, vi);
        for (J = O;; J = (J - 1) & O) {
            int W = O ^ J, j = popc(J), deg = n - 1 - j, a;
            ensure_fw(w, W);
            for (a = 0; a <= deg; a++) {
                int tot = j + a, len = 2 * (deg - a + 1);
                acb_ptr c = (tot == 0) ? c0 : (tot == 1) ? w->cs[f][1 + (a ? vi : __builtin_ctz(J))] : w->sc;
                arb_srcptr fw = acb_realref(w->F[W][a]);   /* interleaved (Re, Im) from p = a */
                arb_dot(acb_realref(c), NULL, 0, w->Vc[vi][a][2 * a], 1, fw, 1, len, prec);
                arb_dot(acb_imagref(c), NULL, 0, w->Vs[vi][a][2 * a], 1, fw, 1, len, prec);
                if (tot == 0) continue;
                acb_get_mag(w->m1, c);
                mag_mul(w->m1, w->m1, w->rhoJ[J]);
                if (a) mag_mul(w->m1, w->m1, w->rpow[vi][a]);
                mag_add(tot == 1 ? L : R, tot == 1 ? L : R, w->m1);
            }
            if (J == 0) break;
        }
    }

    /* magnitude (interval) bounds */
    acb_get_mag(w->c0hi, c0);
    acb_get_mag_lower(w->c0lo[f], c0);
    mag_add(w->LR, L, R);
    mag_add(w->Shi[f], w->c0hi, w->LR);            /* sup |T| */
    mag_sub_lower(w->Slo[f], w->c0lo[f], w->LR);   /* inf |T| */
    if (f) {
        int vi = f - 1;
        /* |u_i| in [mlo - rho, mhi + rho] */
        mag_sub_lower(w->m1, w->mlo[vi], w->rho[vi]);
        mag_add(w->m2, w->mhi[vi], w->rho[vi]);
        if (mag_is_zero(w->m1)) mag_inf(w->Shi[f]);
        else {
            mag_pow_ui_lower(w->m1, w->m1, n - 1);
            mag_div(w->Shi[f], w->Shi[f], w->m1);
        }
        mag_pow_ui(w->m2, w->m2, n - 1);
        mag_div_lower(w->Slo[f], w->Slo[f], w->m2);
    }
}

static int model_aff(ws_t *w, int f)
{
    slong prec = w->prec;
    int n = g_n, v;
    acb_srcptr c0 = w->cs[f][0];
    if (w->astate[f]) return w->aff[f];
    model(w, f);
    w->astate[f] = 1;
    w->aff[f] = 0;
    w->naff++;
    if (mag_is_zero(w->c0lo[f])) return 0;
    mag_add(w->LR, w->Lf[f], w->Rf[f]);
    mag_div(w->m1, w->LR, w->c0lo[f]);              /* q >= |w| where T = c0 (1 + w) */
    if (mag_cmp_2exp_si(w->m1, 0) >= 0) return 0;
    if (f) {
        int vi = f - 1;
        if (mag_is_zero(w->mlo[vi])) return 0;
        mag_div(w->m3, w->rho[vi], w->mlo[vi]);      /* p >= |delta_i / m_i| */
        if (mag_cmp_2exp_si(w->m3, 0) >= 0) return 0;
    }
    /* e_T = R/|c0| + q^2 / (2 (1-q))  (|log(1+w) - w| <= q^2/(2(1-q)) for |w| <= q < 1) */
    mag_one(w->m2);
    mag_sub_lower(w->m2, w->m2, w->m1);             /* 1 - q, lower bound */
    mag_mul(w->m4, w->m1, w->m1);
    mag_div(w->m4, w->m4, w->m2);
    mag_mul_2exp_si(w->m4, w->m4, -1);
    mag_div(w->e[f], w->Rf[f], w->c0lo[f]);
    mag_add(w->e[f], w->e[f], w->m4);

    arb_log_hypot(w->a[f], acb_realref(c0), acb_imagref(c0), prec);
    arb_sqr(w->t1, acb_realref(c0), prec);
    arb_addmul(w->t1, acb_imagref(c0), acb_imagref(c0), prec);
    arb_inv(w->inv2, w->t1, prec);

    for (v = 0; v < g_nv; v++) {
        acb_srcptr c = w->cs[f][1 + v];
        /* gamma = c / c0 = c conj(c0) / |c0|^2 ;  Re(gamma delta) = Re(gamma) dx - Im(gamma) dy */
        arb_mul(w->t1, acb_realref(c), acb_realref(c0), prec);
        arb_addmul(w->t1, acb_imagref(c), acb_imagref(c0), prec);
        arb_mul(w->G[f][2 * v], w->t1, w->inv2, prec);
        arb_mul(w->t2, acb_realref(c), acb_imagref(c0), prec);
        arb_submul(w->t2, acb_imagref(c), acb_realref(c0), prec);   /* = -Im(c conj c0) */
        arb_mul(w->G[f][2 * v + 1], w->t2, w->inv2, prec);
    }
    if (f) {
        int vi = f - 1;
        if (!w->uc[vi]) {
            acb_srcptr mi = w->m[vi];
            w->uc[vi] = 1;
            arb_log_hypot(w->lu[vi], acb_realref(mi), acb_imagref(mi), prec);
            arb_sqr(w->t1, acb_realref(mi), prec);
            arb_addmul(w->t1, acb_imagref(mi), acb_imagref(mi), prec);
            arb_div(w->gu[vi][0], acb_realref(mi), w->t1, prec);
            arb_div(w->gu[vi][1], acb_imagref(mi), w->t1, prec);
        }
        /* log|u_i| = log|m_i| + Re(delta_i/m_i) +- p^2/(2(1-p)) */
        arb_submul_ui(w->a[f], w->lu[vi], n - 1, prec);
        arb_submul_ui(w->G[f][2 * vi], w->gu[vi][0], n - 1, prec);
        arb_submul_ui(w->G[f][2 * vi + 1], w->gu[vi][1], n - 1, prec);
        mag_one(w->m2);
        mag_sub_lower(w->m2, w->m2, w->m3);
        mag_mul(w->m4, w->m3, w->m3);
        mag_div(w->m4, w->m4, w->m2);
        mag_mul_2exp_si(w->m4, w->m4, -1);
        mag_mul_ui(w->m4, w->m4, n - 1);
        mag_add(w->e[f], w->e[f], w->m4);
    }
    w->aff[f] = 1;
    return 1;
}


/* ------------------------------------------------------------------ */
/* w-form of S_i, i >= 2 (SPEC_W.md)                                    */
/*
 * With w = 1/u_i and O = {2..n} \ {i}:
 *   S_i = sum_{k=0}^{n-1} (-1)^k mu_k w^k e_k({1} u {u_j : j in O}).
 * With u_j = m_j + delta_j (j in O), the coefficient of delta^J (J subset O) is the polynomial
 *   B_J(w) = sum_k (-1)^k mu_k E1(W, k-|J|) w^k,   W = O \ J,  E1(W, r) = e_r({1} u m_W),
 * and its Taylor coefficients at w0 = 1/m_i are
 *   c_{J,a} = sum_{k >= max(a,|J|)} C(k,a) (-1)^k mu_k w0^{k-a} E1(W, k-|J|).
 * With p = n-1-k:  E1(W, k-|J|) = F_W[p]  and  c_{J,a} = sum_{p=0}^{n-1-max(a,|J|)} VW[a][p] F_W[p].
 * The exact centre w0 = 1/m_i is not a double; it is enclosed in a ball, so every coefficient ball
 * contains the exact Taylor coefficient at the exact w0, and all bounds below refer to the exact
 * expansion point.
 */
static arb_t g_KW[MAXN][MAXN];            /* g_KW[a][p] = C(n-1-p, a) (-1)^{n-1-p} mu_{n-1-p} */

/*
 * Per box and variable vi: certify the domain |m_i| > rho_i and compute
 *   rho_w >= rho_i / (|m_i| (|m_i| - rho_i))         (>= |1/u_i - 1/m_i| on the box)
 *   etab  >= rho_i^2 / (|m_i|^2 (|m_i| - rho_i))     (>= |eta|, eta = delta_w + delta_i/m_i^2)
 * using mlo <= |m_i| (both are decreasing in |m_i| for |m_i| > rho_i) and rho >= sup |delta_i|.
 */
static int ensure_w0(ws_t *w, int vi)
{
    int a, p, n = g_n;
    slong prec = w->prec;
    if (w->w0_done & (1u << vi)) return w->wdom[vi];
    ensure_basic(w);
    w->w0_done |= 1u << vi;
    w->wdom[vi] = 0;
    if (mag_cmp(w->mlo[vi], w->rho[vi]) <= 0) return 0;      /* box may meet u_i = 0: w-form not used */
    mag_sub_lower(w->wm1, w->mlo[vi], w->rho[vi]);           /* <= |m_i| - rho_i <= inf |u_i| */
    if (mag_is_zero(w->wm1)) return 0;
    mag_mul_lower(w->wm2, w->mlo[vi], w->wm1);               /* <= |m_i| (|m_i| - rho_i) */
    if (mag_is_zero(w->wm2)) return 0;
    mag_div(w->rhow[vi], w->rho[vi], w->wm2);
    mag_mul_lower(w->wm3, w->wm2, w->mlo[vi]);               /* <= |m_i|^2 (|m_i| - rho_i) */
    if (mag_is_zero(w->wm3)) return 0;
    mag_mul(w->wm4, w->rho[vi], w->rho[vi]);
    mag_div(w->etab[vi], w->wm4, w->wm3);
    if (mag_is_inf(w->rhow[vi]) || mag_is_inf(w->etab[vi])) return 0;
    mag_one(w->rwpow[vi][0]);
    for (a = 1; a < n; a++) mag_mul(w->rwpow[vi][a], w->rwpow[vi][a - 1], w->rhow[vi]);
    acb_inv(w->w0[vi], w->m[vi], prec);                      /* ball containing 1/m_i (m_i exact, != 0) */
    acb_sqr(w->w0sq[vi], w->w0[vi], prec);                   /* ball containing 1/m_i^2 */
    if (!acb_is_finite(w->w0[vi]) || !acb_is_finite(w->w0sq[vi])) return 0;
    acb_one(w->pw[0]);
    for (p = 1; p < n; p++) acb_mul(w->pw[p], w->pw[p - 1], w->w0[vi], prec);
    for (a = 0; a < n; a++)
        for (p = 0; p <= n - 1 - a; p++)
            acb_mul_arb(w->VW[vi][a][p], w->pw[n - 1 - p - a], g_KW[a][p], prec);
    w->wdom[vi] = 1;
    return 1;
}

/*
 * w-form model of S_f (f >= 1): coefficient balls, L = sum_{order 1} |c| rho^alpha,
 * R = sum_{order >= 2} |c| rho^alpha (rho_w for delta_w), and SloW <= |S_f| <= ShiW at every
 * point of the box (S_f is defined on the whole box because |u_i| >= |m_i| - rho_i > 0).
 * Returns 0 if the w-form does not apply on this box.
 */
static int model_w(ws_t *w, int f)
{
    slong prec = w->prec;
    int n = g_n, vi, O, J;
    acb_ptr c0;
    mag_ptr L, R;
    if (f == 0) return 0;
    if (w->wstate[f]) return w->wok[f];
    w->wstate[f] = 1;
    w->wok[f] = 0;
    vi = f - 1;
    if (!ensure_w0(w, vi)) return 0;
    w->nwmodel++;
    O = g_ALL ^ (1 << vi);
    c0 = w->cw[f][0];
    L = w->LW[f];
    R = w->RW[f];
    mag_zero(L);
    mag_zero(R);
    for (J = O;; J = (J - 1) & O) {
        int W = O ^ J, j = popc(J), a;
        ensure_fw(w, W);
        for (a = 0; a < n; a++) {
            int tot = j + a, len = n - (a > j ? a : j);
            acb_ptr c = (tot == 0) ? c0 : (tot == 1) ? w->cw[f][1 + (a ? vi : __builtin_ctz(J))] : w->wc;
            acb_dot(c, NULL, 0, w->VW[vi][a][0], 1, w->F[W][0], 1, len, prec);
            if (tot == 0) continue;
            acb_get_mag(w->wm1, c);
            mag_mul(w->wm1, w->wm1, w->rhoJ[J]);
            if (a) mag_mul(w->wm1, w->wm1, w->rwpow[vi][a]);
            mag_add(tot == 1 ? L : R, tot == 1 ? L : R, w->wm1);
        }
        if (J == 0) break;
    }
    acb_get_mag(w->wm2, c0);
    acb_get_mag_lower(w->c0loW[f], c0);
    mag_add(w->wm3, L, R);
    mag_add(w->ShiW[f], w->wm2, w->wm3);            /* sup |S_f| */
    mag_sub_lower(w->SloW[f], w->c0loW[f], w->wm3); /* inf |S_f| */
    w->wok[f] = 1;
    return 1;
}

/*
 * w-form affine model: for every x in the box,
 *   log|S_f(x)| in aW + GW.(x - m) +- eW,
 *   aW = log|c0|,  eW = R/|c0| + q^2/(2(1-q)) + |gamma_w| etab,  q = (L+R)/|c0| < 1,
 * where gamma = c_lin / c0, and the delta_w term is made affine with
 *   delta_w = -delta_i / m_i^2 + eta,  |eta| <= etab:
 *   Re(gamma_w delta_w) = Re(-gamma_w w0^2 delta_i) +- |gamma_w| etab.
 */
static int model_waff(ws_t *w, int f)
{
    slong prec = w->prec;
    int v, vi;
    acb_srcptr c0;
    if (f == 0) return 0;
    if (w->wastate[f]) return w->waff[f];
    w->wastate[f] = 1;
    w->waff[f] = 0;
    if (!model_w(w, f)) return 0;
    w->nwaff++;
    vi = f - 1;
    c0 = w->cw[f][0];
    if (mag_is_zero(w->c0loW[f])) return 0;
    mag_add(w->wm1, w->LW[f], w->RW[f]);
    mag_div(w->wm1, w->wm1, w->c0loW[f]);            /* q >= |z| where S = c0 (1 + z) */
    if (mag_cmp_2exp_si(w->wm1, 0) >= 0) return 0;   /* need q < 1 */
    mag_one(w->wm2);
    mag_sub_lower(w->wm2, w->wm2, w->wm1);           /* 1 - q, lower bound */
    if (mag_is_zero(w->wm2)) return 0;
    mag_mul(w->wm4, w->wm1, w->wm1);
    mag_div(w->wm4, w->wm4, w->wm2);
    mag_mul_2exp_si(w->wm4, w->wm4, -1);             /* >= q^2 / (2 (1 - q)) */
    mag_div(w->eW[f], w->RW[f], w->c0loW[f]);         /* >= R / |c0| */
    mag_add(w->eW[f], w->eW[f], w->wm4);
    arb_log_hypot(w->aW[f], acb_realref(c0), acb_imagref(c0), prec);
    for (v = 0; v < g_nv; v++) {
        acb_div(w->wg, w->cw[f][1 + v], c0, prec);   /* gamma_v (v = vi: gamma_w) */
        if (v == vi) {
            acb_get_mag(w->wm1, w->wg);
            mag_mul(w->wm1, w->wm1, w->etab[vi]);     /* >= |gamma_w| |eta| */
            mag_add(w->eW[f], w->eW[f], w->wm1);
            acb_mul(w->wg, w->wg, w->w0sq[vi], prec);
            acb_neg(w->wg, w->wg);                    /* coefficient of delta_i: -gamma_w / m_i^2 */
        }
        if (!acb_is_finite(w->wg)) return 0;
        /* Re(g delta) = Re(g) dx - Im(g) dy */
        arb_set(w->GW[f][2 * v], acb_realref(w->wg));
        arb_neg(w->GW[f][2 * v + 1], acb_imagref(w->wg));
    }
    if (!arb_is_finite(w->aW[f]) || mag_is_inf(w->eW[f])) return 0;
    w->waff[f] = 1;
    return 1;
}

/* best certified magnitude bounds of |S_f| available in the current pass (T-form; w-form if usew) */
static mag_srcptr Shi_of(ws_t *w, int f)
{
    if (w->usew && f && model_w(w, f) && mag_cmp(w->ShiW[f], w->Shi[f]) < 0) return w->ShiW[f];
    return w->Shi[f];
}

static mag_srcptr Slo_of(ws_t *w, int f)
{
    if (w->usew && f && model_w(w, f) && mag_cmp(w->SloW[f], w->Slo[f]) > 0) return w->SloW[f];
    return w->Slo[f];
}

/* affine model s of S_f: s = 0 T-form (model_aff), s = 1 w-form (model_waff). Returns 0 if unavailable. */
static int aff_get(ws_t *w, int f, int s, arb_srcptr *a, arb_srcptr *G, mag_srcptr *e)
{
    if (s == 0) {
        if (!model_aff(w, f)) return 0;
        *a = w->a[f]; *G = w->G[f][0]; *e = w->e[f];
    } else {
        if (!model_waff(w, f)) return 0;
        *a = w->aW[f]; *G = w->GW[f][0]; *e = w->eW[f];
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* the claims                                                           */

static int test_F(ws_t *w)
{
    int ord[MAXN], i, j, k, f;
    ensure_basic(w);
    ensure_est(w);
    for (i = 0; i < g_n; i++) ord[i] = i;
    for (i = 1; i < g_n; i++)
        for (j = i; j > 0 && w->fup[ord[j]] < w->fup[ord[j - 1]]; j--) { int t = ord[j]; ord[j] = ord[j - 1]; ord[j - 1] = t; }
    for (i = 0; i < g_n; i++) {
        f = ord[i];
        if (!(w->est[f] < g_logc_d + 1e-12)) continue;
        model(w, f);
        if (!w->usew) {
            /* T-form pass (identical to arbcheck) */
            /* sup |S_f| < c */
            if (mag_cmp(w->Shi[f], g_clo) < 0) return 1;
            /* a + sum |G_k| H_k + e < log c */
            if (model_aff(w, f)) {
                mag_set(w->m1, w->e[f]);
                for (k = 0; k < g_D; k++) {
                    arb_get_mag(w->m2, w->G[f][k]);
                    mag_addmul(w->m1, w->m2, w->Hm[k]);
                }
                arb_set(w->t1, w->a[f]);
                arb_add_error_mag(w->t1, w->m1);
                if (arb_lt(w->t1, g_logc)) return 1;
            }
        } else if (f && model_w(w, f)) {
            /* w-form pass: sup |S_f| < c, or aW + sum |GW_k| H_k + eW < log c */
            if (mag_cmp(w->ShiW[f], g_clo) < 0) return 1;
            if (model_waff(w, f)) {
                mag_set(w->m1, w->eW[f]);
                for (k = 0; k < g_D; k++) {
                    arb_get_mag(w->m2, w->GW[f][k]);
                    mag_addmul(w->m1, w->m2, w->Hm[k]);
                }
                arb_set(w->t1, w->aW[f]);
                arb_add_error_mag(w->t1, w->m1);
                if (arb_lt(w->t1, g_logc)) return 1;
            }
        }
        if (g_debug) {
            fprintf(stderr, "DBG F try f=%d est=%.4f fup=%.4f aff=%d up=%.4f e=%.3g logShi=%.4f logc=%.4f H0=%g\n", f, w->est[f], w->fup[f], w->aff[f],
                w->aff[f] ? arf_get_d(arb_midref(w->t1), ARF_RND_NEAR) + mag_get_d(arb_radref(w->t1)) : NAN,
                w->aff[f] ? mag_get_d(w->e[f]) : NAN, log(mag_get_d(w->Shi[f])), g_logc_d, w->b->H[0]);
        }
    }
    return 0;
}

static int try_pair(ws_t *w, int hi, int lo)
{
    int k, sh, sl;
    model(w, hi);
    model(w, lo);
    /* inf |S_hi| > sup |S_lo| or vice versa (in the w pass: the best of both forms for each function) */
    if (mag_cmp(Slo_of(w, hi), Shi_of(w, lo)) > 0) return 1;
    if (mag_cmp(Slo_of(w, lo), Shi_of(w, hi)) > 0) return 1;
    /* (a_hi - a_lo) +- (sum |G_hi,k - G_lo,k| H_k + e_hi + e_lo) excludes 0.
       T pass: both T-form models.  w pass: every combination that uses at least one w-form model. */
    for (sh = 0; sh < 2; sh++)
        for (sl = 0; sl < 2; sl++) {
            arb_srcptr ah, al, Gh, Gl;
            mag_srcptr eh, el;
            if (!w->usew && (sh || sl)) continue;
            if (w->usew && !sh && !sl) continue;
            if (!aff_get(w, hi, sh, &ah, &Gh, &eh) || !aff_get(w, lo, sl, &al, &Gl, &el)) continue;
            arb_sub(w->t1, ah, al, w->prec);
            mag_add(w->m1, eh, el);
            for (k = 0; k < g_D; k++) {
                arb_sub(w->t2, Gh + k, Gl + k, w->prec);
                arb_get_mag(w->m2, w->t2);
                mag_addmul(w->m1, w->m2, w->Hm[k]);
            }
            arb_add_error_mag(w->t1, w->m1);
            if (arb_is_positive(w->t1) || arb_is_negative(w->t1)) return 1;
        }
    return 0;
}

static int test_E(ws_t *w)
{
    int pi[MAXN * MAXN], pj[MAXN * MAXN];
    double sc[MAXN * MAXN];
    int np = 0, i, j, q;
    ensure_basic(w);
    ensure_est(w);
    for (i = 0; i < g_n; i++)
        for (j = i + 1; j < g_n; j++) {
            double s = -INFINITY, t;
            int k;
            if (isfinite(w->de[i]) && isfinite(w->de[j])) {
                s = fabs(w->est[i] - w->est[j]) - w->de[i] - w->de[j];
                for (k = 0; k < g_D; k++) s -= fabs(w->Gd[i][k] - w->Gd[j][k]) * w->b->H[k];
            }
            t = w->lSlo[i] - w->lShi[j]; if (t > s) s = t;
            t = w->lSlo[j] - w->lShi[i]; if (t > s) s = t;
            if (isnan(s)) s = -INFINITY;
            pi[np] = i; pj[np] = j; sc[np] = s;
            for (q = np; q > 0 && sc[q] > sc[q - 1]; q--) {
                double ts = sc[q]; sc[q] = sc[q - 1]; sc[q - 1] = ts;
                int ti = pi[q]; pi[q] = pi[q - 1]; pi[q - 1] = ti;
                ti = pj[q]; pj[q] = pj[q - 1]; pj[q - 1] = ti;
            }
            np++;
        }
    for (q = 0; q < np; q++)
        if (try_pair(w, pi[q], pj[q])) return 1;
    return 0;
}

static int test_L_w(ws_t *w);

static int test_L(ws_t *w)
{
    int f, k;
    double s = 0;
    ensure_basic(w);
    ensure_est(w);
    for (f = 0; f < g_n; f++) s += w->est[f];
    if (!(s < g_n * g_logc_d + 1e-12)) return 0;
    if (w->usew) return test_L_w(w);
    for (f = 0; f < g_n; f++) {
        model(w, f);
        if (mag_is_inf(w->Shi[f])) return 0;
    }
    /* pure magnitude bound: prod Shi < c^n */
    mag_one(w->m1);
    for (f = 0; f < g_n; f++) mag_mul(w->m1, w->m1, w->Shi[f]);
    if (mag_cmp(w->m1, g_cnlo) < 0) return 1;
    /* sum of affine models where available (gradients summed before taking |.|),
       magnitude bounds log Shi otherwise.  (model_aff uses the temporaries: call it first.) */
    for (f = 0; f < g_n; f++) model_aff(w, f);
    arb_zero(w->t1);
    mag_zero(w->m1);
    for (k = 0; k < g_D; k++) arb_zero(w->GS[k]);
    for (f = 0; f < g_n; f++) {
        if (w->aff[f]) {
            arb_add(w->t1, w->t1, w->a[f], w->prec);
            mag_add(w->m1, w->m1, w->e[f]);
            for (k = 0; k < g_D; k++) arb_add(w->GS[k], w->GS[k], w->G[f][k], w->prec);
        } else {
            if (mag_is_zero(w->Shi[f])) return 0; /* T identically 0: impossible; be conservative */
            arf_set_mag(arb_midref(w->t2), w->Shi[f]);
            mag_zero(arb_radref(w->t2));
            arb_log(w->t2, w->t2, w->prec);
            arb_add(w->t1, w->t1, w->t2, w->prec);
        }
    }
    for (k = 0; k < g_D; k++) {
        arb_get_mag(w->m2, w->GS[k]);
        mag_addmul(w->m1, w->m2, w->Hm[k]);
    }
    arb_add_error_mag(w->t1, w->m1);
    return arb_lt(w->t1, g_nlogc);
}

/*
 * L test, w pass. For each f one of: 0 = T-form affine model, 1 = w-form affine model,
 * 2 = log of the best magnitude bound (min of the T-form and w-form sup |S_f|).
 * The combination is chosen by a double-precision estimate (choice only); the chosen
 * combination is then decided in arb exactly as in the T pass.
 */
static int test_L_w(ws_t *w)
{
    int f, k, opt[MAXN], best[MAXN], av[MAXN][3], ncomb = 1, ci;
    double bestv = INFINITY;
    double ad[MAXN][3], ed[MAXN][3], Gdd[MAXN][2][MAXDIM];
    for (f = 0; f < g_n; f++) {
        mag_srcptr sh;
        model(w, f);
        sh = Shi_of(w, f);
        if (mag_is_inf(sh)) return 0;
        if (mag_is_zero(sh)) return 0;                /* S identically 0: impossible; be conservative */
    }
    /* pure magnitude bound with the best bounds: prod Shi < c^n */
    mag_one(w->m1);
    for (f = 0; f < g_n; f++) mag_mul(w->m1, w->m1, Shi_of(w, f));
    if (mag_cmp(w->m1, g_cnlo) < 0) return 1;
    for (f = 0; f < g_n; f++) {
        int s;
        for (s = 0; s < 3; s++) {
            arb_srcptr a, G;
            mag_srcptr e;
            av[f][s] = 0;
            if (s < 2) {
                if (!aff_get(w, f, s, &a, &G, &e)) continue;
                ad[f][s] = arf_get_d(arb_midref(a), ARF_RND_NEAR);
                ed[f][s] = mag_get_d(e) + mag_get_d(arb_radref(a));
                for (k = 0; k < g_D; k++) Gdd[f][s][k] = arf_get_d(arb_midref(G + k), ARF_RND_NEAR);
            } else {
                ad[f][s] = log(mag_get_d(Shi_of(w, f)));
                ed[f][s] = 0;
            }
            av[f][s] = isfinite(ad[f][s]) && isfinite(ed[f][s]);
        }
        if (!av[f][2] && !av[f][0] && !av[f][1]) return 0;
        ncomb *= 3;
    }
    /* enumerate the combinations (at most 3^6) in double precision */
    for (ci = 0; ci < ncomb; ci++) {
        int c = ci, okc = 1;
        double v = 0, GSd[MAXDIM] = {0};
        for (f = 0; f < g_n; f++) { opt[f] = c % 3; c /= 3; if (!av[f][opt[f]]) { okc = 0; break; } }
        if (!okc) continue;
        for (f = 0; f < g_n; f++) {
            v += ad[f][opt[f]] + ed[f][opt[f]];
            if (opt[f] < 2) for (k = 0; k < g_D; k++) GSd[k] += Gdd[f][opt[f]][k];
        }
        for (k = 0; k < g_D; k++) v += fabs(GSd[k]) * w->b->H[k];
        if (v < bestv) { bestv = v; memcpy(best, opt, sizeof(best)); }
    }
    if (!(bestv < INFINITY)) return 0;
    /* decide the chosen combination in arb */
    arb_zero(w->t1);
    mag_zero(w->m1);
    for (k = 0; k < g_D; k++) arb_zero(w->GS[k]);
    for (f = 0; f < g_n; f++) {
        if (best[f] < 2) {
            arb_srcptr a, G;
            mag_srcptr e;
            if (!aff_get(w, f, best[f], &a, &G, &e)) return 0;    /* cached; cannot fail here */
            arb_add(w->t1, w->t1, a, w->prec);
            mag_add(w->m1, w->m1, e);
            for (k = 0; k < g_D; k++) arb_add(w->GS[k], w->GS[k], G + k, w->prec);
        } else {
            arf_set_mag(arb_midref(w->t2), Shi_of(w, f));
            mag_zero(arb_radref(w->t2));
            arb_log(w->t2, w->t2, w->prec);
            arb_add(w->t1, w->t1, w->t2, w->prec);
        }
    }
    for (k = 0; k < g_D; k++) {
        arb_get_mag(w->m2, w->GS[k]);
        mag_addmul(w->m1, w->m2, w->Hm[k]);
    }
    arb_add_error_mag(w->t1, w->m1);
    return arb_lt(w->t1, g_nlogc);
}

/* lower bound of |x| on [C-H, C+H] (nonnegative ball) */
static void lo_abs(arb_t r, arb_t tmp, double C, double H, slong prec)
{
    arb_set_d(r, fabs(C));
    arb_set_d(tmp, H);
    arb_sub(r, r, tmp, prec);
    if (!arb_is_positive(r)) arb_zero(r);
}

static void hi_abs(arb_t r, arb_t tmp, double C, double H, slong prec)
{
    arb_set_d(r, fabs(C));
    arb_set_d(tmp, H);
    arb_add(r, r, tmp, prec);
}

/* inf over the box of |u_v|^2 */
static void inf_abs2(ws_t *w, arb_t r, int v)
{
    slong prec = w->prec;
    lo_abs(w->t3, w->t4, w->b->C[2 * v], w->b->H[2 * v], prec);
    arb_sqr(r, w->t3, prec);
    lo_abs(w->t3, w->t4, w->b->C[2 * v + 1], w->b->H[2 * v + 1], prec);
    arb_addmul(r, w->t3, w->t3, prec);
}

static void sup_abs2(ws_t *w, arb_t r, int v)
{
    slong prec = w->prec;
    hi_abs(w->t3, w->t4, w->b->C[2 * v], w->b->H[2 * v], prec);
    arb_sqr(r, w->t3, prec);
    hi_abs(w->t3, w->t4, w->b->C[2 * v + 1], w->b->H[2 * v + 1], prec);
    arb_addmul(r, w->t3, w->t3, prec);
}

static int test_O(ws_t *w)
{
    int v;
    for (v = 0; v < g_nv; v++) {
        inf_abs2(w, w->t1, v);
        if (arb_gt(w->t1, g_one)) return 1;
    }
    return 0;
}

static int test_S(ws_t *w)
{
    int v;
    slong prec = w->prec;
    /* Im u_2 < 0 on the box: C_1 + H_1 < 0 */
    arb_set_d(w->t1, w->b->C[1]);
    arb_set_d(w->t2, w->b->H[1]);
    arb_add(w->t1, w->t1, w->t2, prec);
    if (arb_is_negative(w->t1)) return 1;
    for (v = 0; v + 1 < g_nv; v++) {
        sup_abs2(w, w->t1, v);
        inf_abs2(w, w->t2, v + 1);
        if (arb_lt(w->t1, w->t2)) return 1;
    }
    return 0;
}

static int test_X(ws_t *w)
{
    int p, k;
    slong prec = w->prec;
    const box_t *b = w->b;
    for (p = 0; p < g_npts; p++) {
        /* quick non-rigorous filter: centre distance must be < 1/20 */
        double d2 = 0;
        for (k = 0; k < g_D; k++) {
            double t = b->C[k] - g_peq_d[p * g_D + k];
            d2 += t * t;
        }
        if (d2 > 0.0026) continue;
        arb_zero(w->t1);
        for (k = 0; k < g_D; k++) {
            arb_set_d(w->t2, b->C[k]);
            arb_sub(w->t2, w->t2, g_peq + p * g_D + k, prec);
            arb_abs(w->t2, w->t2);
            arb_set_d(w->t3, b->H[k]);
            arb_add(w->t2, w->t2, w->t3, prec);
            arb_addmul(w->t1, w->t2, w->t2, prec);
        }
        if (arb_le(w->t1, g_r2excl)) return 1;
    }
    return 0;
}

static int run_test(ws_t *w, int code)
{
    switch (code) {
        case C_F: return test_F(w);
        case C_E: return test_E(w);
        case C_L: return test_L(w);
        case C_O: return test_O(w);
        case C_S: return test_S(w);
        case C_X: return test_X(w);
        default: return 0;
    }
}

/*
 * Prove the claim `code` on box b. First the T-form pass (exactly the tests of arbcheck);
 * if that fails and the claim is F, E or L, the w-form pass.
 * Returns 0 (not proved), 1 (proved by the T-form pass), 2 (proved only in the w-form pass).
 */
static int claim(ws_t *w, const box_t *b, int code)
{
    int r;
    ws_set_box(w, b);
    w->usew = 0;
    if (run_test(w, code)) return 1;
    if (!g_wform || (code != C_F && code != C_E && code != C_L)) return 0;
    w->usew = 1;
    r = run_test(w, code);
    w->usew = 0;
    return r ? 2 : 0;
}

/* any claim of the table, the claimed one first; T-form pass for all claims before the w-form pass */
static int any_claim(ws_t *w, const box_t *b, int code)
{
    static const int order[6] = {C_O, C_S, C_X, C_F, C_E, C_L};
    int i, r;
    if ((r = claim(w, b, code))) return r;     /* claim() attached b; caches stay valid */
    w->usew = 0;
    for (i = 0; i < 6; i++) {
        int c = order[i];
        if (c == code) continue;
        if (run_test(w, c)) return 1;
    }
    if (!g_wform) return 0;
    w->usew = 1;
    for (i = 3; i < 6; i++) {                  /* F, E, L */
        int c = order[i];
        if (c == code) continue;
        if (run_test(w, c)) { w->usew = 0; return 2; }
    }
    w->usew = 0;
    return 0;
}

/* exact halving of coordinate k; returns 0 if not exact (never expected) */
static int exact_split(double C, double H, double *lo, double *hi, double *h2)
{
    double h = H * 0.5, s, bb, err;
    if (!(h >= 2.2250738585072014e-308) || h * 2.0 != H) return 0;
    s = C - h; bb = s - C; err = (C - (s - bb)) + (-h - bb);
    if (err != 0) return 0;
    *lo = s;
    s = C + h; bb = s - C; err = (C - (s - bb)) + (h - bb);
    if (err != 0) return 0;
    *hi = s;
    *h2 = h;
    return 1;
}

static int refine(ws_t *w, const box_t *b, int code, int depth)
{
    box_t c;
    int k, best = 0, side, r;
    double lo, hi, h2;
    for (k = 1; k < g_D; k++) if (b->H[k] > b->H[best]) best = k;
    if (!exact_split(b->C[best], b->H[best], &lo, &hi, &h2)) return 0;
    for (side = 0; side < 2; side++) {
        c = *b;
        c.C[best] = side ? hi : lo;
        c.H[best] = h2;
        w->subboxes++;
        if ((r = any_claim(w, &c, code))) {
            if (r == 2) w->wsub++;
            continue;
        }
        if (depth < g_maxref && refine(w, &c, code, depth + 1)) continue;
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* tasks and tree file                                                  */

typedef struct {
    uint32_t id;
    uint64_t off;      /* offset of node data */
    uint64_t nbytes;
} rec_t;

static uint32_t g_ntask = 0, g_maxid = 0;
static double *g_centres = NULL;     /* (g_maxid+1) * D */
static unsigned char *g_hastask = NULL;
static double g_hw = 0;

static int read_tasks(const char *path)
{
    FILE *fp = fopen(path, "r");
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    uint32_t cap_ids = 0;
    int lineno = 0;
    if (!fp) { perror(path); return 0; }
    while ((len = getline(&line, &cap, fp)) > 0) {
        lineno++;
        if (lineno == 1) {
            char *p = strstr(line, "half_width=");
            if (!p) { fprintf(stderr, "tasks: no half_width in header\n"); return 0; }
            g_hw = strtod(p + 11, NULL);
            if (!(g_hw > 0) || !isfinite(g_hw)) { fprintf(stderr, "tasks: bad half_width\n"); return 0; }
            continue;
        }
        {
            char *p = line, *q;
            unsigned long id;
            int k;
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\n' || *p == 0) continue;
            id = strtoul(p, &q, 10);
            if (q == p || id > 0xffffffffUL) { fprintf(stderr, "tasks: bad line %d\n", lineno); return 0; }
            if (id >= cap_ids) {
                uint32_t nc = cap_ids ? cap_ids : 1024;
                while (nc <= id) nc *= 2;
                g_centres = realloc(g_centres, (size_t) nc * g_D * sizeof(double));
                g_hastask = realloc(g_hastask, nc);
                memset(g_hastask + cap_ids, 0, nc - cap_ids);
                cap_ids = nc;
            }
            if (g_hastask[id]) { fprintf(stderr, "tasks: duplicate id %lu\n", id); return 0; }
            p = q;
            for (k = 0; k < g_D; k++) {
                double x = strtod(p, &q);
                if (q == p) { fprintf(stderr, "tasks: line %d has too few coordinates\n", lineno); return 0; }
                if (!isfinite(x)) { fprintf(stderr, "tasks: line %d has a non-finite coordinate\n", lineno); return 0; }
                g_centres[(size_t) id * g_D + k] = x;
                p = q;
            }
            while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
            if (*p) { fprintf(stderr, "tasks: line %d has trailing data\n", lineno); return 0; }
            g_hastask[id] = 1;
            g_ntask++;
            if (id > g_maxid) g_maxid = (uint32_t) id;
        }
    }
    free(line);
    fclose(fp);
    return 1;
}

/* per-task statistics */
typedef struct {
    uint64_t nodes, leaves[NCODE], direct, refined, subboxes, fail;
    uint64_t wdirect[NCODE], wsub;   /* leaves proved directly only by the w-form pass; refinement sub-boxes likewise */
    int maxdepth;
    int parse_err;             /* 0 ok, 1 missing bytes, 2 leftover bytes, 3 invalid code, 4 inexact split, 5 too deep */
} tstat_t;

static const char *PERR[] = {"ok", "missing_bytes", "leftover_bytes", "invalid_code", "inexact_split", "too_deep"};

typedef struct {
    FILE *fp;
    unsigned char *buf;
    size_t pos, len, bufsize;
    uint64_t remaining;
} reader_t;

static inline int rd_byte(reader_t *r)
{
    if (r->pos == r->len) {
        size_t want, got;
        if (r->remaining == 0) return -1;
        want = r->remaining < r->bufsize ? (size_t) r->remaining : r->bufsize;
        got = fread(r->buf, 1, want, r->fp);
        if (got == 0) return -1;
        r->remaining -= got;
        r->pos = 0;
        r->len = got;
    }
    return r->buf[r->pos++];
}

static pthread_mutex_t g_out_lock = PTHREAD_MUTEX_INITIALIZER;
static int g_failprinted = 0;

typedef struct {
    ws_t *w;
    reader_t rd;
    tstat_t *st;
    uint32_t id;
    uint64_t tns[NCODE], tcnt[NCODE], tmod[NCODE];
} ctx_t;

static uint64_t now_ns(void)
{
#ifdef __APPLE__
    return clock_gettime_nsec_np(CLOCK_UPTIME_RAW);
#else
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t) t.tv_sec * 1000000000u + (uint64_t) t.tv_nsec;
#endif
}

static int g_refprinted = 0;
static void report_fail(ctx_t *cx, const box_t *b, int code, const char *why)
{
    int k, isnote = why[0] == 'n';
    int *cnt = isnote ? &g_refprinted : &g_failprinted;
    pthread_mutex_lock(&g_out_lock);
    if (*cnt < g_maxfailprint) {
        (*cnt)++;
        fprintf(stderr, "%s id=%u code=%s (%s) box:", isnote ? "REFINED" : "FAIL", cx->id, code >= 0 && code < NCODE ? CODE_NAME[code] : "?", why);
        for (k = 0; k < g_D; k++) fprintf(stderr, " %a+-%a", b->C[k], b->H[k]);
        fprintf(stderr, "\n");
    }
    pthread_mutex_unlock(&g_out_lock);
}

static void verify_leaf(ctx_t *cx, const box_t *b, int code)
{
    tstat_t *st = cx->st;
    uint64_t t0 = g_timing ? now_ns() : 0, nm0 = cx->w->nmodel;
    int r;
    st->leaves[code]++;
    if (code == C_U) {
        st->fail++;
        report_fail(cx, b, code, "unresolved leaf");
    } else if ((r = claim(cx->w, b, code))) {
        st->direct++;
        if (r == 2) st->wdirect[code]++;
    } else {
        uint64_t sb = cx->w->subboxes, ws0 = cx->w->wsub;
        int ok = g_maxref > 0 && refine(cx->w, b, code, 1);
        st->subboxes += cx->w->subboxes - sb;
        st->wsub += cx->w->wsub - ws0;
        if (ok) {
            st->refined++;
            report_fail(cx, b, code, "note: proved only after refinement");
        }
        else {
            st->fail++;
            report_fail(cx, b, code, g_maxref > 0 ? "claim not proved, refinement failed" : "claim not proved, refinement disabled");
        }
    }
    if (g_timing) {
        cx->tns[code] += now_ns() - t0;
        cx->tcnt[code]++;
        cx->tmod[code] += cx->w->nmodel - nm0;
    }
}

static int walk(ctx_t *cx, box_t *b, int depth)
{
    int c = rd_byte(&cx->rd), k;
    double C, H, lo, hi, h2;
    tstat_t *st = cx->st;
    if (c < 0) { st->parse_err = 1; return -1; }
    st->nodes++;
    if (depth > st->maxdepth) st->maxdepth = depth;
    if (c >= 16) {
        k = c - 16;
        if (k >= g_D) { st->parse_err = 3; return -1; }
        if (depth >= MAXDEPTH_TREE) { st->parse_err = 5; return -1; }
        C = b->C[k]; H = b->H[k];
        if (!exact_split(C, H, &lo, &hi, &h2)) { st->parse_err = 4; return -1; }
        b->C[k] = lo; b->H[k] = h2;
        if (walk(cx, b, depth + 1) < 0) return -1;
        b->C[k] = hi; b->H[k] = h2;
        if (walk(cx, b, depth + 1) < 0) return -1;
        b->C[k] = C; b->H[k] = H;
        return 0;
    }
    if (c >= NCODE) { st->parse_err = 3; return -1; }
    verify_leaf(cx, b, c);
    return 0;
}

/* ------------------------------------------------------------------ */
/* driver                                                               */

static const char *g_tree_path;
static rec_t *g_jobs;
static size_t g_njobs;
static size_t g_next = 0;
static tstat_t g_tot;
static uint64_t g_tns[NCODE], g_tcnt[NCODE], g_tmod[NCODE];
static uint64_t g_fail_records = 0;

static void *worker(void *arg)
{
    ctx_t cx;
    (void) arg;
    memset(&cx, 0, sizeof(cx));
    cx.w = ws_new(g_prec);
    cx.rd.fp = fopen(g_tree_path, "rb");
    cx.rd.bufsize = 1 << 20;
    cx.rd.buf = malloc(cx.rd.bufsize);
    if (!cx.rd.fp) { perror(g_tree_path); exit(2); }
    for (;;) {
        size_t j = __atomic_fetch_add(&g_next, 1, __ATOMIC_RELAXED);
        rec_t *rc;
        tstat_t st;
        box_t b;
        int k, rfail;
        if (j >= g_njobs) break;
        rc = &g_jobs[j];
        memset(&st, 0, sizeof(st));
        cx.st = &st;
        cx.id = rc->id;
        for (k = 0; k < g_D; k++) { b.C[k] = g_centres[(size_t) rc->id * g_D + k]; b.H[k] = g_hw; }
        if (fseeko(cx.rd.fp, (off_t) rc->off, SEEK_SET) != 0) { perror("fseeko"); exit(2); }
        cx.rd.pos = cx.rd.len = 0;
        cx.rd.remaining = rc->nbytes;
        if (walk(&cx, &b, 0) == 0) {
            if (cx.rd.pos != cx.rd.len || cx.rd.remaining != 0) st.parse_err = 2;
        }
        rfail = st.fail > 0 || st.parse_err != 0;
        pthread_mutex_lock(&g_out_lock);
        if (!g_quiet || rfail)
            printf("task %u nodes=%" PRIu64 " maxdepth=%d F=%" PRIu64 " E=%" PRIu64 " L=%" PRIu64 " O=%" PRIu64
                   " S=%" PRIu64 " X=%" PRIu64 " U=%" PRIu64 " direct=%" PRIu64 " refined=%" PRIu64
                   " subboxes=%" PRIu64 " fail=%" PRIu64 " wdirect=%" PRIu64 " wsub=%" PRIu64 " parse=%s\n",
                   rc->id, st.nodes, st.maxdepth, st.leaves[0], st.leaves[1], st.leaves[2], st.leaves[3],
                   st.leaves[4], st.leaves[5], st.leaves[6], st.direct, st.refined, st.subboxes,
                   st.fail + (st.parse_err ? 1 : 0), st.wdirect[0] + st.wdirect[1] + st.wdirect[2], st.wsub, PERR[st.parse_err]);
        fflush(stdout);
        g_tot.nodes += st.nodes;
        for (k = 0; k < NCODE; k++) g_tot.leaves[k] += st.leaves[k];
        g_tot.direct += st.direct; g_tot.refined += st.refined; g_tot.subboxes += st.subboxes;
        for (k = 0; k < NCODE; k++) g_tot.wdirect[k] += st.wdirect[k];
        g_tot.wsub += st.wsub;
        g_tot.fail += st.fail + (st.parse_err ? 1 : 0);
        if (st.maxdepth > g_tot.maxdepth) g_tot.maxdepth = st.maxdepth;
        if (rfail) g_fail_records++;
        pthread_mutex_unlock(&g_out_lock);
    }
    pthread_mutex_lock(&g_out_lock);
    for (int k = 0; k < NCODE; k++) { g_tns[k] += cx.tns[k]; g_tcnt[k] += cx.tcnt[k]; g_tmod[k] += cx.tmod[k]; }
    pthread_mutex_unlock(&g_out_lock);
    fclose(cx.rd.fp);
    free(cx.rd.buf);
    ws_free(cx.w);
    flint_cleanup();
    return NULL;
}

static void init_constants(void)
{
    int k, i;
    slong prec = g_prec + 32;
    for (k = 0; k <= MAXN; k++) {
        arb_init(g_smu[k]);
        arb_one(g_smu[k]);
        arb_div_ui(g_smu[k], g_smu[k], (ulong) (k + 1) * (k + 2), prec);
        if (k & 1) arb_neg(g_smu[k], g_smu[k]);
        g_smu_d[k] = ((k & 1) ? -1.0 : 1.0) / ((k + 1) * (k + 2));
    }
    for (k = 0; k < MAXN; k++)
        for (i = 0; i < MAXN; i++) arb_init(g_K[k][i]);
    for (k = 0; k < g_n; k++)            /* k = a */
        for (i = k; i < g_n; i++) {      /* i = p */
            fmpz_t bin;
            fmpz_init(bin);
            fmpz_bin_uiui(bin, i, k);
            arb_mul_fmpz(g_K[k][i], g_smu[g_n - 1 - i], bin, prec);
            fmpz_clear(bin);
        }
    /* w-form: g_KW[a][p] = C(n-1-p, a) (-1)^{n-1-p} mu_{n-1-p}, p = 0..n-1-a */
    for (k = 0; k < MAXN; k++)
        for (i = 0; i < MAXN; i++) arb_init(g_KW[k][i]);
    for (k = 0; k < g_n; k++)            /* k = a */
        for (i = 0; i <= g_n - 1 - k; i++) {   /* i = p */
            fmpz_t bin;
            fmpz_init(bin);
            fmpz_bin_uiui(bin, g_n - 1 - i, k);
            arb_mul_fmpz(g_KW[k][i], g_smu[g_n - 1 - i], bin, prec);
            fmpz_clear(bin);
        }
    arb_init(g_logc); arb_init(g_nlogc); arb_init(g_r2excl); arb_init(g_one);
    arb_one(g_one);
    arb_set_ui(g_logc, g_d - 1);
    arb_div_ui(g_logc, g_logc, g_d, prec);
    arb_log(g_logc, g_logc, prec);
    arb_mul_ui(g_nlogc, g_logc, g_n, prec);
    g_logc_d = log((double) (g_d - 1) / g_d);
    arb_one(g_r2excl);
    arb_div_ui(g_r2excl, g_r2excl, 400, prec);
    mag_init(g_clo); mag_init(g_cnlo);
    {
        mag_t a, b;
        mag_init(a); mag_init(b);
        mag_set_ui_lower(a, g_d - 1);
        mag_set_ui(b, g_d);
        mag_div_lower(g_clo, a, b);
        mag_pow_ui_lower(g_cnlo, g_clo, g_n);
        mag_clear(a); mag_clear(b);
    }
    /* equality points: u_j = conj(omega^{a_j}), (a_2..a_n) a permutation of (1..n-1) */
    {
        int perm[MAXNV], np = 1, c[MAXNV], idx;
        arb_t s, co;
        fmpq_t x;
        for (i = 2; i <= g_nv; i++) np *= i;
        g_npts = np;
        g_peq = _arb_vec_init((slong) np * g_D);
        g_peq_d = malloc(sizeof(double) * np * g_D);
        arb_init(s); arb_init(co); fmpq_init(x);
        for (i = 0; i < g_nv; i++) { perm[i] = i + 1; c[i] = 0; }
        /* Heap's algorithm */
        idx = 0;
        for (;;) {
            for (i = 0; i < g_nv; i++) {
                fmpq_set_si(x, 2 * perm[i], g_n);
                arb_sin_cos_pi_fmpq(s, co, x, prec);
                arb_set(g_peq + idx * g_D + 2 * i, co);
                arb_neg(g_peq + idx * g_D + 2 * i + 1, s);
                g_peq_d[idx * g_D + 2 * i] = arf_get_d(arb_midref(co), ARF_RND_NEAR);
                g_peq_d[idx * g_D + 2 * i + 1] = -arf_get_d(arb_midref(s), ARF_RND_NEAR);
            }
            idx++;
            i = 1;
            while (i < g_nv && c[i] >= i) { c[i] = 0; i++; }
            if (i >= g_nv) break;
            if (i % 2 == 0) { int t = perm[0]; perm[0] = perm[i]; perm[i] = t; }
            else { int t = perm[c[i]]; perm[c[i]] = perm[i]; perm[i] = t; }
            c[i]++;
        }
        if (idx != np) { fprintf(stderr, "internal: permutation count %d != %d\n", idx, np); exit(2); }
        arb_clear(s); arb_clear(co); fmpq_clear(x);
    }
}

/* ------------------------------------------------------------------ */
/* self test: compare the expansions with a direct double evaluation   */

static double complex direct_S(int f, const double complex *u)
{
    /* integrand polynomial in t, integrate (1-t) t^k exactly */
    double complex p[MAXN + 2];
    int deg = 0, j, k;
    double complex res = 0;
    p[0] = 1;
    if (f == 0) {
        for (j = 0; j < g_nv; j++) {
            /* multiply by (1 - t u_j) */
            p[deg + 1] = 0;
            for (k = deg + 1; k >= 1; k--) p[k] = p[k] - u[j] * p[k - 1];
            deg++;
        }
    } else {
        int vi = f - 1;
        /* (u_i - t) */
        p[1] = -1; p[0] = u[vi]; deg = 1;
        for (j = 0; j < g_nv; j++) {
            if (j == vi) continue;
            p[deg + 1] = 0;
            for (k = deg + 1; k >= 1; k--) p[k] = u[vi] * p[k] - u[j] * p[k - 1];
            p[0] = u[vi] * p[0];
            deg++;
        }
    }
    for (k = 0; k <= deg; k++) res += p[k] / ((k + 1.0) * (k + 2.0));
    if (f) res /= cpow(u[f - 1], g_n - 1);
    return res;
}

/* direct double evaluation of the w-form integral of SPEC_W.md: int (1-t)(1-tw) prod_{j != i} (1 - t u_j w) dt, w = 1/u_i */
static double complex direct_Sw(int f, const double complex *u)
{
    double complex p[MAXN + 2], wv, res = 0;
    int deg = 1, j, k, vi = f - 1;
    wv = 1.0 / u[vi];
    p[0] = 1; p[1] = -wv;
    for (j = 0; j < g_nv; j++) {
        if (j == vi) continue;
        p[deg + 1] = 0;
        for (k = deg + 1; k >= 1; k--) p[k] = p[k] - u[j] * wv * p[k - 1];
        deg++;
    }
    for (k = 0; k <= deg; k++) res += p[k] / ((k + 1.0) * (k + 2.0));
    return res;
}

/*
 * Ball evaluation at an exact point (doubles) of both sides of the w-form identity, from the definitions:
 *   lhs = T_i(u) / u_i^{n-1},  T_i = int (1-t)(u_i - t) prod_{j != i} (u_i - t u_j) dt
 *   rhs = int (1-t)(1 - t w) prod_{j != i} (1 - t u_j w) dt,  w = 1/u_i
 * Returns 1 if the two balls overlap (they must, as they enclose the same number).
 */
static int identity_check(int f, const double *x, slong prec, double *reldiff)
{
    acb_ptr p = _acb_vec_init(MAXN + 2), u = _acb_vec_init(MAXNV);
    acb_t lhs, rhs, wv, t, mu;
    int deg, j, k, vi = f - 1, ok;
    acb_init(lhs); acb_init(rhs); acb_init(wv); acb_init(t); acb_init(mu);
    for (j = 0; j < g_nv; j++) acb_set_d_d(u + j, x[2 * j], x[2 * j + 1]);
    /* T_i */
    _acb_vec_zero(p, MAXN + 2);
    acb_set(p + 0, u + vi); acb_set_si(p + 1, -1); deg = 1;
    for (j = 0; j < g_nv; j++) {
        if (j == vi) continue;
        for (k = deg + 1; k >= 1; k--) {           /* p <- p (u_i - t u_j) */
            acb_mul(p + k, p + k, u + vi, prec);
            acb_submul(p + k, p + k - 1, u + j, prec);
        }
        acb_mul(p + 0, p + 0, u + vi, prec);
        deg++;
    }
    acb_zero(lhs);
    for (k = 0; k <= deg; k++) { acb_set_ui(mu, (ulong) (k + 1) * (k + 2)); acb_div(t, p + k, mu, prec); acb_add(lhs, lhs, t, prec); }
    acb_pow_ui(t, u + vi, g_n - 1, prec);
    acb_div(lhs, lhs, t, prec);
    /* w-form */
    acb_inv(wv, u + vi, prec);
    _acb_vec_zero(p, MAXN + 2);
    acb_one(p + 0); acb_neg(p + 1, wv); deg = 1;
    for (j = 0; j < g_nv; j++) {
        if (j == vi) continue;
        acb_mul(t, u + j, wv, prec);
        for (k = deg + 1; k >= 1; k--) acb_submul(p + k, p + k - 1, t, prec);   /* p <- p (1 - t u_j w) */
        deg++;
    }
    acb_zero(rhs);
    for (k = 0; k <= deg; k++) { acb_set_ui(mu, (ulong) (k + 1) * (k + 2)); acb_div(t, p + k, mu, prec); acb_add(rhs, rhs, t, prec); }
    ok = acb_overlaps(lhs, rhs) && acb_is_finite(lhs) && acb_is_finite(rhs);
    acb_sub(t, lhs, rhs, prec);
    {
        double dn = mag_get_d(arb_radref(acb_realref(t))) + mag_get_d(arb_radref(acb_imagref(t)))
                    + hypot(arf_get_d(arb_midref(acb_realref(t)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_NEAR));
        double dd = hypot(arf_get_d(arb_midref(acb_realref(lhs)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(lhs)), ARF_RND_NEAR));
        *reldiff = dd > 0 ? dn / dd : dn;
    }
    acb_clear(lhs); acb_clear(rhs); acb_clear(wv); acb_clear(t); acb_clear(mu);
    _acb_vec_clear(p, MAXN + 2); _acb_vec_clear(u, MAXNV);
    return ok;
}

static int selftest(void)
{
    ws_t *w = ws_new(g_prec);
    flint_rand_t st;
    int trial, f, bad = 0, naff = 0, k;
    int nwm = 0, nwaff = 0, nwdomfail = 0, nwpts = 0, nid = 0, ntarget = 0;
    double maxratio = 0, maxratio_w = 0, maxid = 0, maxid_d = 0, maxc0 = 0;
    flint_rand_init(st);
    for (trial = 0; trial < 4000; trial++) {
        box_t b;
        double complex u[MAXNV];
        double scale = ldexp(1.0, -(int) (n_randint(st, 12) + 2));
        for (k = 0; k < g_D; k++) {
            b.C[k] = ldexp((double) ((slong) n_randint(st, 2048) - 1024), -10);
            b.H[k] = scale;
        }
        if (trial >= 2000) {
            /* boxes near the edge of the w-form domain: |m_vi| = r rho_vi with r in (1, 3) */
            int vi = (int) n_randint(st, g_nv);
            double rho = sqrt(2.0) * scale, r = 1.0 + 2.0 * (double) (n_randint(st, 1000000) + 1) / 1000001.0;
            double th = 2 * M_PI * (double) n_randint(st, 1000000) / 1000000.0;
            if (trial >= 3000) r = 1.0 + ldexp((double) (n_randint(st, 1000) + 1), -14);   /* very close to the edge */
            b.C[2 * vi] = r * rho * cos(th);
            b.C[2 * vi + 1] = r * rho * sin(th);
            ntarget++;
        }
        ws_set_box(w, &b);
        ensure_basic(w);
        ensure_est(w);
        for (f = 0; f < g_n; f++) {
            int s, hasw = 0, hwaff = 0;
            model(w, f);
            if (model_aff(w, f)) naff++;
            if (f) {
                w->usew = 1;
                hasw = model_w(w, f);
                hwaff = hasw && model_waff(w, f);
                w->usew = 0;
                if (hasw) nwm++; else nwdomfail++;
                if (hwaff) nwaff++;
                if (hasw) {
                    /* the w-form constant term is S_f at the centre */
                    double complex mc[MAXNV], Sm, c0;
                    for (k = 0; k < g_nv; k++) mc[k] = b.C[2 * k] + b.C[2 * k + 1] * _Complex_I;
                    Sm = direct_S(f, mc);
                    c0 = arf_get_d(arb_midref(acb_realref(w->cw[f][0])), ARF_RND_NEAR)
                         + arf_get_d(arb_midref(acb_imagref(w->cw[f][0])), ARF_RND_NEAR) * _Complex_I;
                    if (cabs(c0 - Sm) / (1e-300 + cabs(Sm)) > maxc0) maxc0 = cabs(c0 - Sm) / (1e-300 + cabs(Sm));
                    if (cabs(c0 - Sm) > 1e-9 * (1 + cabs(Sm))) { printf("w-form c0 mismatch f=%d %g %g\n", f, cabs(c0), cabs(Sm)); bad++; }
                }
            }
            for (s = 0; s < 24; s++) {
                double x[MAXDIM], val, lin;
                for (k = 0; k < g_D; k++) {
                    double r = (double) n_randint(st, 1000001) / 500000.0 - 1.0;
                    if (s == 0) r = 0;
                    if (s >= 1 && s <= 4) r = n_randint(st, 2) ? 1.0 : -1.0;   /* box corners */
                    x[k] = b.C[k] + r * b.H[k];
                }
                for (k = 0; k < g_nv; k++) u[k] = x[2 * k] + x[2 * k + 1] * _Complex_I;
                val = log(cabs(direct_S(f, u)));
                if (s == 0 && fabs(val - w->est[f]) > 1e-9 * (1 + fabs(val))) {
                    printf("est mismatch f=%d %g %g\n", f, val, w->est[f]); bad++;
                }
                /* the w-form identity S_i = T_i / u_i^{n-1}, at points with u_i != 0 */
                if (f && cabs(u[f - 1]) > 0 && s < 6) {
                    double rd, dd;
                    double complex a1 = direct_S(f, u), a2 = direct_Sw(f, u);
                    nid++;
                    if (!identity_check(f, x, 128, &rd)) { printf("w-form identity fails (balls) f=%d\n", f); bad++; }
                    if (rd > maxid) maxid = rd;
                    dd = cabs(a1 - a2) / (1e-300 + cabs(a1));
                    if (dd > maxid_d) maxid_d = dd;
                    if (rd > 1e-20) { printf("w-form identity: relative difference %g f=%d\n", rd, f); bad++; }
                }
                if (!isfinite(val)) continue;
                if (mag_get_d(w->Shi[f]) * (1 + 1e-12) < exp(val) || mag_get_d(w->Slo[f]) > exp(val) * (1 + 1e-12)) {
                    printf("magnitude bound violated f=%d |S|=%g [%g,%g]\n", f, exp(val), mag_get_d(w->Slo[f]), mag_get_d(w->Shi[f]));
                    bad++;
                }
                if (w->aff[f]) {
                    lin = arf_get_d(arb_midref(w->a[f]), ARF_RND_NEAR);
                    for (k = 0; k < g_D; k++) lin += arf_get_d(arb_midref(w->G[f][k]), ARF_RND_NEAR) * (x[k] - b.C[k]);
                    if (mag_get_d(w->e[f]) > 0 && fabs(val - lin) / mag_get_d(w->e[f]) > maxratio) maxratio = fabs(val - lin) / mag_get_d(w->e[f]);
                    if (fabs(val - lin) > mag_get_d(w->e[f]) * (1 + 1e-9) + 1e-12) {
                        printf("affine bound violated f=%d val=%.17g lin=%.17g e=%g\n", f, val, lin, mag_get_d(w->e[f]));
                        bad++;
                    }
                }
                /* w-form interval bounds and affine model */
                if (hasw) {
                    nwpts++;
                    if (mag_get_d(w->ShiW[f]) * (1 + 1e-12) < exp(val) || mag_get_d(w->SloW[f]) > exp(val) * (1 + 1e-12)) {
                        printf("w-form magnitude bound violated f=%d |S|=%g [%g,%g]\n", f, exp(val), mag_get_d(w->SloW[f]), mag_get_d(w->ShiW[f]));
                        bad++;
                    }
                }
                if (hwaff) {
                    double ew = mag_get_d(w->eW[f]);
                    lin = arf_get_d(arb_midref(w->aW[f]), ARF_RND_NEAR);
                    for (k = 0; k < g_D; k++) lin += arf_get_d(arb_midref(w->GW[f][k]), ARF_RND_NEAR) * (x[k] - b.C[k]);
                    if (ew > 0 && fabs(val - lin) / ew > maxratio_w) maxratio_w = fabs(val - lin) / ew;
                    if (fabs(val - lin) > ew * (1 + 1e-9) + 1e-12) {
                        printf("w-form affine bound violated f=%d val=%.17g lin=%.17g e=%g\n", f, val, lin, ew);
                        bad++;
                    }
                }
            }
        }
    }
    /* every equality point must have |S_f| = c for all f, and the points must be distinct */
    {
        int p, q2, v;
        double maxdev = 0;
        for (p = 0; p < g_npts; p++) {
            double complex u[MAXNV];
            for (v = 0; v < g_nv; v++) u[v] = g_peq_d[p * g_D + 2 * v] + g_peq_d[p * g_D + 2 * v + 1] * _Complex_I;
            for (f = 0; f < g_n; f++) {
                double dev = fabs(cabs(direct_S(f, u)) - (double) (g_d - 1) / g_d);
                if (dev > maxdev) maxdev = dev;
            }
            for (q2 = 0; q2 < p; q2++) {
                double dd = 0;
                for (k = 0; k < g_D; k++) dd += fabs(g_peq_d[p * g_D + k] - g_peq_d[q2 * g_D + k]);
                if (dd < 1e-9) { printf("equality points %d and %d coincide\n", p, q2); bad++; }
            }
        }
        printf("equality points: %d, max | |S_f(p)| - c | = %.2e\n", g_npts, maxdev);
        if (maxdev > 1e-12) bad++;
    }
    printf("w-form identity S_i = T_i/u_i^(n-1): %d points, balls (128 bit) overlap at all; max relative difference %.2e (arb), %.2e (double)\n",
           nid, maxid, maxid_d);
    printf("w-form models: %d checked (%d boxes near the domain edge among 4000), %d not applicable (box meets |u_i| <= rho_i); "
           "%d affine; %d point checks; max |c0 - S(m)|/|S(m)| = %.2e; max |log|S|-affine|/e = %.3f\n",
           nwm, ntarget, nwdomfail, nwaff, nwpts, maxc0, maxratio_w);
    printf("selftest d=%d: %d problems, %d affine models checked, max |log|S|-affine|/e = %.3f\n", g_d, bad, naff, maxratio);
    flint_rand_clear(st);
    ws_free(w);
    return bad == 0;
}

static void bench(void)
{
    ws_t *w = ws_new(g_prec);
    flint_rand_t st;
    int N = 20000, i, k, f;
    box_t *bx = malloc(sizeof(box_t) * N);
    uint64_t t0, tes = 0, test_ = 0, tmod[MAXN] = {0}, tF = 0;
    flint_rand_init(st);
    for (i = 0; i < N; i++)
        for (k = 0; k < g_D; k++) { bx[i].C[k] = ldexp((double) ((slong) n_randint(st, 64) - 32), -5) + ldexp(1, -6); bx[i].H[k] = ldexp(1, -6); }
    for (i = 0; i < N; i++) { ws_set_box(w, &bx[i]); ensure_basic(w); t0 = now_ns(); ensure_est(w); test_ += now_ns() - t0; }
    for (i = 0; i < N; i++) { ws_set_box(w, &bx[i]); ensure_basic(w); t0 = now_ns(); ensure_es_mask(w, g_ALL); tes += now_ns() - t0;
        for (f = 0; f < g_n; f++) { t0 = now_ns(); model(w, f); model_aff(w, f); tmod[f] += now_ns() - t0; } }
    for (i = 0; i < N; i++) { t0 = now_ns(); claim(w, &bx[i], C_F); tF += now_ns() - t0; }
    printf("est %.0f ns, es(all) %.0f ns, F-claim %.0f ns\n", (double) test_ / N, (double) tes / N, (double) tF / N);
    for (f = 0; f < g_n; f++) printf("model f=%d %.0f ns\n", f, (double) tmod[f] / N);
    free(bx);
    ws_free(w);
}

/* ------------------------------------------------------------------ */

static uint64_t splitmix(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static int cmp_size_desc(const void *a, const void *b)
{
    const rec_t *x = a, *y = b;
    if (x->nbytes != y->nbytes) return x->nbytes < y->nbytes ? 1 : -1;
    return x->id < y->id ? -1 : x->id > y->id;
}

static void usage(void)
{
    fprintf(stderr,
        "usage: arbcheck2 [options] d tasks_file tree_file\n"
        "  -t N            threads (default 1)\n"
        "  --shard S/N     only records with id %% N == S\n"
        "  --ids A-B       only records with A <= id <= B\n"
        "  --sample K      only a pseudo-random 1/K of the records (by hash of id)\n"
        "  --seed X        seed for --sample (default 1)\n"
        "  --limit N       only the first N selected records in file order\n"
        "  --largest N     also check the N largest available records\n"
        "  --prec P        arb precision in bits (default 64)\n"
        "  --maxref R      maximal extra refinement depth (default 12)\n"
        "  --timing        time every leaf, report ns/leaf by type\n"
        "  --quiet         print only failing task lines and the summary\n"
        "  --selftest      compare the expansions with direct evaluation, then exit\n"
        "  --no-wform      do not use the w-form (T-form only, as arbcheck)\n"
        "Sample/limit/largest modes skip the missing-id check.\n");
    exit(2);
}

int main(int argc, char **argv)
{
    int nthreads = 1, shard = 0, nshard = 1, i, do_self = 0;
    uint64_t idlo = 0, idhi = UINT64_MAX, sample = 1, seed = 1, limit = 0, largest = 0;
    int sampling = 0;
    const char *tasks_path, *tree_path;
    int pos = 0;
    const char *posv[3];
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(a, "--shard") && i + 1 < argc) { if (sscanf(argv[++i], "%d/%d", &shard, &nshard) != 2 || nshard < 1 || shard < 0 || shard >= nshard) usage(); }
        else if (!strcmp(a, "--ids") && i + 1 < argc) { if (sscanf(argv[++i], "%" SCNu64 "-%" SCNu64, &idlo, &idhi) != 2) usage(); }
        else if (!strcmp(a, "--sample") && i + 1 < argc) { sample = strtoull(argv[++i], NULL, 10); sampling = 1; if (!sample) usage(); }
        else if (!strcmp(a, "--seed") && i + 1 < argc) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(a, "--limit") && i + 1 < argc) { limit = strtoull(argv[++i], NULL, 10); sampling = 1; }
        else if (!strcmp(a, "--largest") && i + 1 < argc) { largest = strtoull(argv[++i], NULL, 10); sampling = 1; }
        else if (!strcmp(a, "--prec") && i + 1 < argc) g_prec = atoi(argv[++i]);
        else if (!strcmp(a, "--maxref") && i + 1 < argc) g_maxref = atoi(argv[++i]);
        else if (!strcmp(a, "--timing")) g_timing = 1;
        else if (!strcmp(a, "--quiet")) g_quiet = 1;
        else if (!strcmp(a, "--selftest")) do_self = 1;
        else if (!strcmp(a, "--debug")) g_debug = 1;
        else if (!strcmp(a, "--no-wform")) g_wform = 0;
        else if (!strcmp(a, "--bench")) do_self = 2;
        else if (a[0] == '-' && a[1]) usage();
        else if (pos < 3) posv[pos++] = a;
        else usage();
    }
    if (pos < 1) usage();
    g_d = atoi(posv[0]);
    if (g_d < 4 || g_d > MAXNV + 2) { fprintf(stderr, "d must be in 4..%d\n", MAXNV + 2); return 2; }
    g_n = g_d - 1; g_nv = g_d - 2; g_D = 2 * g_nv; g_ALL = (1 << g_nv) - 1;
    if (g_prec < 32) g_prec = 32;
    if (nthreads < 1) nthreads = 1;
    init_constants();
    printf("VERSION %s, built with FLINT %s, w-form %s\n", ARBCHECK2_VERSION, flint_version, g_wform ? "on" : "off");
    fflush(stdout);
    if (do_self == 2) { bench(); return 0; }
    if (do_self) return selftest() ? 0 : 1;
    if (pos != 3) usage();
    tasks_path = posv[1]; tree_path = posv[2];
    g_tree_path = tree_path;

    if (!read_tasks(tasks_path)) return 2;

    /* scan record headers (only complete records are used) */
    rec_t *recs = NULL;
    size_t nrec = 0, caprec = 0;
    uint64_t nbadhdr = 0, ndup = 0, nunknown = 0, tail_incomplete = 0, nmissing = 0;
    unsigned char *seen = calloc((size_t) g_maxid + 1, 1);
    {
        FILE *fp = fopen(tree_path, "rb");
        struct stat sb;
        uint64_t off = 0, fsize;
        if (!fp || fstat(fileno(fp), &sb) != 0) { perror(tree_path); return 2; }
        fsize = (uint64_t) sb.st_size;
        for (;;) {
            unsigned char h[14];
            uint32_t id;
            uint64_t nb;
            size_t got;
            if (off == fsize) break;
            if (fsize - off < 14) { tail_incomplete = 1; break; }
            if (fseeko(fp, (off_t) off, SEEK_SET) != 0) { perror("fseeko"); return 2; }
            got = fread(h, 1, 14, fp);
            if (got != 14) { tail_incomplete = 1; break; }
            if (h[0] != 'T' || h[1] != 'K') {
                fprintf(stderr, "tree: bad record magic at offset %" PRIu64 "\n", off);
                nbadhdr++;
                break;
            }
            id = (uint32_t) h[2] | (uint32_t) h[3] << 8 | (uint32_t) h[4] << 16 | (uint32_t) h[5] << 24;
            nb = 0;
            for (i = 0; i < 8; i++) nb |= (uint64_t) h[6 + i] << (8 * i);
            if (nb > fsize - off - 14) { tail_incomplete = 1; break; }
            if (id > g_maxid || !g_hastask[id]) {
                fprintf(stderr, "tree: record id %u not in tasks file\n", id);
                nunknown++;
            } else if (seen[id]) {
                fprintf(stderr, "tree: duplicate record for id %u\n", id);
                ndup++;
            } else {
                seen[id] = 1;
                if (nrec == caprec) { caprec = caprec ? 2 * caprec : 4096; recs = realloc(recs, caprec * sizeof(rec_t)); }
                recs[nrec].id = id; recs[nrec].off = off + 14; recs[nrec].nbytes = nb;
                nrec++;
            }
            off += 14 + nb;
        }
        fclose(fp);
    }

    /* selection */
    g_jobs = malloc((nrec + 1) * sizeof(rec_t));
    g_njobs = 0;
    {
        unsigned char *pick = calloc(nrec + 1, 1);
        size_t r, cnt = 0;
        for (r = 0; r < nrec; r++) {
            uint32_t id = recs[r].id;
            if (id % (uint32_t) nshard != (uint32_t) shard || id < idlo || id > idhi) continue;
            if (sample > 1 && splitmix((uint64_t) id ^ (seed * 0x5851f42d4c957f2dULL)) % sample != 0) continue;
            if (limit && cnt >= limit) continue;
            pick[r] = 1;
            cnt++;
        }
        if (largest) {
            size_t *ix = malloc(nrec * sizeof(size_t) + 1);
            for (r = 0; r < nrec; r++) ix[r] = r;
            /* partial selection of the largest */
            for (r = 0; r < largest && r < nrec; r++) {
                size_t best = r, q;
                for (q = r + 1; q < nrec; q++) if (recs[ix[q]].nbytes > recs[ix[best]].nbytes) best = q;
                size_t t = ix[r]; ix[r] = ix[best]; ix[best] = t;
                pick[ix[r]] = 1;
            }
            free(ix);
        }
        for (r = 0; r < nrec; r++) if (pick[r]) g_jobs[g_njobs++] = recs[r];
        free(pick);
    }
    if (!sampling) {
        uint32_t id;
        for (id = 0; id <= g_maxid; id++) {
            if (!g_hastask[id]) continue;
            if (id % (uint32_t) nshard != (uint32_t) shard || id < idlo || id > idhi) continue;
            if (!seen[id]) {
                if (nmissing < 20) fprintf(stderr, "tree: no record for task id %u\n", id);
                nmissing++;
            }
        }
    }
    qsort(g_jobs, g_njobs, sizeof(rec_t), cmp_size_desc);
    {
        uint64_t tb = 0; size_t r;
        for (r = 0; r < g_njobs; r++) tb += g_jobs[r].nbytes;
        fprintf(stderr, "arbcheck2 d=%d prec=%ld maxref=%d threads=%d wform=%d: %" PRIu32 " tasks, %zu complete records%s, checking %zu records (%" PRIu64 " nodes)\n",
                g_d, (long) g_prec, g_maxref, nthreads, g_wform, g_ntask, nrec, tail_incomplete ? " (+ incomplete tail)" : "", g_njobs, tb);
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    {
        pthread_t *th = malloc(sizeof(pthread_t) * nthreads);
        for (i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, worker, NULL);
        for (i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
        free(th);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double wall = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    double cpu = ru.ru_utime.tv_sec + 1e-6 * ru.ru_utime.tv_usec + ru.ru_stime.tv_sec + 1e-6 * ru.ru_stime.tv_usec;

    uint64_t leaves = 0;
    for (i = 0; i < NCODE; i++) leaves += g_tot.leaves[i];
    int bad = g_tot.fail > 0 || nmissing > 0 || ndup > 0 || nunknown > 0 || nbadhdr > 0 || (!sampling && tail_incomplete);
    printf("SUMMARY d=%d records=%zu nodes=%" PRIu64 " leaves=%" PRIu64 " F=%" PRIu64 " E=%" PRIu64 " L=%" PRIu64
           " O=%" PRIu64 " S=%" PRIu64 " X=%" PRIu64 " U=%" PRIu64 " direct=%" PRIu64 " refined=%" PRIu64
           " subboxes=%" PRIu64 " fail=%" PRIu64 " failed_records=%" PRIu64 " missing=%" PRIu64 "%s dup=%" PRIu64
           " unknown=%" PRIu64 " badheader=%" PRIu64 " tail_incomplete=%" PRIu64 " wform_direct=%" PRIu64
           " (F=%" PRIu64 " E=%" PRIu64 " L=%" PRIu64 ") wform_subboxes=%" PRIu64 " maxdepth=%d wall=%.1fs cpu=%.1fs status=%s\n",
           g_d, g_njobs, g_tot.nodes, leaves, g_tot.leaves[0], g_tot.leaves[1], g_tot.leaves[2], g_tot.leaves[3],
           g_tot.leaves[4], g_tot.leaves[5], g_tot.leaves[6], g_tot.direct, g_tot.refined, g_tot.subboxes,
           g_tot.fail, g_fail_records, nmissing, sampling ? "(not checked: sample mode)" : "", ndup, nunknown, nbadhdr,
           tail_incomplete, g_tot.wdirect[0] + g_tot.wdirect[1] + g_tot.wdirect[2], g_tot.wdirect[0], g_tot.wdirect[1],
           g_tot.wdirect[2], g_tot.wsub, g_tot.maxdepth, wall, cpu, bad ? "FAIL" : "PASS");
    if (g_timing) {
        for (i = 0; i < NCODE; i++)
            if (g_tcnt[i])
                printf("TIMING %s leaves=%" PRIu64 " ns/leaf=%.0f leaves/s/thread=%.0f models/leaf=%.2f\n", CODE_NAME[i], g_tcnt[i],
                       (double) g_tns[i] / g_tcnt[i], 1e9 * g_tcnt[i] / (double) g_tns[i], (double) g_tmod[i] / g_tcnt[i]);
    }
    return bad ? 1 : 0;
}
