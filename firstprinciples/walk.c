/* walk.c -- independent tree walker for the first-principles check.
 *
 * Written from the record format description only.  It does NOT evaluate any
 * critical value; it only (a) parses every record and checks that it is one
 * complete pre-order binary tree, (b) reconstructs leaf boxes exactly (all
 * box coordinates are dyadic, so double arithmetic is exact while the split
 * depth per coordinate stays well below 50 -- checked), (c) checks the purely
 * geometric leaf claims (outside / symmetry / excluded) for EVERY such leaf in
 * double precision with a tolerance, writing borderline and failing leaves to
 * a file for an exact re-check in Python, and (d) in mode "sample" writes a
 * hash-selected subset of leaves (plus all leaves at depth >= a threshold) to
 * a binary file for evaluation of the true V_i in Python.
 *
 * usage:
 *   walk stats  TREE TASKBIN EQPTS OUT_PREFIX
 *   walk sample TREE TASKBIN EQPTS OUT_PREFIX PROBFILE DEEP_THRESHOLD SEED
 *
 * TASKBIN: int32 ntask, int32 D, double h, then ntask x (int32 id, double c[D])
 *          (ids are 0..ntask-1 and stored in order).
 * EQPTS:   int32 m, then m x double[D]  (equality points, double-rounded).
 * PROBFILE: ntask doubles, the per-task selection probability (sample mode).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAXD 16
#define MAXDEPTH 4096

static int D, n;            /* D = 2(n-1) */
static double H0;
static int ntask;
static double *tcen;        /* ntask x D */
static int neq;
static double *eqp;         /* neq x D */

static uint64_t splitmix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static double minmod2(double a1, double a2, double b1, double b2) {
    double mx = (a1 <= 0 && 0 <= a2) ? 0 : fmin(fabs(a1), fabs(a2));
    double my = (b1 <= 0 && 0 <= b2) ? 0 : fmin(fabs(b1), fabs(b2));
    return mx * mx + my * my;
}
static double maxmod2(double a1, double a2, double b1, double b2) {
    double mx = fmax(fabs(a1), fabs(a2)), my = fmax(fabs(b1), fabs(b2));
    return mx * mx + my * my;
}

/* margins: >0 means claim holds (by that much, in double) */
static double margin_outside(const double *lo, const double *hi) {
    double best = -1e300;
    for (int v = 0; v < n - 1; v++) {
        double m = minmod2(lo[2*v], hi[2*v], lo[2*v+1], hi[2*v+1]) - 1.0;
        if (m > best) best = m;
    }
    return best;
}
static double margin_symmetry(const double *lo, const double *hi) {
    double best = -hi[1];                      /* Im u_2 < 0 on all of B */
    for (int v = 0; v + 1 < n - 1; v++) {      /* |u_v| < |u_{v+1}| */
        double m = minmod2(lo[2*v+2], hi[2*v+2], lo[2*v+3], hi[2*v+3])
                 - maxmod2(lo[2*v], hi[2*v], lo[2*v+1], hi[2*v+1]);
        if (m > best) best = m;
    }
    return best;
}
static double margin_excluded(const double *lo, const double *hi, int *which) {
    double best = -1e300; *which = -1;
    for (int e = 0; e < neq; e++) {
        const double *p = eqp + (size_t)e * D;
        double s = 0;
        for (int k = 0; k < D; k++) {
            double t = fmax(fabs(lo[k] - p[k]), fabs(hi[k] - p[k]));
            s += t * t;
        }
        double m = 1.0 / 400.0 - s;
        if (m > best) { best = m; *which = e; }
    }
    return best;
}
static int contains_zero_u(const double *lo, const double *hi) {
    for (int v = 0; v < n - 1; v++)
        if (lo[2*v] <= 0 && 0 <= hi[2*v] && lo[2*v+1] <= 0 && 0 <= hi[2*v+1]) return 1;
    return 0;
}

typedef struct { int k; double c, h; int state; } Frame;

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 2; }
    int sample = strcmp(argv[1], "sample") == 0;
    const char *treep = argv[2], *taskp = argv[3], *eqpp = argv[4], *outp = argv[5];
    FILE *f = fopen(taskp, "rb");
    int32_t i32;
    fread(&i32, 4, 1, f); ntask = i32; fread(&i32, 4, 1, f); D = i32; fread(&H0, 8, 1, f);
    n = D / 2 + 1;
    tcen = malloc(sizeof(double) * (size_t)ntask * D);
    for (int t = 0; t < ntask; t++) {
        fread(&i32, 4, 1, f);
        if (i32 != t) { fprintf(stderr, "task id order\n"); return 2; }
        fread(tcen + (size_t)t * D, 8, D, f);
    }
    fclose(f);
    f = fopen(eqpp, "rb"); fread(&i32, 4, 1, f); neq = i32;
    eqp = malloc(sizeof(double) * (size_t)neq * D);
    fread(eqp, 8, (size_t)neq * D, f); fclose(f);

    double *prob = NULL; int deepT = 1 << 30; uint64_t seed = 1;
    if (sample) {
        prob = malloc(sizeof(double) * ntask);
        f = fopen(argv[6], "rb"); fread(prob, 8, ntask, f); fclose(f);
        deepT = atoi(argv[7]); seed = strtoull(argv[8], 0, 10);
    }
    char path[4096];
    snprintf(path, sizeof path, "%s.%s.tasks.txt", outp, sample ? "sample" : "stats");
    FILE *ftask = fopen(path, "w");
    snprintf(path, sizeof path, "%s.%s.geom_flag.bin", outp, sample ? "sample" : "stats");
    FILE *fgeom = sample ? NULL : fopen(path, "wb");
    snprintf(path, sizeof path, "%s.sample.leaves.bin", outp);
    FILE *fs = sample ? fopen(path, "wb") : NULL;

    FILE *ft = fopen(treep, "rb");
    static unsigned char buf[1 << 20];
    unsigned char *seen = calloc(ntask, 1);
    long long depth_hist[7][MAXDEPTH]; memset(depth_hist, 0, sizeof depth_hist);
    long long tot[7] = {0};
    double gmin[3] = {1e300, 1e300, 1e300};   /* smallest positive margin per geometric type */
    long long gfail[3] = {0}, gnear[3] = {0}, zero_in[7] = {0};
    long long nrec = 0, bad = 0, nsel = 0;
    int maxsplit_global = 0, maxdepth_global = 0;
    const double TOL = getenv("WALK_TOL") ? atof(getenv("WALK_TOL")) : 1e-12;

    for (;;) {
        unsigned char hdr[14];
        size_t got = fread(hdr, 1, 14, ft);
        if (got == 0) break;
        if (got != 14 || hdr[0] != 'T' || hdr[1] != 'K') { fprintf(stderr, "bad header at record %lld\n", nrec); bad++; break; }
        uint32_t id; uint64_t nb;
        memcpy(&id, hdr + 2, 4); memcpy(&nb, hdr + 6, 8);   /* little endian host */
        nrec++;
        if (id >= (uint32_t)ntask) { fprintf(stderr, "id out of range %u\n", id); bad++; break; }
        if (seen[id]) { fprintf(stderr, "duplicate id %u\n", id); bad++; }
        seen[id] = 1;
        double cen[MAXD], half[MAXD], lo[MAXD], hi[MAXD];
        int nsplit[MAXD] = {0};
        for (int k = 0; k < D; k++) { cen[k] = tcen[(size_t)id * D + k]; half[k] = H0; }
        static Frame st[MAXDEPTH]; int sp = 0;
        long long cnt[7] = {0}; long long leafidx = 0; int complete = 0, extra = 0, badcode = 0;
        int maxdepth = 0, maxsplit = 0;
        double p = sample ? prob[id] : 0;
        uint64_t thr = (p >= 1.0) ? UINT64_MAX : (uint64_t)(p * 18446744073709551616.0);
        uint64_t left = nb;
        while (left > 0) {
            size_t want = left < sizeof buf ? (size_t)left : sizeof buf;
            size_t g = fread(buf, 1, want, ft);
            if (g != want) { fprintf(stderr, "truncated record %u\n", id); bad++; left = 0; break; }
            left -= g;
            for (size_t b = 0; b < g; b++) {
                int code = buf[b];
                if (complete) { extra++; continue; }
                if (code >= 16 && code < 16 + D) {
                    int k = code - 16;
                    if (sp >= MAXDEPTH) { fprintf(stderr, "too deep\n"); return 3; }
                    st[sp].k = k; st[sp].c = cen[k]; st[sp].h = half[k]; st[sp].state = 0; sp++;
                    nsplit[k]++; if (nsplit[k] > maxsplit) maxsplit = nsplit[k];
                    half[k] = half[k] / 2; cen[k] = cen[k] - half[k];
                    continue;
                }
                if (code > 6) { badcode++; code = 6; }
                /* leaf */
                cnt[code]++;
                int depth = sp; if (depth > maxdepth) maxdepth = depth;
                depth_hist[code][depth < MAXDEPTH ? depth : MAXDEPTH - 1]++;
                for (int k = 0; k < D; k++) { lo[k] = cen[k] - half[k]; hi[k] = cen[k] + half[k]; }
                if (code <= 5 && code != 3 && code != 4 && contains_zero_u(lo, hi)) zero_in[code]++;
                if (!sample && code >= 3 && code <= 5) {
                    int g3 = code - 3, which = -1; double m;
                    if (code == 3) m = margin_outside(lo, hi);
                    else if (code == 4) m = margin_symmetry(lo, hi);
                    else m = margin_excluded(lo, hi, &which);
                    if (m > 0 && m < gmin[g3]) gmin[g3] = m;
                    if (m <= TOL) {           /* borderline or failing: exact recheck later */
                        if (m < -TOL) gfail[g3]++; else gnear[g3]++;
                        uint32_t rec[2] = {id, (uint32_t)code};
                        fwrite(rec, 4, 2, fgeom); fwrite(&leafidx, 8, 1, fgeom);
                        fwrite(cen, 8, D, fgeom); fwrite(half, 8, D, fgeom);
                    }
                }
                if (sample && code <= 2) {
                    uint64_t hsh = splitmix(splitmix(((uint64_t)id << 40) ^ (uint64_t)leafidx) ^ seed);
                    if (hsh < thr || depth >= deepT) {
                        uint32_t rec[2] = {id, (uint32_t)code};
                        uint32_t dd[2] = {(uint32_t)depth, (uint32_t)(depth >= deepT)};
                        fwrite(rec, 4, 2, fs); fwrite(&leafidx, 8, 1, fs); fwrite(dd, 4, 2, fs);
                        fwrite(cen, 8, D, fs); fwrite(half, 8, D, fs);
                        nsel++;
                    }
                }
                leafidx++;
                /* pop */
                for (;;) {
                    if (sp == 0) { complete = 1; break; }
                    Frame *fr = &st[sp - 1];
                    if (fr->state == 0) {
                        fr->state = 1; half[fr->k] = fr->h / 2; cen[fr->k] = fr->c + half[fr->k];
                        break;
                    }
                    cen[fr->k] = fr->c; half[fr->k] = fr->h; nsplit[fr->k]--; sp--;
                }
            }
        }
        if (maxdepth > maxdepth_global) maxdepth_global = maxdepth;
        if (maxsplit > maxsplit_global) maxsplit_global = maxsplit;
        for (int c = 0; c < 7; c++) tot[c] += cnt[c];
        if (!complete || extra || badcode) bad++;
        fprintf(ftask, "%u %llu %d %d %d %d", id, (unsigned long long)nb, complete, extra, badcode, maxdepth);
        for (int c = 0; c < 7; c++) fprintf(ftask, " %lld", cnt[c]);
        fprintf(ftask, " %d\n", maxsplit);
    }
    fclose(ft);
    long long missing = 0; for (int t = 0; t < ntask; t++) if (!seen[t]) missing++;
    snprintf(path, sizeof path, "%s.%s.summary.txt", outp, sample ? "sample" : "stats");
    FILE *fsum = fopen(path, "w");
    fprintf(fsum, "records %lld bad %lld missing_ids %lld maxdepth %d maxsplit_per_coord %d selected %lld\n",
            nrec, bad, missing, maxdepth_global, maxsplit_global, nsel);
    fprintf(fsum, "leaves F %lld E %lld L %lld outside %lld symmetry %lld excluded %lld unresolved %lld\n",
            tot[0], tot[1], tot[2], tot[3], tot[4], tot[5], tot[6]);
    fprintf(fsum, "F/E/L/excluded leaves whose box contains some u_j=0: F %lld E %lld L %lld excl %lld\n",
            zero_in[0], zero_in[1], zero_in[2], zero_in[5]);
    if (!sample) {
        const char *nm[3] = {"outside", "symmetry", "excluded"};
        for (int g = 0; g < 3; g++)
            fprintf(fsum, "geom %s: double-fail %lld borderline %lld smallest_positive_margin %.6e\n",
                    nm[g], gfail[g], gnear[g], gmin[g]);
        fprintf(fsum, "depth_hist (depth: F E L out sym excl unres)\n");
        for (int dpt = 0; dpt < MAXDEPTH; dpt++) {
            long long s = 0; for (int c = 0; c < 7; c++) s += depth_hist[c][dpt];
            if (!s) continue;
            fprintf(fsum, "%d:", dpt);
            for (int c = 0; c < 7; c++) fprintf(fsum, " %lld", depth_hist[c][dpt]);
            fprintf(fsum, "\n");
        }
    }
    fclose(fsum); fclose(ftask); if (fgeom) fclose(fgeom); if (fs) fclose(fs);
    return bad ? 1 : 0;
}
