/* export_tree_v3.c — (copy of export_tree.c for c/smale_bb_v3.c) re-run the smale_bb_v3 subdivision for the tasks of a completed run and write the subdivision tree,
 * one byte per node, for an independent checker (c/arbcheck/).  The evaluator is c/smale_bb_v3.c, included verbatim
 * (its sha256 is bound into the run records); this file only adds a driver that records the tree.
 *
 * usage: export_tree D R_EXCL TASKS_FILE OUT_PREFIX NTHREADS
 *   TASKS_FILE  e.g. runs/d6_v2/d6.tasks (first line = header, then "id c_1 .. c_D" in hex floats; half-width from header)
 *   writes OUT_PREFIX.tree   : records  "TK" | uint32 id | uint64 nbytes | nbytes node codes   (little endian)
 *          OUT_PREFIX.counts : one line per task "id processed F E L outside sym excl unresolved maxdepth"
 * Node codes, depth-first pre-order, lower child (c - h/2) before upper child (c + h/2), exactly as bb_box:
 *   0 = F, 1 = E, 2 = L, 3 = outside, 4 = symmetry, 5 = excluded, 6 = unresolved, 16 + k = bisect coordinate k.
 * Only excl_mode 2.  Build: see c/arbcheck/README.md. */
#define main smale_main
#include "smale_bb_v3.c"
#undef main

typedef struct { unsigned char *b; size_t n, cap; } buf_t;
static void put(buf_t *B, unsigned char x) {
    if (B->n == B->cap) { B->cap = B->cap ? 2 * B->cap : 1 << 16; B->b = realloc(B->b, B->cap); if (!B->b) { perror("realloc"); exit(3); } }
    B->b[B->n++] = x;
}

static void tree_box(const prob_t *P, const double *c0, const double *h0, stats_t *st, buf_t *B) {
    int D = P->D;
    if (!in_domain(D, c0, h0)) { fprintf(stderr, "FATAL: top box outside the domain (A5)\n"); exit(4); }
    size_t cap = 4096, sp = 0;
    box_t *stk = malloc(cap * sizeof(box_t));
    memcpy(stk[0].c, c0, D * sizeof(double)); memcpy(stk[0].h, h0, D * sizeof(double)); stk[0].depth = 0; sp = 1;
    while (sp) {
        box_t b = stk[--sp];
        st->processed++;
        int g = geom_discard(P, b.c, b.h);
        if (g) { if (g == 1) { st->outside++; put(B, 3); } else { st->sym++; put(B, 4); } continue; }
        if (excluded(P, b.c, b.h)) { st->excl++; put(B, 5); continue; }
        int sk = 0;
        int r = eval_box(P, b.c, b.h, &sk);
        if (r) { if (r == 1) { st->F++; put(B, 0); } else if (r == 2) { st->E++; put(B, 1); } else { st->L++; put(B, 2); } continue; }
        if (b.h[sk] < 1e-9) { st->unresolved++; put(B, 6); continue; }
        put(B, (unsigned char)(16 + sk));
        if (sp + 2 > cap) { cap *= 2; stk = realloc(stk, cap * sizeof(box_t)); if (!stk) { perror("realloc"); exit(3); } }
        box_t b1 = b, b2 = b;
        double hh = b.h[sk] * 0.5;
        if (!(hh >= HMIN)) { fprintf(stderr, "FATAL: bisection below 2^-40\n"); exit(4); }
        b1.h[sk] = hh; b2.h[sk] = hh; b1.c[sk] = b.c[sk] - hh; b2.c[sk] = b.c[sk] + hh;
        b1.depth = b2.depth = b.depth + 1;
        if (b1.depth > st->maxdepth) st->maxdepth = b1.depth;
        stk[sp++] = b2; stk[sp++] = b1;
    }
    free(stk);
}

typedef struct {
    prob_t *P; int ntask; int *ids; double *tc; double th; int next;
    pthread_mutex_t mtx; FILE *tree, *counts;
} job_t;

static void *work(void *arg) {
    job_t *J = arg;
    int D = J->P->D;
    double h[DMAX];
    for (int k = 0; k < D; k++) h[k] = J->th;
    buf_t B = {0};
    for (;;) {
        pthread_mutex_lock(&J->mtx);
        int t = J->next < J->ntask ? J->next++ : -1;
        pthread_mutex_unlock(&J->mtx);
        if (t < 0) break;
        stats_t st; memset(&st, 0, sizeof st);
        B.n = 0;
        tree_box(J->P, J->tc + (size_t)t * D, h, &st, &B);
        uint32_t id = (uint32_t)J->ids[t]; uint64_t nb = B.n;
        pthread_mutex_lock(&J->mtx);
        fwrite("TK", 1, 2, J->tree); fwrite(&id, 4, 1, J->tree); fwrite(&nb, 8, 1, J->tree); fwrite(B.b, 1, B.n, J->tree);
        fprintf(J->counts, "%u %llu %llu %llu %llu %llu %llu %llu %llu %d\n", id, (unsigned long long)st.processed,
                (unsigned long long)st.F, (unsigned long long)st.E, (unsigned long long)st.L, (unsigned long long)st.outside,
                (unsigned long long)st.sym, (unsigned long long)st.excl, (unsigned long long)st.unresolved, st.maxdepth);
        pthread_mutex_unlock(&J->mtx);
    }
    free(B.b);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: export_tree D R_EXCL TASKS_FILE OUT_PREFIX NTHREADS\n"); return 2; }
    int d = atoi(argv[1]); double r = atof(argv[2]); int nth = atoi(argv[5]);
    prob_t P; setup(&P, d, r, 2);
    int D = P.D;
    FILE *f = fopen(argv[3], "r"); if (!f) { perror(argv[3]); return 1; }
    char *line = NULL; size_t cap = 0;
    if (getline(&line, &cap, f) < 0) return 1;
    char *hw = strstr(line, "half_width="); if (!hw) { fprintf(stderr, "no half_width in header\n"); return 1; }
    double th = strtod(hw + 11, NULL);
    int ntask = 0, tcap = 1024;
    int *ids = malloc(sizeof(int) * tcap); double *tc = malloc(sizeof(double) * tcap * D);
    while (getline(&line, &cap, f) > 0) {
        char *p = line, *e;
        long id = strtol(p, &e, 10); if (e == p) continue; p = e;
        if (ntask == tcap) { tcap *= 2; ids = realloc(ids, sizeof(int) * tcap); tc = realloc(tc, sizeof(double) * tcap * D); }
        ids[ntask] = (int)id;
        for (int k = 0; k < D; k++) { tc[(size_t)ntask * D + k] = strtod(p, &e); if (e == p) { fprintf(stderr, "bad task line\n"); return 1; } p = e; }
        ntask++;
    }
    fclose(f);
    char fn[1024];
    job_t J = {&P, ntask, ids, tc, th, 0, PTHREAD_MUTEX_INITIALIZER, NULL, NULL};
    snprintf(fn, sizeof fn, "%s.tree", argv[4]); J.tree = fopen(fn, "wb");
    snprintf(fn, sizeof fn, "%s.counts", argv[4]); J.counts = fopen(fn, "w");
    if (!J.tree || !J.counts) { perror("open output"); return 1; }
    pthread_t th_[64];
    if (nth < 1 || nth > 64) nth = 1;
    for (int i = 0; i < nth; i++) pthread_create(&th_[i], NULL, work, &J);
    for (int i = 0; i < nth; i++) pthread_join(th_[i], NULL);
    fclose(J.tree); fclose(J.counts);
    fprintf(stderr, "exported %d tasks\n", ntask);
    return 0;
}
