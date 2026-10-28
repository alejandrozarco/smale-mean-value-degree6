/* run_tree_v3.c — a certification run of c/smale_bb_v3.c that also writes the subdivision tree, so that a large run
 * (d = 7) need not be repeated by c/export_tree_v3.c for the independent check (c/arbcheck/arbcheck2).
 * The evaluator is c/smale_bb_v3.c, included verbatim: build with -DSRC_SHA256_RAW=<sha256 of c/smale_bb_v3.c>, so the
 * configuration string, its hash and every checkpoint record are exactly those of `smale_bb_v3 run`.  This file only
 * adds a driver.
 *
 * usage: run_tree_v3 D R_EXCL S NTHREADS OUTPREFIX     (excl_mode 2, seed 1; partition with SMALE_PART="k1,k2,.../N")
 *   OUTPREFIX.tasks       must exist (written by `smale_bb_v3 run` for the same configuration); its header is verified
 *   OUTPREFIX.done        checkpoint, same format as `smale_bb_v3 run` (shared: tasks done by either program count)
 *   OUTPREFIX.unresolved  as `smale_bb_v3 run`
 *   OUTPREFIX.tree        records "TK" | uint32 id | uint64 nbytes | node codes (format of c/export_tree_v3.c), one per
 *                         task finished by THIS program.  Trees of tasks finished by `smale_bb_v3 run` are written later
 *                         by c/export_tree_v3.c from a tasks file listing only those ids.
 * Per task: the tree record is appended and synced, then the .done record (one write(2), O_APPEND, synced).  On
 * resume, the .done file is validated as by `run` (checksums, configuration, ids, uniqueness, tree count), a torn tail
 * is truncated, and the .tree file is rewritten to keep only complete records whose id has a .done record (first copy).
 * The per-box logic is that of bb_box (stats, volume, unresolved lines, 2^-40 limit) plus the tree codes of
 * c/export_tree_v3.c; c/regress/cmp_export_counts.py and c/check_done.py check the outputs against each other. */
#define main smale_main
#include "smale_bb_v3.c"
#undef main

typedef struct { unsigned char *b; size_t n, cap; } tbuf_t;
static void tput(tbuf_t *B, unsigned char x) {
    if (B->n == B->cap) { B->cap = B->cap ? 2 * B->cap : 1 << 16; B->b = realloc(B->b, B->cap); if (!B->b) { perror("realloc"); exit(3); } }
    B->b[B->n++] = x;
}

typedef struct {
    const prob_t *P; int ntask; double *tc; double th; int *done; int next;
    pthread_mutex_t mtx; int done_fd, tree_fd, unres_fd; const char *cfg16; long ndone_here;
} rjob_t;

static void tree_bb(const prob_t *P, const double *c0, const double *h0, stats_t *st, tbuf_t *B, rjob_t *J) {
    int D = P->D;
    if (!in_domain(D, c0, h0)) { fprintf(stderr, "FATAL: top box outside the domain hypotheses (A5)\n"); exit(4); }
    size_t cap = 4096, sp = 1;
    box_t *stk = malloc(cap * sizeof(box_t));
    if (!stk) { perror("malloc"); exit(3); }
    memcpy(stk[0].c, c0, D * sizeof(double)); memcpy(stk[0].h, h0, D * sizeof(double)); stk[0].depth = 0;
    while (sp) {
        box_t b = stk[--sp];
        st->processed++;
        double vol = 1; for (int k = 0; k < D; k++) vol *= 2 * b.h[k];
        int g = geom_discard(P, b.c, b.h);
        if (g) { if (g == 1) { st->outside++; tput(B, 3); } else { st->sym++; tput(B, 4); } st->vol_done += vol; continue; }
        if (excluded(P, b.c, b.h)) { st->excl++; tput(B, 5); st->vol_done += vol; continue; }
        int sk = 0;
        int r = eval_box(P, b.c, b.h, &sk);
        if (r) {
            if (r == 1) { st->F++; tput(B, 0); } else if (r == 2) { st->E++; tput(B, 1); } else { st->L++; tput(B, 2); }
            st->vol_done += vol; continue;
        }
        if (b.h[sk] < 1e-9) {
            st->unresolved++; st->vol_done += vol; tput(B, 6);
            char line[1024]; int L = snprintf(line, sizeof line, "UNRESOLVED cfg=%s", J->cfg16);
            for (int k = 0; k < D; k++) L += snprintf(line + L, sizeof line - L, " %a", b.c[k]);
            for (int k = 0; k < D; k++) L += snprintf(line + L, sizeof line - L, " %a", b.h[k]);
            L += snprintf(line + L, sizeof line - L, "\n");
            pthread_mutex_lock(&J->mtx); write_all(J->unres_fd, line, (size_t)L); full_sync(J->unres_fd); pthread_mutex_unlock(&J->mtx);
            continue;
        }
        tput(B, (unsigned char)(16 + sk));
        if (sp + 2 > cap) { cap *= 2; stk = realloc(stk, cap * sizeof(box_t)); if (!stk) { perror("realloc"); exit(3); } }
        box_t b1 = b, b2 = b;
        double hh = b.h[sk] * 0.5;
        if (!(hh >= HMIN)) { fprintf(stderr, "FATAL: bisection below 2^-40\n"); exit(4); }
        b1.h[sk] = hh; b2.h[sk] = hh; b1.c[sk] = b.c[sk] - hh; b2.c[sk] = b.c[sk] + hh;
        b1.depth = b2.depth = b.depth + 1;
        if (b1.depth > st->maxdepth) st->maxdepth = b1.depth;
        stk[sp++] = b2; stk[sp++] = b1;     /* lower child on top: popped first (pre-order, lower subtree first) */
    }
    free(stk);
}

static void *rworker(void *arg) {
    rjob_t *J = arg;
    const prob_t *P = J->P;
    double h[DMAX];
    for (int k = 0; k < P->D; k++) h[k] = J->th;
    tbuf_t B = {0};
    for (;;) {
        pthread_mutex_lock(&J->mtx);
        while (J->next < J->ntask && J->done[J->next]) J->next++;
        int t = J->next < J->ntask ? J->next++ : -1;
        pthread_mutex_unlock(&J->mtx);
        if (t < 0) break;
        stats_t st; memset(&st, 0, sizeof st);
        B.n = 0;
        double t1 = now();
        tree_bb(P, J->tc + (size_t)t * P->D, h, &st, &B, J);
        double el = now() - t1;
        char payload[512], line[600];
        snprintf(payload, sizeof payload, "R %d %llu %llu %llu %llu %llu %llu %llu %llu %d %.3f %.17g cfg=%s", t,
                 (unsigned long long)st.processed, (unsigned long long)st.F, (unsigned long long)st.E,
                 (unsigned long long)st.L, (unsigned long long)st.outside, (unsigned long long)st.sym,
                 (unsigned long long)st.excl, (unsigned long long)st.unresolved, st.maxdepth, el, st.vol_done, J->cfg16);
        int L = make_line(line, sizeof line, payload);
        unsigned char hdr[14]; uint32_t id = (uint32_t)t; uint64_t nb = B.n;
        hdr[0] = 'T'; hdr[1] = 'K'; memcpy(hdr + 2, &id, 4); memcpy(hdr + 6, &nb, 8);   /* little endian hosts only */
        pthread_mutex_lock(&J->mtx);
        write_all(J->tree_fd, (const char *)hdr, 14); write_all(J->tree_fd, (const char *)B.b, B.n);
        full_sync(J->tree_fd);
        write_all(J->done_fd, line, (size_t)L);
        full_sync(J->done_fd);
        J->ndone_here++;
        pthread_mutex_unlock(&J->mtx);
    }
    free(B.b);
    return NULL;
}

static int parse_part(const char *ps, int ntask, int *done) {
    const char *sl = strchr(ps, '/');
    if (!sl) return -1;
    char *e; long N = strtol(sl + 1, &e, 10);
    if (e == sl + 1 || *e || N < 1 || N > 4096) return -1;
    char *in = calloc((size_t)N, 1); const char *q = ps; int nk = 0;
    if (q == sl) { free(in); return -1; }
    while (q < sl) {
        long k = strtol(q, &e, 10);
        if (e == q || k < 0 || k >= N || e > sl) { free(in); return -1; }
        in[k] = 1; nk++;
        if (e == sl) break;
        if (*e != ',') { free(in); return -1; }
        q = e + 1;
        if (q == sl) break;                      /* tolerate one trailing comma */
    }
    for (int t = 0; t < ntask; t++) if (!in[t % N] && !done[t]) done[t] = 2;
    free(in);
    return nk;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: run_tree_v3 D R_EXCL S NTHREADS OUTPREFIX   [env SMALE_PART=k1,k2,.../N]\n"); return 2; }
    int d = atoi(argv[1]); double r = atof(argv[2]); int s = atoi(argv[3]), nth = atoi(argv[4]);
    const char *pre = argv[5];
    int em = 2; unsigned long long seed = 1;
    if (!(r > 0 && r <= 0.05)) { fprintf(stderr, "REFUSED: r_excl must be in (0, 0.05]\n"); return 5; }
    if (s < 1 || s > (1 << 20) || (s & (s - 1))) { fprintf(stderr, "REFUSED: s must be a power of two\n"); return 5; }
    if (nth < 1 || nth > 64) { fprintf(stderr, "nthreads must be in 1..64\n"); return 1; }
    prob_t P; setup(&P, d, r, em);
    int D = P.D;
    char cfgs[512], cfghex[65], cfg16[17];
    cfg_string(cfgs, sizeof cfgs, d, r, em, s, seed);
    sha256_hex(cfgs, strlen(cfgs), cfghex);
    memcpy(cfg16, cfghex, 16); cfg16[16] = 0;
    char fn[1024];
    /* tasks file: must exist, header must match this configuration; ids 0..ntask-1 in order */
    snprintf(fn, sizeof fn, "%s.tasks", pre);
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); return 3; }
    char *ln = NULL; size_t lcap = 0;
    if (getline(&ln, &lcap, f) < 0) { fprintf(stderr, "empty tasks file\n"); return 3; }
    char want[1024]; int ntask = 0; double th = 0;
    {
        char hw[64];
        if (sscanf(ln, "# smale_bb_v3 tasks cfg=%*s %*s ntask=%d half_width=%63s", &ntask, hw) != 2) { fprintf(stderr, "REFUSED: tasks header\n"); return 6; }
        th = strtod(hw, NULL);
        if (th != 1.0 / s) { fprintf(stderr, "REFUSED: tasks half_width %a != 1/s\n", th); return 6; }   /* review finding 5 */
        snprintf(want, sizeof want, "# smale_bb_v3 tasks cfg=%s %s ntask=%d half_width=%a\n", cfg16, cfgs, ntask, th);
        if (strcmp(ln, want)) { fprintf(stderr, "REFUSED: %s header is for a different configuration\n", fn); return 6; }
    }
    double *tc = malloc(sizeof(double) * (size_t)ntask * D);
    for (int t = 0; t < ntask; t++) {
        if (getline(&ln, &lcap, f) <= 0) { fprintf(stderr, "REFUSED: tasks file truncated\n"); return 6; }
        char *p = ln, *e; long id = strtol(p, &e, 10);
        if (e == p || id != t) { fprintf(stderr, "REFUSED: tasks file id order\n"); return 6; }
        p = e;
        for (int k = 0; k < D; k++) { tc[(size_t)t * D + k] = strtod(p, &e); if (e == p) { fprintf(stderr, "REFUSED: bad task line\n"); return 6; } p = e; }
    }
    fclose(f);
    /* checkpoint: validate as `run` does */
    int *done = calloc((size_t)ntask, sizeof(int));
    snprintf(fn, sizeof fn, "%s.done", pre);
    int fd = open(fn, O_RDWR | O_CREAT | O_APPEND, 0644);
    if (fd < 0) { perror(fn); return 3; }
    if (flock(fd, LOCK_EX | LOCK_NB) != 0) { fprintf(stderr, "REFUSED: %s is locked by another process\n", fn); return 6; }
    struct stat sb; fstat(fd, &sb);
    size_t fl = (size_t)sb.st_size;
    char *fb = malloc(fl + 1);
    if (fl && pread(fd, fb, fl, 0) != (ssize_t)fl) { perror("pread"); return 3; }
    size_t good = 0; int ndone = 0, have_header = 0, lineno = 0;
    char payload[1024];
    while (good < fl) {
        char *nl = memchr(fb + good, '\n', fl - good);
        if (!nl) break;
        size_t len = (size_t)(nl - (fb + good));
        lineno++;
        if (!check_line(fb + good, len, payload, sizeof payload)) { fprintf(stderr, "REFUSED: %s line %d: bad checksum\n", fn, lineno); return 7; }
        if (lineno == 1) {
            snprintf(want, sizeof want, "H smale_bb_v3 cfg=%s %s", cfg16, cfgs);
            if (strcmp(payload, want)) { fprintf(stderr, "REFUSED: %s header is for a different configuration\n", fn); return 6; }
            have_header = 1;
        } else {
            int t; unsigned long long pr, F_, E_, L_, o, sy, ex, un; int md; double el, vd; char cf[64]; int pos = 0;
            if (sscanf(payload, "R %d %llu %llu %llu %llu %llu %llu %llu %llu %d %lf %lf cfg=%63s%n", &t, &pr, &F_, &E_, &L_,
                       &o, &sy, &ex, &un, &md, &el, &vd, cf, &pos) != 13 || payload[pos] != 0) { fprintf(stderr, "REFUSED: %s line %d: malformed\n", fn, lineno); return 7; }
            if (strcmp(cf, cfg16) || t < 0 || t >= ntask || done[t]) { fprintf(stderr, "REFUSED: %s line %d: bad record\n", fn, lineno); return 7; }
            if (pr != 2 * (F_ + E_ + L_ + o + sy + ex + un) - 1) { fprintf(stderr, "REFUSED: %s line %d: tree count\n", fn, lineno); return 7; }
            done[t] = 1; ndone++;
        }
        good += len + 1;
    }
    if (good < fl) { fprintf(stderr, "note: truncating torn trailing line of %s\n", fn); if (ftruncate(fd, (off_t)good)) { perror("ftruncate"); return 3; } full_sync(fd); }
    free(fb);
    if (!have_header) {
        char pl[1024], line[1100];
        snprintf(pl, sizeof pl, "H smale_bb_v3 cfg=%s %s", cfg16, cfgs);
        int L = make_line(line, sizeof line, pl);
        write_all(fd, line, (size_t)L); full_sync(fd);
    }
    /* tree file: keep complete records with a .done record, first copy only; rewrite via a temporary file */
    snprintf(fn, sizeof fn, "%s.tree", pre);
    {
        char tmp[1100]; snprintf(tmp, sizeof tmp, "%s.tmp", fn);
        FILE *in = fopen(fn, "rb");
        long kept = 0, dropped = 0;
        if (in) {
            FILE *out = fopen(tmp, "wb"); if (!out) { perror(tmp); return 3; }
            char *seen = calloc((size_t)ntask, 1);
            unsigned char hdr[14];
            unsigned char *buf = NULL; size_t bcap = 0;
            while (fread(hdr, 1, 14, in) == 14) {
                uint32_t id; uint64_t nb; memcpy(&id, hdr + 2, 4); memcpy(&nb, hdr + 6, 8);
                if (hdr[0] != 'T' || hdr[1] != 'K') { fprintf(stderr, "REFUSED: corrupt tree file\n"); return 7; }
                if (nb > bcap) { bcap = nb; buf = realloc(buf, bcap); if (!buf) { perror("realloc"); return 3; } }
                if (fread(buf, 1, nb, in) != nb) break;                         /* torn last record */
                if (id < (uint32_t)ntask && done[id] == 1 && !seen[id]) { seen[id] = 1; fwrite(hdr, 1, 14, out); fwrite(buf, 1, nb, out); kept++; }
                else dropped++;
            }
            fclose(in);
            if (fflush(out) || fsync(fileno(out))) { perror("tree rewrite"); return 3; }
            fclose(out);
            if (rename(tmp, fn)) { perror("rename"); return 3; }
            free(seen); free(buf);
        }
        fprintf(stderr, "tree file: kept %ld records, dropped %ld\n", kept, dropped);
    }
    int tfd = open(fn, O_WRONLY | O_CREAT | O_APPEND, 0644);
    if (tfd < 0) { perror(fn); return 3; }
    int npart = 0;
    if (getenv("SMALE_PART")) {
        npart = parse_part(getenv("SMALE_PART"), ntask, done);
        if (npart < 0) { fprintf(stderr, "REFUSED: SMALE_PART must be k1,k2,.../N with 0 <= k < N <= 4096\n"); return 5; }
    }
    snprintf(fn, sizeof fn, "%s.unresolved", pre);
    int ufd = open(fn, O_WRONLY | O_CREAT | O_APPEND, 0644);
    if (ufd < 0) { perror(fn); return 3; }
    fprintf(stderr, "run_tree_v3 %s d=%d cfg=%s: %d tasks, %d already done, partition %s\n", SRC_SHA256, d, cfg16, ntask,
            ndone, getenv("SMALE_PART") ? getenv("SMALE_PART") : "(none)");
    rjob_t J = {&P, ntask, tc, th, done, 0, PTHREAD_MUTEX_INITIALIZER, fd, tfd, ufd, cfg16, 0};
    pthread_t th_[64];
    for (int k = 0; k < nth; k++) pthread_create(&th_[k], NULL, rworker, &J);
    for (int k = 0; k < nth; k++) pthread_join(th_[k], NULL);
    full_sync(fd); close(fd); full_sync(tfd); close(tfd); close(ufd);
    fprintf(stderr, "run_tree_v3: %ld tasks finished in this session\n", J.ndone_here);
    struct stat us; stat(fn, &us);                                   /* fn = OUTPREFIX.unresolved */
    return us.st_size > 0 ? 2 : 0;                                   /* review finding 5: as `smale_bb_v3 run` */
}
