/* prim_test_v2.c — adversarial inputs for the v2 primitives (cb_mul/cb_addto/cb_subto, cmod_bounds, cmod_hi_fast,
   powu, cr_log).  Output (hex) is checked with exact rationals / mpmath by c/regress/prim_check_v2.py.
   Build: see README.md (Reproduce).  Inputs cover all binades incl. subnormals, and cancellation. */
#define main smale_main
#include "../smale_bb_v2.c"
#undef main
static uint64_t s_ = 88172645463325252ULL;
static uint64_t rn(void){ s_ ^= s_ << 13; s_ ^= s_ >> 7; s_ ^= s_ << 17; return s_; }
static double rd(int maxexp){
  int k = rn() % 6; double x;
  switch (k) {
    case 0: { uint64_t u = rn() & 0x000fffffffffffffULL; memcpy(&x, &u, 8); break; }          /* subnormal */
    case 1: x = ldexp((rn() >> 11) * 0x1p-53, -(int)(rn() % 1100)); break;
    case 2: x = ldexp((rn() >> 11) * 0x1p-53, -(int)(rn() % 40)); break;
    case 3: x = 1.0 + ldexp((double)(rn() % 1000), -52); break;
    case 4: x = ldexp(1.0, -(int)(rn() % 1074)); break;
    default: x = ldexp((rn() >> 11) * 0x1p-53, (int)(rn() % (maxexp + 1))); break;
  }
  return (rn() & 1) ? -x : x;
}
int main(int argc, char **argv){
  int N = argc > 1 ? atoi(argv[1]) : 200000;
  for (int k = 0; k < N; k++) {
    cb a = {rd(0), rd(0), 0}, b = {rd(0), rd(0), 0};
    if (k % 3 == 1) { b.re = -a.re * (1 + ldexp((double)(rn()%8), -52)); b.im = -a.im; } /* cancellation in add */
    cb r = cb_mul(a, b);
    printf("mul %a %a %a %a %a %a\nR %a\n", a.re, a.im, b.re, b.im, r.re, r.im, r.rad);
    cb q = a; cb_addto(&q, b); printf("add %a %a %a %a %a %a\nR %a\n", a.re, a.im, b.re, b.im, q.re, q.im, q.rad);
    q = a; cb_subto(&q, b); printf("sub %a %a %a %a %a %a\nR %a\n", a.re, a.im, b.re, b.im, q.re, q.im, q.rad);
    double x = rd(400), y = rd(400), lo, hi;
    cmod_bounds(x, y, &lo, &hi);
    printf("mod %a %a %a %a %a\n", x, y, lo, hi, cmod_hi_fast(x, y));
    double z = fabs(rd(1)); int e = rn() % 9;
    printf("pow %a %d %a\n", z, e, powu(z, e));
    double w = fabs(rd(600)); if (w == 0) w = 1;
    printf("log %a %a\n", w, cr_log(w));
  }
  return 0;
}
