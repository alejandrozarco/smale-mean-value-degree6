/* fpenv_probe.c — record the floating-point environment a run relies on (c/ERROR_ANALYSIS_v2.md A.1): rounding mode
 * round-to-nearest, no flush-to-zero / denormals-are-zero, FLT_EVAL_METHOD 0, IEEE sqrt and gradual underflow.
 * Build with the flags of c/build_v3.sh; run on the machine of the run; exit status 0 iff all checks pass.
 * Note: threads inherit the floating-point control register from the creating thread; smale_bb_v3 does not change it. */
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
int main(void) {
    int ok = 1;
    printf("FLT_EVAL_METHOD %d\n", (int)FLT_EVAL_METHOD); ok &= FLT_EVAL_METHOD == 0;
    int rm = fegetround(); printf("fegetround %s\n", rm == FE_TONEAREST ? "FE_TONEAREST" : "OTHER"); ok &= rm == FE_TONEAREST;
#if defined(__aarch64__)
    uint64_t fpcr; __asm__ volatile("mrs %0, fpcr" : "=r"(fpcr));
    printf("FPCR 0x%llx (FZ bit 24 = %d, RMode bits 22-23 = %d)\n", (unsigned long long)fpcr, (int)((fpcr >> 24) & 1), (int)((fpcr >> 22) & 3));
    ok &= ((fpcr >> 24) & 1) == 0 && ((fpcr >> 22) & 3) == 0;
#elif defined(__x86_64__)
    unsigned int mxcsr; __asm__ volatile("stmxcsr %0" : "=m"(mxcsr));
    printf("MXCSR 0x%x (FTZ bit 15 = %u, DAZ bit 6 = %u, RC bits 13-14 = %u)\n", mxcsr, (mxcsr >> 15) & 1, (mxcsr >> 6) & 1, (mxcsr >> 13) & 3);
    ok &= ((mxcsr >> 15) & 1) == 0 && ((mxcsr >> 6) & 1) == 0 && ((mxcsr >> 13) & 3) == 0;
#endif
    volatile double tiny = 0x1p-1074, two = 2.0, half = 0.5;
    double g = tiny * two;                      /* gradual underflow: 2^-1073 exactly */
    printf("gradual underflow 2^-1074*2 == 2^-1073: %d\n", g == 0x1p-1073); ok &= g == 0x1p-1073;
    double hl = 0x1p-1022 * half;               /* subnormal result, not flushed */
    printf("2^-1022/2 == 2^-1023 (not flushed): %d\n", hl == 0x1p-1023); ok &= hl == 0x1p-1023;
    volatile double x = 2.0; double s = sqrt(x);
    printf("sqrt(2) == 0x1.6a09e667f3bcdp+0: %d\n", s == 0x1.6a09e667f3bcdp+0); ok &= s == 0x1.6a09e667f3bcdp+0;
    volatile double a = 1.0, b = 0x1p-53; double c = a + b;   /* ties-to-even: 1 + 2^-53 rounds to 1 */
    printf("1 + 2^-53 == 1 (ties to even): %d\n", c == 1.0); ok &= c == 1.0;
    printf("RESULT %s\n", ok ? "OK" : "FAIL");
    return ok ? 0 : 1;
}
