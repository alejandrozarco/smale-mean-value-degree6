/* flint_compat.h: force-included by the Makefile (-include), so arbcheck.c stays byte-identical to the checked source.
 * FLINT 3.0 names the random state functions flint_randinit / flint_randclear; FLINT 3.1 renamed them to
 * flint_rand_init / flint_rand_clear (used in arbcheck.c, self-test only). FLINT 2.x is not supported. */
#include <flint/flint.h>
#if !defined(__FLINT_RELEASE) || __FLINT_RELEASE < 30000
#error "arbcheck needs FLINT >= 3.0 (arb is part of FLINT from 3.0 on)"
#endif
#if __FLINT_RELEASE < 30100
#define flint_rand_init flint_randinit
#define flint_rand_clear flint_randclear
#endif
