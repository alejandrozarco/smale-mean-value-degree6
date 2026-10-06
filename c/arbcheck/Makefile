# Build the independent arb checker. FLINT >= 3.0 (with arb/acb), GMP and MPFR are required.
# Linux: FLINT_PREFIX=/usr (e.g. Ubuntu 24.04: apt install libflint-dev libgmp-dev libmpfr-dev, FLINT 3.0.1).
FLINT_PREFIX ?= /opt/homebrew
CC ?= cc
# -fcx-limited-range only affects the double-precision predictor (C99 complex), never a claim.
CFLAGS ?= -O2 -Wall -ffp-contract=off -fcx-limited-range

# GCC 13 prints -Wstringop-overflow/-overread warnings inside inlined FLINT functions (arb_add_error_mag, mag_addmul,
# arb_log). The arguments are ordinary arb_t/mag_t struct fields, clang gives no such warning, and the d = 5 check on
# Linux (GCC 13.3, FLINT 3.0.1) gives the same counts and PASS; we take them to be false positives of GCC's size analysis.
# flint_compat.h is force-included: FLINT 3.0 compatibility and a clear error for FLINT 2.x. arbcheck.c is unchanged.
arbcheck: arbcheck.c flint_compat.h
	$(CC) $(CFLAGS) -include flint_compat.h -I$(FLINT_PREFIX)/include -L$(FLINT_PREFIX)/lib -o $@ $< -lflint -lgmp -lmpfr -lpthread -lm

# arbcheck2: arbcheck + the w-form of SPEC_W.md (README2.md). Same flags and force-included flint_compat.h.
arbcheck2: arbcheck2.c flint_compat.h
	$(CC) $(CFLAGS) -include flint_compat.h -I$(FLINT_PREFIX)/include -L$(FLINT_PREFIX)/lib -o $@ $< -lflint -lgmp -lmpfr -lpthread -lm

selftest2: arbcheck2
	./arbcheck2 --selftest 4 && ./arbcheck2 --selftest 5 && ./arbcheck2 --selftest 6 && ./arbcheck2 --selftest 7

selftest: arbcheck
	./arbcheck --selftest 4 && ./arbcheck --selftest 5 && ./arbcheck --selftest 6 && ./arbcheck --selftest 7

clean:
	rm -f arbcheck arbcheck2

.PHONY: selftest selftest2 clean
