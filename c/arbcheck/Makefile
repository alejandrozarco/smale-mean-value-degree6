# Build the independent arb checker. FLINT >= 3 (with arb/acb), GMP and MPFR are required.
FLINT_PREFIX ?= /opt/homebrew
CC ?= cc
# -fcx-limited-range only affects the double-precision predictor (C99 complex), never a claim.
CFLAGS ?= -O2 -Wall -ffp-contract=off -fcx-limited-range

arbcheck: arbcheck.c
	$(CC) $(CFLAGS) -I$(FLINT_PREFIX)/include -L$(FLINT_PREFIX)/lib -o $@ $< -lflint -lgmp -lmpfr -lpthread

selftest: arbcheck
	./arbcheck --selftest 4 && ./arbcheck --selftest 5 && ./arbcheck --selftest 6 && ./arbcheck --selftest 7

clean:
	rm -f arbcheck

.PHONY: selftest clean
