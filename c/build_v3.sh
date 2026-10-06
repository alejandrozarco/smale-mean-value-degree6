#!/bin/bash
# Build c/smale_bb_v3.c + CORE-MATH cr_log and record provenance.
# usage: c/build_v3.sh OUTDIR      -> OUTDIR/smale_bb_v3.bin, OUTDIR/BUILD.txt
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:?usage: c/build_v3.sh OUTDIR}
mkdir -p "$OUT"
SRC=c/smale_bb_v3.c
CM=third_party/core-math-log
# sha256: shasum (macOS, most Linux) or sha256sum (coreutils)
if command -v shasum >/dev/null 2>&1; then SHA256="shasum -a 256"; else SHA256="sha256sum"; fi
SRC_SHA=$($SHA256 "$SRC" | cut -d' ' -f1)
CFLAGS="-O2 -ffp-contract=off -fno-fast-math -std=gnu11 -Wall -Wno-unused-function"
CMD="cc $CFLAGS -DSRC_SHA256_RAW=$SRC_SHA -I$CM -o $OUT/smale_bb_v3.bin $SRC $CM/log.c -lpthread -lm"
eval "$CMD"
BIN_SHA=$($SHA256 "$OUT/smale_bb_v3.bin" | cut -d' ' -f1)
# fused multiply-add scan: arm64 fmadd/fmsub/fnmadd/fnmsub/fmla/fmls, x86-64 vfmadd*/vfmsub*/vfnmadd*/vfnmsub*
if command -v objdump >/dev/null 2>&1; then
  NFMA=$(objdump -d "$OUT/smale_bb_v3.bin" | awk '/>:$/{fn=$NF} /\t(f(n?m(add|sub)|ml[as])|vfn?m(add|sub))/{print fn}' | sort | uniq -c | tr '\n' ' ')
else
  NFMA="(objdump not found; scan skipped)"
fi
{
  echo "# smale_bb_v3 build record"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(uname -a)"
  echo "compiler: $(cc --version | head -1)"
  cc --version | sed 's/^/  /'
  echo "build command (run from repo root):"
  echo "  $CMD"
  echo "source sha256: $SRC_SHA  $SRC"
  echo "core-math sha256:"
  $SHA256 $CM/log.c $CM/dint.h | sed 's/^/  /'
  echo "core-math upstream commit: $(cat $CM/UPSTREAM_COMMIT)"
  echo "binary sha256: $BIN_SHA  $OUT/smale_bb_v3.bin"
  echo "fused multiply-add instructions by function: $NFMA (expected: only cr_log, shown as <_cr_log> on macOS and <cr_log> on Linux; CORE-MATH explicit __builtin_fma; -ffp-contract=off)"
  echo "sha256 self-test: $("$OUT/smale_bb_v3.bin" sha256 "$SRC" | cut -d' ' -f1) (must equal source sha256)"
} > "$OUT/BUILD.txt"
cat "$OUT/BUILD.txt"
