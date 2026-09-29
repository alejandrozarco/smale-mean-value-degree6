#!/bin/bash
# Build c/smale_bb_v2.c + CORE-MATH cr_log and record provenance.
# usage: c/build_v2.sh OUTDIR      -> OUTDIR/smale_bb_v2.bin, OUTDIR/BUILD.txt
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:?usage: c/build_v2.sh OUTDIR}
mkdir -p "$OUT"
SRC=c/smale_bb_v2.c
CM=third_party/core-math-log
SRC_SHA=$(shasum -a 256 "$SRC" | cut -d' ' -f1)
CFLAGS="-O2 -ffp-contract=off -fno-fast-math -std=gnu11 -Wall -Wno-unused-function"
CMD="cc $CFLAGS -DSRC_SHA256_RAW=$SRC_SHA -I$CM -o $OUT/smale_bb_v2.bin $SRC $CM/log.c -lpthread"
eval "$CMD"
BIN_SHA=$(shasum -a 256 "$OUT/smale_bb_v2.bin" | cut -d' ' -f1)
NFMA=$(objdump -d "$OUT/smale_bb_v2.bin" | awk '/>:$/{fn=$NF} /\tf(n?m(add|sub)|ml[as])/{print fn}' | sort | uniq -c | tr '\n' ' ')
{
  echo "# smale_bb_v2 build record"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(uname -a)"
  echo "compiler: $(cc --version | head -1)"
  cc --version | sed 's/^/  /'
  echo "build command (run from repo root):"
  echo "  $CMD"
  echo "source sha256: $SRC_SHA  $SRC"
  echo "core-math sha256:"
  shasum -a 256 $CM/log.c $CM/dint.h | sed 's/^/  /'
  echo "core-math upstream commit: $(cat $CM/UPSTREAM_COMMIT)"
  echo "binary sha256: $BIN_SHA  $OUT/smale_bb_v2.bin"
  echo "fused multiply-add instructions by function: $NFMA (expected: only <_cr_log>, CORE-MATH explicit __builtin_fma; -ffp-contract=off)"
  echo "sha256 self-test: $("$OUT/smale_bb_v2.bin" sha256 "$SRC" | cut -d' ' -f1) (must equal source sha256)"
} > "$OUT/BUILD.txt"
cat "$OUT/BUILD.txt"
