#!/bin/bash
# d7_tree_stream.sh MANIFEST SHALOG -- write the decompressed d = 7 tree streams to standard output, one after the
# other, and append "date label sha256" of each decompressed stream to SHALOG (compare with runs/d7/verify/streams.txt).
# MANIFEST lines: LABEL FILE.zst [FILE.zst ...]  (a stream is the decompressed concatenation of its files, in order).
# Used as FP_TREE_CMD for sample_tree.py at d = 7.
set -o pipefail
while read -r lab files; do
  case "$lab" in ''|\#*) continue ;; esac
  { sha=$(cat $files | zstd -dcq | tee /dev/fd/3 | sha256sum | cut -c1-64); } 3>&1 || { echo "$(date +%FT%T) $lab FAILED" >> "$2"; exit 1; }
  echo "$(date +%FT%T) $lab $sha" >> "$2"
done < "$1"
