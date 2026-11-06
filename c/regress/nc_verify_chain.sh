#!/bin/bash
# Negative controls for the bound verification chain (c/verify_chunks.py, c/cover_check2.py) on the d = 5 v3 run, plus
# two d = 7 records end to end. Self-checking: prints PASS/WRONG per control, exit 1 if any control behaves wrongly.
# usage: c/regress/nc_verify_chain.sh WORKDIR [--no-d7]     (log it to runs/regress_v3/nc_verify_chain.log)
# Needs c/arbcheck/arbcheck2, runs/d5_v3/, runs/d6_v3/d6.tasks and, for the d = 7 part, the d = 7 tasks file and the
# tree of machine A (env D7TASKS, D7TREE; defaults: this repository's run records, or runs/d7/d7.tasks and A_tree.zst).
set -uo pipefail
W=$(cd "$(dirname "$0")/../.." && pwd); WD=${1:?usage: nc_verify_chain.sh WORKDIR [--no-d7]}; D7=${2:-}
mkdir -p "$WD"; WD=$(cd "$WD" && pwd); cd "$WD" || exit 1
if command -v shasum >/dev/null 2>&1; then SHA="shasum -a 256"; else SHA=sha256sum; fi
sha() { $SHA "$1" | cut -d' ' -f1; }
T=$W/runs/d5_v3/d5.tasks; TS=$(sha $T); AB=$W/c/arbcheck/arbcheck2; ABS=$(sha $AB)
CFG=$(head -1 $T | sed 's/.*cfg=\([0-9a-f]*\).*/\1/'); P=$W/runs/d5_v3/d5; CH=200000
VC="nice -n 19 python3 -I $W/c/verify_chunks.py"; CC="python3 -I $W/c/cover_check2.py"; NT="python3 -I $W/c/regress/nc_tools.py"
wrong=0; n=0
check() { n=$((n+1)); if [ "$2" = "$3" ]; then echo "PASS  $1"; else echo "WRONG $1 (expected $2, got $3)"; wrong=$((wrong+1)); fi; }
vrun() { $VC "$1" "${5:-$T}" "${6:-$TS}" $CFG $AB $ABS 2 $CH "$3" "$4" 5 < "$2" >> nc.out 2>&1; echo $?; }  # label stream sha wd [tasks sha]
cover() { $CC "$1" "${ABX:-$ABS}" $CFG "$2" "${@:3}" >> nc.out 2>&1; echo $?; }                            # prefix streams dirs
streams() { echo "$1 $2" > "$3"; }

echo "chain: verify_chunks $(sha $W/c/verify_chunks.py), cover_check2 $(sha $W/c/cover_check2.py), arbcheck2 $ABS"
echo "date: $(date -u +%FT%TZ), host: $(uname -sm)"
if [ ! -f good.tree ]; then
  cc -O2 -ffp-contract=off -fno-fast-math -std=gnu11 -w -DSRC_SHA256_RAW=$(sha $W/c/smale_bb_v3.c) -I$W/c -I$W/third_party/core-math-log \
     -o export_tree_v3 $W/c/export_tree_v3.c $W/third_party/core-math-log/log.c -lpthread -lm || exit 1
  nice -n 19 ./export_tree_v3 5 0.05 $T good 2 > /dev/null 2>&1 || exit 1
fi
GOOD=good.tree; GS=$(sha $GOOD); streams c0 $GS s0
BIG=$($NT largest $GOOD | cut -d' ' -f1); X=$($NT firstx $GOOD)
rm -rf w? w?? ; 

check "C0 positive: correct inputs verify" 0 $(vrun c0 $GOOD $GS w0)
check "C0 positive: COVER COMPLETE" 0 $(cover $P s0 w0)

$NT rootleaf $GOOD t1.tree $BIG; s=$(sha t1.tree); streams c1 $s s1
check "NC1 record replaced by one F leaf: verify rejects" 2 $(vrun c1 t1.tree $s w1)
check "NC1 cover rejects" 1 $(cover $P s1 w1)
check "NC1 failed receipt written" 1 $(ls w1 | grep -c 'FAILED$')

$NT rootleaf $GOOD t1x.tree $X; s=$(sha t1x.tree)
check "NC1x record near an equality point replaced by one F leaf: verify rejects" 2 $(vrun c1x t1x.tree $s w1x)

$NT swapcols $T t2.tasks 0 1; mkdir -p cd2; cp $P.done cd2/d5.done; cp t2.tasks cd2/d5.tasks; : > cd2/d5.unresolved; T2S=$(sha t2.tasks); streams c2 $GS s2
check "NC2 non-canonical tasks with a consistent hash: verify cannot know (passes)" 0 $(vrun c2 $GOOD $GS w2 t2.tasks $T2S)
check "NC2 cover rejects (check_done on the same file)" 1 $(cover cd2/d5 s2 w2 )

cp -r w0 w3; s=$(sha t1.tree); streams c0 $s s3
check "NC3 tampered stream into a finished workdir: verify refuses" 2 $(vrun c0 t1.tree $s w3)
check "NC3 stale receipts, expected hash of the new stream: cover rejects" 1 $(cover $P s3 w3)
check "NC3 stale receipts, expected hash of the old stream: cover accepts (genuine)" 0 $(cover $P s0 w3)

$NT drop $GOOD t4.tree $BIG; s=$(sha t4.tree); streams c4 $s s4
check "NC4 missing task: verify passes the rest" 0 $(vrun c4 t4.tree $s w4)
check "NC4 cover rejects" 1 $(cover $P s4 w4)

$NT dup $GOOD t5.tree 7; s=$(sha t5.tree); streams c5 $s s5
check "NC5 duplicate record: cover rejects" 1 $(vrun c5 t5.tree $s w5 >/dev/null; cover $P s5 w5)

$NT unres $GOOD t6.tree 7; s=$(sha t6.tree); streams c6 $s s6
check "NC6 unresolved leaf: verify rejects" 2 $(vrun c6 t6.tree $s w6)
check "NC6 cover rejects" 1 $(cover $P s6 w6)

check "NC7 wrong expected stream hash: verify rejects" 2 $(vrun c7 $GOOD 0000000000000000000000000000000000000000000000000000000000000000 w7)
streams c7 $GS s7
check "NC7 cover rejects" 1 $(cover $P s7 w7)

check "NC8 tasks file of another degree: verify rejects" 2 $(vrun c8 $GOOD $GS w8 $W/runs/d6_v3/d6.tasks $(sha $W/runs/d6_v3/d6.tasks))

$NT tear $GOOD t10.tree; s=$(sha t10.tree); streams c10 $s s10
check "NC10 torn tail: verify rejects" 2 $(vrun c10 t10.tree $s w10)
check "NC10 cover rejects" 1 $(cover $P s10 w10)

mkdir w11; head -c 100000 $GOOD > w11/c11.0000.tree; echo junk > w11/c11.0000.tasks; echo junk > w11/c11.0000.ids.tmp; streams c11 $GS s11
check "NC11 leftovers of a killed run, no receipt: recomputed" 0 $(vrun c11 $GOOD $GS w11)
check "NC11 cover accepts" 0 $(cover $P s11 w11)

cp -r w0 w12; $NT editids w0/c0.0000.ids w12/c0.0000.ids
check "NC12 ids file edited afterwards: cover rejects" 1 $(cover $P s0 w12)

check "NC13 receipts from another arbcheck2 hash: cover rejects" 1 $(ABX=deadbeef cover $P s0 w0)

echo "other 0000000000000000000000000000000000000000000000000000000000000000" | cat s0 - > s14
check "NC14 stream in STREAMS without a receipt: cover rejects" 1 $(cover $P s14 w0)
echo "other 0" > s15
check "NC15 receipt for a stream not in STREAMS: cover rejects" 1 $(cover $P s15 w0)

cp -r w0 w16; for f in w16/*.ids; do zstd -q --rm $f; done
check "NC16 compressed ids files: cover accepts" 0 $(cover $P s0 w16)
cp -r w0 w17; $NT editids w0/c0.0000.ids w17/c0.0000.ids; for f in w17/*.ids; do zstd -q --rm $f; done
check "NC17 edited compressed ids file: cover rejects" 1 $(cover $P s0 w17)

CT=${D7TASKS:-$W/runs/d7_v3/merge/d7.tasks}; [ -f $CT ] || CT=$W/runs/d7/d7.tasks
TR=${D7TREE:-$W/runs/d7_v3/d7.tree.zst}; [ -f $TR ] || TR=$W/A_tree.zst
if [ "$D7" != "--no-d7" ] && { [ ! -f $CT ] || [ ! -f $TR ]; }; then echo "SKIP  d = 7 part (set D7TASKS and D7TREE, e.g. A_tree.zst of release v1.4)"; D7=--no-d7; fi
if [ "$D7" != "--no-d7" ]; then
  CTS=$(sha $CT); CF7=$(head -1 $CT | sed 's/.*cfg=\([0-9a-f]*\).*/\1/')
  zstd -dc $TR 2>/dev/null | $NT pickx d7.tree 2 400000000
  ID=$(python3 -I $W/c/regress/nc_tools.py firstx d7.tree); $NT rootleaf d7.tree d7bad.tree $ID
  v7() { $VC "$1" $CT $CTS $CF7 $AB $ABS 2 $CH "$(sha $2)" "$3" 7 < "$2" >> nc.out 2>&1; echo $?; }
  check "D7 two real records near equality points verify" 0 $(v7 d7ok d7.tree w71)
  check "D7 one of them replaced by one F leaf: verify rejects" 2 $(v7 d7bad d7bad.tree w72)
fi
echo "$n controls, $wrong wrong"
[ $wrong = 0 ]
