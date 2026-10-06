#!/bin/bash
# d = 7 runs of the local certificate (n = 6, c = 6/7, omega = exp(2 pi i/6)); see README.md, section "d = 7".
# Usage:  bash run_d7.sh     (writes run_d7_T1_19.log, run_d7_T1_12.log, run_d7_T1_10.log, run_d7_T1_8.log, run_d7_T1_5.log)
# Every job runs at nice 15 on one thread; the longest single job (N = 8) takes about 1 minute and 0.3 GB.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
P="nice -n 15 python3"
overall=0
run() {   # run <log> <commands...>: each command is one string
  local LOG=$1; shift
  {
    echo "# run started $(date)   python-flint $(python3 -c 'import flint; print(flint.__version__)')"
    local start=$(date +%s) status=0
    for c in "$@"; do
      echo; echo "########## $c"
      $P $c || status=1
    done
    local end=$(date +%s)
    echo; echo "# overall status: $([ $status -eq 0 ] && echo ALL PASSED || echo FAILURE)   wall time $((end-start)) s"
    return $status
  } 2>&1 | tee "$LOG"
  local ps=("${PIPESTATUS[@]}")
  [ "${ps[0]}" -eq 0 ] && [ "${ps[1]}" -eq 0 ] || overall=1
}
run run_d7_T1_19.log "exact_qomega.py --d=7" "certify.py --d=7 --joint 1/19 6" "handover.py --d=7 1/19 1/20" \
    "handover.py --d=7 --rmax 1/19" "crosscheck.py --d=7 --joint 1/19 6"
run run_d7_T1_12.log "certify.py --d=7 --joint 1/12 6" "handover.py --d=7 --rmax 1/12"
run run_d7_T1_10.log "certify.py --d=7 --joint 1/10 6" "handover.py --d=7 --rmax 1/10"
run run_d7_T1_8.log  "certify.py --d=7 --joint 1/8 6" "certify.py --d=7 --joint 1/8 8" "handover.py --d=7 --rmax 1/8" \
    "crosscheck.py --d=7 --joint 1/8 6"
run run_d7_T1_5.log  "certify.py --d=7 --joint 1/6 8" "certify.py --d=7 --joint 1/5 8" "handover.py --d=7 --rmax 1/6 1/5"
echo "run_d7.sh: $([ $overall -eq 0 ] && echo ALL PASSED || echo FAILURE)"
exit $overall
