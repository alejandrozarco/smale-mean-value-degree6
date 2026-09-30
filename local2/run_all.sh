#!/bin/bash
# Full reproducible run of the independent local certificate (d = 6).
# Usage:  bash run_all.sh            (writes run_T1_19.log next to this script)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1   # single-threaded BLAS
LOG=run_T1_19.log
{
  echo "# run started $(date)   python-flint $(python3 -c 'import flint; print(flint.__version__)')"
  start=$(date +%s)
  status=0
  echo; echo "########## main certificate: T = 1/19, Taylor degree 6"
  nice -n 15 python3 certify.py 1/19 6 || status=1
  echo; echo "########## hand-over (SPEC section 3)"
  nice -n 15 python3 handover.py || status=1
  echo; echo "########## extended certificate: T = 1/8, Taylor degree 6 (not needed for the hand-over)"
  nice -n 15 python3 certify.py 1/8 6 || status=1
  echo; echo "########## non-rigorous float cross-check of the pipeline (not part of the proof)"
  nice -n 15 python3 crosscheck.py 1/19 6 || status=1
  end=$(date +%s)
  echo; echo "# overall status: $([ $status -eq 0 ] && echo ALL PASSED || echo FAILURE)   wall time $((end-start)) s"
} 2>&1 | tee "$LOG"
