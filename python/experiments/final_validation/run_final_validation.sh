#!/bin/bash
# Final validation: (1) replication [controls at LO'=0.35, cases, Nd8001 checks], then (2) regeneration.
# Separate scripts, output directories and logs. GPU-gated so it never starts next to another user's job.
cd "$(dirname "$0")"
FVR=../../../results/final_validation
while [ "$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)" -ge 10000 ]; do
  echo "waiting for GPU $(date +%H:%M:%S)" > $FVR/replication/progress.txt; sleep 60; done
python3 fv_controls.py 0.35 2>&1 | grep --line-buffered -v "+++>" > $FVR/replication/logs_controls_LO0.35.log
python3 fv_replication.py 2>&1 | grep --line-buffered -v "+++>" > $FVR/replication/logs_replication.log
echo REPLICATION_DONE >> $FVR/replication/logs_replication.log
python3 regen_stage05.py 2>&1 | grep --line-buffered -v "+++>" > $FVR/regen_stage05/logs_regen.log
echo REGEN_DONE >> $FVR/regen_stage05/logs_regen.log
