#!/bin/bash
# Stage-0.6 follow-up (pre-registered in results/stage06_dwell_physics/11_followup_preregistration.json).
# Waits for the GPU to be free of other users' jobs before starting (the CUDA kernel sizes its
# buffers to ~95% of free memory, so it must not start next to another GPU job).
cd "$(dirname "$0")"
R=../../../results/stage06_dwell_physics
L=$R/logs_followup.log
while true; do
  used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
  [ "$used" -lt 10000 ] && break
  echo "follow-up (111 jobs) queued; waiting for GPU: ${used} MiB in use by another process ($(date +%H:%M:%S))" > $R/progress_stage06.txt
  sleep 60
done
START=$(date +%Y-%m-%dT%H:%M:%S)
python3 make_followup_plan.py
nohup python3 progress06.py "$START" plan_followup.json > $R/logs_progress_followup.log 2>&1 &
{
echo "GPU free at $START"
python3 fu_cw_steps.py
python3 run06.py --design dwell --reals 0-3 --powers 3,0 --cw-carrier --tag fu_cwcarrier --traces 0
python3 run06.py --design dwell --reals 0-3 --powers 3 --n-samples 2400000 --xis 1,4,8,16,32 --tag fu_long_qpsk --traces 0
python3 run06.py --design dwell --reals 0-3 --powers 3 --n-samples 2400000 --xis 1,4,8,16,32 --cw-carrier --tag fu_long_cw --traces 0
echo FOLLOWUP_DONE
} 2>&1 | grep --line-buffered -v "+++>" > $L
