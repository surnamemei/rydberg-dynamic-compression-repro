#!/bin/bash
# Stage-0.6 decisive campaign (sequential: one GPU).
cd "$(dirname "$0")"
L=../../../results/stage06_dwell_physics/logs_campaign.log
{
python3 stage05_recheck.py
python3 run06.py --design dwell   --reals 0-7 --powers 3,0,6 --traces 0
python3 run06.py --design shuffle --reals 0-7 --powers 3,0,6 --traces 0
python3 run06.py --design dwell   --reals 0-7 --powers -6 --traces 0
python3 run06.py --design shuffle --reals 0-7 --powers -6 --traces 0
# Step 14: numerical stability on decisive cases (realizations 0,1; +3 dB)
python3 run06.py --design dwell   --reals 0-1 --powers 3 --variants xi_0.05,xi_0.25,xi_1,xi_4 --nd 8001 --tag num_Nd8001
python3 run06.py --design dwell   --reals 0-1 --powers 3 --variants xi_0.05,xi_0.25,xi_1,xi_4 --dt 5e-10 --tag num_dt0p5
python3 run06.py --design shuffle --reals 0-1 --powers 3 --variants ORIGINAL,BLOCK_SHUFFLED_CYCLES --nd 8001 --tag num_Nd8001
python3 run06.py --design shuffle --reals 0-1 --powers 3 --variants ORIGINAL,BLOCK_SHUFFLED_CYCLES --dt 5e-10 --tag num_dt0p5
echo CAMPAIGN_DONE
} 2>&1 | grep --line-buffered -v "+++>" > $L
