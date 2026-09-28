#!/bin/bash
# Controls at converged Nd, CW-only checks, then tau_atom.
set -e
cd "$(dirname "$0")"
L=../../../results/stage06_dwell_physics
python3 build_controls06.py 4001 > $L/logs_controls_Nd4001.log 2>&1
python3 build_controls06.py 8001 --no-kernel > $L/logs_controls_Nd8001.log 2>&1
python3 build_controls06.py 1501 --no-kernel > $L/logs_controls_Nd1501.log 2>&1
python3 tau_atom.py > $L/logs_tau_atom.log 2>&1
