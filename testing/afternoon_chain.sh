#!/bin/bash
# Waits for the 9-camera rerun (pid $1) to finish, then runs the 15-CCD and
# 30-CCD random campaigns back to back. Detached so it survives independent
# of the driving Claude Code session across the whole node's walltime window.
set -x
RERUN_PID=$1
cd /global/cfs/cdirs/desicollab/users/cdwarner/code/specex

while kill -0 $RERUN_PID 2>/dev/null; do sleep 5; done
echo "9-camera rerun (pid $RERUN_PID) done, starting 15-CCD random campaign at $(date)"

mkdir -p /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_15
python testing/full_ccd_campaign.py \
  --cases-file /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_15/all_cases.jsonl \
  --outdir /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_15 \
  > /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_15/campaign.log 2>&1

echo "15-CCD campaign done at $(date), starting 30-CCD random campaign"

mkdir -p /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_30
python testing/full_ccd_campaign.py \
  --cases-file /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_30/all_cases.jsonl \
  --outdir /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_30 \
  > /pscratch/sd/c/cdwarner/specex/testing/random_full_ccd_30/campaign.log 2>&1

echo "30-CCD campaign done at $(date). CHAIN_COMPLETE"
