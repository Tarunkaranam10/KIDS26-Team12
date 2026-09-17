#!/bin/bash
# Submit the shared-scratch staging job to the biohackathon LSF queue.
set -euo pipefail

REPO="/research/groups/shelagrp/home/esavage/KIDS26/KIDS26-Team12"
mkdir -p "$REPO/logs"

bsub \
  -q biohackathon \
  -J kids26_t12_share \
  -n 2 \
  -R "rusage[mem=4000]" \
  -W 480 \
  -o "$REPO/logs/share_%J.out" \
  -e "$REPO/logs/share_%J.err" \
  bash "$REPO/scripts/share_to_scratch.sh"
