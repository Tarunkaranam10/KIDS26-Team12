#!/usr/bin/env bash
# =============================================================================
# scripts/watch_loco_array.sh - Poll the LOCO array, resubmit failures, merge.
# =============================================================================
#
# USAGE
#   scripts/watch_loco_array.sh <JOBID> [RUN_DIR] [POLL_SECONDS]
#   e.g.  nohup scripts/watch_loco_array.sh 323078995 results/loco_run01 300 \
#           >> logs/watchdog_323078995.log 2>&1 &
#
# WHY THIS EXISTS
# ---------------
# The array is 30 folds at ~55 min each, 4 concurrent => ~8 h wall. Nobody
# should sit on it. This polls, resubmits any index that EXITs, and runs the
# merge exactly once when all 30 fold outputs are on disk.
#
# WHAT IT DELIBERATELY DOES NOT DO
# --------------------------------
# It never interprets results and never relaxes a gate. loco_merge.R refuses
# partial cohorts; this script does not work around that, it waits for the
# cohort to be complete. If a fold fails repeatedly it stops and says so rather
# than merging 29 folds and quietly calling it a run.
#
# Resubmission spec mirrors scripts/lsf_loco_array.bsub exactly:
#   -n 4, rusage[mem=60GB] PER SLOT (=240GB/task), -M 240GB, -W 8:00.
# rusage[mem] is PER SLOT on this cluster - see commit 5cc611d. Do not "fix"
# the 60GB to 240GB here; that reserves 960GB/task and pends forever.
# =============================================================================
set -uo pipefail

JOBID="${1:?usage: watch_loco_array.sh <JOBID> [RUN_DIR] [POLL_SECONDS]}"
RUN_DIR="${2:-results/loco_run01}"
POLL="${3:-300}"
FOLD_DIR="$RUN_DIR/folds"
MASTER="data/processed/master_samples.tsv"
BETA="data/processed/beta.tsv"
N_FOLDS=30
MAX_RETRIES=2

declare -A RETRIES
log() { echo "[$(date '+%F %T')] $*"; }

# A fold counts as complete only if its metrics file exists AND is non-empty.
# A truncated file from a killed task must not be mistaken for a finished fold.
completed_folds() {
  find "$FOLD_DIR" -name 'metrics_*.tsv' -size +0c 2>/dev/null | wc -l
}

# Indices LSF reports as EXIT for this array.
failed_indices() {
  bjobs -a -noheader -o "jobindex stat" "$JOBID" 2>/dev/null \
    | awk '$2=="EXIT"{print $1}'
}

pending_or_running() {
  bjobs -noheader -o "stat" "$JOBID" 2>/dev/null \
    | grep -cE 'PEND|RUN'
}

resubmit() {
  local idx="$1"
  local n="${RETRIES[$idx]:-0}"
  if [ "$n" -ge "$MAX_RETRIES" ]; then
    log "fold $idx has failed $n times; NOT resubmitting. Inspect logs/loco_${JOBID}_${idx}.err"
    return 1
  fi
  RETRIES[$idx]=$((n+1))
  log "resubmitting fold $idx (attempt $((n+1))/$MAX_RETRIES)"
  bsub -q biohackathon -n 4 -R "span[hosts=1] rusage[mem=60GB]" -M 240GB -W 8:00 \
       -J "kids26_loco_rerun_${idx}" \
       -o "logs/loco_rerun_${idx}_%J.out" -e "logs/loco_rerun_${idx}_%J.err" \
       "cd $PWD && Rscript scripts/loco_one_fold.R $BETA $MASTER $RUN_DIR $idx"
}

log "watching array $JOBID -> $FOLD_DIR (poll ${POLL}s, expecting $N_FOLDS folds)"

while true; do
  done_n=$(completed_folds)
  active=$(pending_or_running)

  if [ "$done_n" -ge "$N_FOLDS" ]; then
    log "all $N_FOLDS folds present; running merge"
    Rscript scripts/loco_merge.R "$RUN_DIR" "$MASTER" \
      && log "MERGE COMPLETE -> inspect pooled skill_vs_tissue_mean and lambda_at_boundary in $RUN_DIR" \
      || log "MERGE FAILED - see output above"
    exit 0
  fi

  # Chase any index LSF reports as EXIT. resubmit() caps attempts per index.
  stuck=0
  for idx in $(failed_indices); do
    resubmit "$idx" || stuck=1
  done

  if [ "$active" -eq 0 ] && [ "$done_n" -lt "$N_FOLDS" ] && [ "$stuck" -eq 1 ]; then
    log "STOPPING: nothing in flight, only $done_n/$N_FOLDS folds done, and a fold exceeded retries."
    log "The cohort is incomplete - do NOT merge. loco_merge.R will refuse anyway."
    exit 1
  fi

  log "progress: $done_n/$N_FOLDS folds complete, $active tasks pend/run"
  sleep "$POLL"
done
