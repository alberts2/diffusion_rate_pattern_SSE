#!/bin/bash

# Number of batches
NBATCHES=10
BATCHSIZE=10
MISSING=20

# Path to your R script with sim_bisse_slow and sim_bisse_fast defined
# RSCRIPT="sim_bisse_missing.R"
RSCRIPT="sim_bisse_missing_no_rejection.R"

# Path to store logs
LOGDIR="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/Log"
# LOGDIR="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/Log"
mkdir -p "$LOGDIR"

# Loop over batches
for b in $(seq 1 $NBATCHES)
do
  echo "Starting batch $b"
  
 Rscript $RSCRIPT $BATCHSIZE $b $MISSING > "$LOGDIR/batch_${b}.log" 2>&1 &
  
  # optional: small delay to avoid starting all jobs at the exact same time
  sleep 1
done

# Wait for all background jobs to finish
wait
echo "All batches finished!"
