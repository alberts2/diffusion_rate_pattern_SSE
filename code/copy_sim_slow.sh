#!/bin/bash

# Paths to source directories
# SRC1="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/25_taxa"
# SRC2="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/50_taxa"
# SRC3="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/100_taxa"
# SRC4="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/200_taxa"
# SRC5="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/400_taxa"
# SRC6="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/800_taxa"
#
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/800_taxa"
#
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/800_taxa"
# 
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/800_taxa"
#
SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/25_taxa"
SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/50_taxa"
SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/100_taxa"
SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/200_taxa"
SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/400_taxa"
SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/800_taxa"




# Path to destination directory
#DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates/combined_slow"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_20percent_miss/combined_slow"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_40percent_miss/combined_slow"
DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/slow_rates_80percent_miss/combined_slow"

# Copy files from both source folders to the destination
cp -r "$SRC1"/* "$DEST"/
cp -r "$SRC2"/* "$DEST"/
cp -r "$SRC3"/* "$DEST"/
cp -r "$SRC4"/* "$DEST"/
cp -r "$SRC5"/* "$DEST"/
cp -r "$SRC6"/* "$DEST"/

echo "All files have been copied to $DEST"
