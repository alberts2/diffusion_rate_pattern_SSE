#!/bin/bash

# Paths to source directories
# SRC1="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/25_taxa"
# SRC2="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/50_taxa"
# SRC3="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/100_taxa"
# SRC4="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/200_taxa"
# SRC5="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/400_taxa"
# SRC6="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/800_taxa"
# #
# SRC1="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/25_taxa"
# SRC2="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/50_taxa"
# SRC3="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/100_taxa"
# SRC4="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/200_taxa"
# SRC5="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/400_taxa"
# SRC6="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/800_taxa"
#
SRC1="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/25_taxa"
SRC2="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/50_taxa"
SRC3="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/100_taxa"
SRC4="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/200_taxa"
SRC5="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/400_taxa"
SRC6="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/800_taxa"
#
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/800_taxa"
#
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/800_taxa"
# 
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/800_taxa"
#
# SRC1="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/25_taxa"
# SRC2="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/50_taxa"
# SRC3="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/100_taxa"
# SRC4="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/200_taxa"
# SRC5="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/400_taxa"
# SRC6="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/800_taxa"

# Path to destination directory
DEST="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/no_rejection_states/fast_rates_80percent_miss/combined_fast"
# DEST="/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/combined_fast"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates/combined_fast"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_20percent_miss/combined_fast"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_40percent_miss/combined_fast"
# DEST="/storage/albert/rate_pattern_diffusion_SSE/data/Simulation/fast_rates_80percent_miss/combined_fast"



# Copy files from both source folders to the destination
cp -r "$SRC1"/* "$DEST"/
cp -r "$SRC2"/* "$DEST"/
cp -r "$SRC3"/* "$DEST"/
cp -r "$SRC4"/* "$DEST"/
cp -r "$SRC5"/* "$DEST"/
cp -r "$SRC6"/* "$DEST"/

echo "All files have been copied to $DEST"
