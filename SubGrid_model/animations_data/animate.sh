#!/bin/bash

# Generate a simulation based on frames
rm -rf /raw_data/.DS_Store

python3 animate.py
rate1=15
rate2=30

ffmpeg -r $rate1 -start_number 0 -i temp_frames/%04d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim_anim.mp4

rm -rf temp_frames/0*  % remove temp frames
rm -rf raw_data/0* % remove raw data
echo "animation complete:"
echo "REMOVE UNWANTED: raw data DATA & temp_frames!!!!"


