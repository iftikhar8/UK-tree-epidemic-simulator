#!/bin/bash

rate1=15
rate2=30

ffmpeg -r $rate1 -start_number 0 -i figs/temp_frames/%04d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim_anim.mp4

rm -rf figs/temp_frames/0*

echo "animation complete"

# Generate a simulation based on frames
