#!/bin/bash

rate1=1
rate2=15

ffmpeg -r $rate1 -start_number 0 -i %02d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim_anim.mp4

rm -rf figs/temp_frames/0*

echo "animation complete"

# Generate a simulation based on frames
