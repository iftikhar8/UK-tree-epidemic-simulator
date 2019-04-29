#!/bin/bash

rate1=15
rate2=30

# convert .npy files to matplotlib generated .png's
python3 animate.py

# convert images to MP4 file
ffmpeg -r $rate1 -start_number 0 -i frames_2_anim/img-%05d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim_anim.mp4

# remove temp data
rm -rf frames_2_anim/0*
rm -rf dat_2_anim/0*

echo "animation complete"
