#r
!/bin/bash

rate1=1
rate2=1

ffmpeg -r $rate1 -start_number 0 -i temp_figs/%03d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim_anim.mp4

echo "animation complete"

# Generate a simulation based on frames

