To make a movie type in terminal:
ffmpeg -r 15 -start_number 0 -i C:\myimages\img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4



