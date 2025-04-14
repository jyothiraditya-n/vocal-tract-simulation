#! /bin/bash

ffmpeg -y -framerate "24" -i "$1/%05d.png" -an -c:v libx264 \
	-pix_fmt yuv420p -- "$1.mp4"
