
movie:
	glslViewer -w 1280 -h 720 -E sequence,0,10 --headless shader.frag
	ffmpeg -r:v 30 -i "%05d.png" -codec:v libx264 -preset veryslow -vb 40M -pix_fmt yuv444p10le -crf 1 -an "movie.mp4"
	rm *.png
