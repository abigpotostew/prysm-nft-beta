DURATION ?= 10
WIDTH ?= 1280
HEIGHT ?= 720
OUTPUT ?= movie.mp4

movie:
	glslViewer -w $(WIDTH) -h $(HEIGHT) -E sequence,0,$(DURATION) --headless shader.frag
	ffmpeg -r:v 30 -i "%05d.png" -codec:v libx264 -preset veryslow -vb 40M -pix_fmt yuv444p10le -crf 1 -an "$(OUTPUT)"
	rm *.png
