DURATION ?= 10
WIDTH ?= 1280
HEIGHT ?= 720
OUTPUT ?= movie.mp4
COMPATIBLE ?= false
PIX_FMT ?= yuv444p10le

ifeq ($(COMPATIBLE),true)
	PIX_FMT = yuv420p
endif


test:
	echo $(PIX_FMT)

# yuv420p is compatible but lossy

movie:
	glslViewer -w $(WIDTH) -h $(HEIGHT) -E sequence,0,$(DURATION) --headless shader.frag
	ffmpeg -r:v 30 -framerate 30 -i "%05d.png" -codec:v libx264 -preset veryslow -vb 40M -pix_fmt $(PIX_FMT) -crf 1 -an "$(OUTPUT)"
	rm *.png
