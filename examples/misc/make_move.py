import os, glob, sys, multiprocessing
from functools import partial

def processImgs(imgs, numbers, imgFiles, imgtxt, imgtxtloc, i):

	imgf = []

	for c,img in enumerate(imgs):
		imgf.append('tmp/' + img + numbers[i] + '.png')

		if imgtxt:
			os.system(""" convert -pointsize 80 -fill black -draw 'text {} "{}" ' {} {} """.format(imgtxtloc[c], imgtxt[c], imgFiles[c][i], imgf[-1]))
		else:
			os.system(""" cp {} {} """.format(imgFiles[c][i], imgf[-1]))

	print 'Processing image ' + numbers[i]
	
	imgout = 'tmp/imgout{}.png'.format(numbers[i])

	os.system(""" convert -pointsize 80 -fill black +append {} {} """.format(' '.join(imgf), imgout))

def make_movie(dir, imgs, imgtxtloc = None, imgtxt = None, nFrames = None, frameRate=10, output = 'out.mp4', skip = None, startfrom = 0, parallel = 0):

	if not os.path.exists('tmp'):
		os.system('mkdir tmp')
	else:
		print 'path tmp already exists'
		sys.exit()

	for img in imgs:
		os.system('cp {}/{}*.png tmp/'.format(dir, img))

	imgFiles = [glob.glob("{}/{}*.png".format(dir, img)) for img in imgs]

	for i in range(len(imgs)):
		imgFiles[i] = imgFiles[i][:nFrames]

	if skip:
		for i in range(len(imgs)):
			imgFiles[i] = imgFiles[i][skip]

	count = 0

	for i in range(len(imgFiles)):
		imgFiles[i].sort(key=os.path.getmtime)

	numbers = []
	nFiles = min([len(file) for file in imgFiles])
	digits = str(nFiles)

	for i in range(startfrom, nFiles + startfrom):
		num = str(i)
		numbers.append('0' * ( len(digits) - len(num) ) + num)

	if parallel:
		pool = multiprocessing.Pool(processes=parallel)
		func = partial(processImgs, imgs, numbers, imgFiles, imgtxt, imgtxtloc)
		pool.map(func, range(nFiles))
		pool.close()
		pool.join()
	else:
		for i in range(nFiles):
			processImgs(i, imgs, numbers, imgFiles, imgtxt, imgtxtloc)

	os.system('ffmpeg -start_number {} -framerate {} -i tmp/imgout%0{}d.png -q:v 4 {}'.format(startfrom, frameRate, len(digits), output))
	os.system('rm -r tmp')

#make_movie(dir='movie2', imgs = ['compress'])

make_movie(dir='movie2', imgs = ['compress'], parallel = 8)

# classical example 
# mixed up friction vs nofriction images with ovito ~ oops!
# make_movie(dir='movie', imgs = ['friction', 'nofriction', 'stages'], \
#	imgtxt = ['no friction', 'friction', 'friction (3 stages)'], imgtxtloc = ['150,990', '150, 990', '150, 990'])

# e.g. on using 'startfrom' to specify where image sequence starts from
# make_movie(dir='movie', imgs = ['beverloo-bottom', 'beverloo-side', 'beverloo-core'], \
#	imgtxt = ['bottom view', 'side view', 'core view'], imgtxtloc = ['150,990', '150,990', '150,990'], startfrom = 80)

#import ovito
#import inspect

# The following function is called by OVITO to let the script
# draw arbitrary graphics content into the viewport.
# It is passed a QPainter (see http://qt-project.org/doc/qt-5/qpainter.html).
def render(painter, **args):
	# This demo code prints the current animation frame
	# into the upper left corner of the viewport.
	font = painter.font()
	font.setPointSize(20.0)
	painter.setFont(font)
	
	xpos = 10
	ypos = 10 + painter.fontMetrics().ascent()
	text = "Time: {:.5f} s".format(ovito.dataset.anim.current_frame * 1e-5)
	painter.drawText(xpos, ypos, text)
	# The following code prints the current number of particles
	# into the lower left corner of the viewport.
	
	xpos = 10
	ypos = painter.window().height() - 10
	if ovito.dataset.selected_node:
		num_particles = ovito.dataset.selected_node.compute().number_of_particles
		text = "{} particles".format(num_particles)
	else:
		text = "no particles"
	#painter.drawText(xpos, ypos, text)
	font.setPointSize(15.0)
	painter.setFont(font)
	
	painter.drawText(xpos + 35, ypos - 110, '0.5 mm')

	# resolution : 800 x 1080