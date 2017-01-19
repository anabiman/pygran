import os, glob, sys

def make_movie(dir, imgs, imgtxtloc = None, imgtxt = None, nFrames = None, frameRate=10, output = 'out.mp4'):

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

	files = []
	count = 0

	for file in imgFiles:
		files.append(sorted(file, key=lambda item: int( (item.split(dir + '/' + imgs[count])[1]).split('.png')[0] )))
		count += 1

	numbers = []
	nFiles = min([len(file) for file in files])
	digits = str(nFiles)

	for i in range(nFiles):
		num = str(i)
		numbers.append('0' * ( len(digits) - len(num) ) + num)

	for i in range(nFiles):
		imgf = []

		for c,img in enumerate(imgs):
			imgf.append('tmp/' + img + numbers[i] + '.png')
			os.system(""" convert -pointsize 40 -fill white -draw 'text {} "{}" ' {} {} """.format(imgtxtloc[c], imgtxt[c], imgf[-1], imgf[-1]))

		print 'Processing image ' + numbers[i]
		
		imgout = 'tmp/imgout{}.png'.format(numbers[i])
		os.system(""" convert -pointsize 40 -fill white +append -draw 'text 10,80 "Frame: {}" ' {} {} """.format(i, ' '.join(imgf), imgout))
		
	os.system('ffmpeg -framerate {} -i tmp/imgout%0{}d.png {}'.format(frameRate, len(digits), output))

	os.system('rm -r tmp')

# mixed up friction vs nofriction images with ovito ~ oops!
make_movie(dir='movie', imgs = ['friction', 'nofriction', 'stages'], \
	imgtxt = ['no friction', 'friction', 'friction (3 stages)'], imgtxtloc = ['150,990', '150, 990', '150, 990'])
