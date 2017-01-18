import os, glob, sys

def make_movie(dir, imgs, nFrames = None, frameRate=10, output = 'out.mp4'):

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

		for img in imgs:
			imgf.append('tmp/' + img + numbers[i] + '.png')

		print 'Processing image ' + numbers[i]
		
		imgout = 'tmp/imgout{}.png'.format(numbers[i])
		os.system(""" convert -pointsize 40 -fill white -draw 'text 5,100 "Frame: {}" text 5, 10 " No friction" ' {} {} """.format(i, imgf[0], imgout))

		os.system(""" convert -pointsize 40 -fill white +append -draw 'text 5,100 "Frame: {}" ' {} {} {} """.format(i, ' '.join(imgf[1:]), imgout, imgout))
		
	os.system('ffmpeg -framerate {} -i tmp/imgout%0{}d.png {}'.format(frameRate, len(digits), output))

	os.system('rm -r tmp')

make_movie(dir='movie2', imgs = ['friction', 'nofriction', 'stages'], nFrames = 5)
