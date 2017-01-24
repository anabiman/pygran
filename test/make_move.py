import os, glob, sys

def make_movie(dir, imgs, imgtxtloc = None, imgtxt = None, nFrames = None, frameRate=10, output = 'out.mp4', skip = None, startfrom = 0):

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

	files = []
	count = 0

	for file in imgFiles:
		files.append(sorted(file, key=lambda item: int( (item.split(dir + '/' + imgs[count])[1]).split('.png')[0] )))
		count += 1

	numbers = []
	nFiles = min([len(file) for file in files])
	digits = str(nFiles)

	for i in range(startfrom, nFiles + startfrom):
		num = str(i)
		numbers.append('0' * ( len(digits) - len(num) ) + num)

	for i in range(nFiles):
		imgf = []

		for c,img in enumerate(imgs):
			imgf.append('tmp/' + img + numbers[i] + '.png')
			os.system(""" convert -pointsize 80 -fill white -draw 'text {} "{}" ' {} {} """.format(imgtxtloc[c], imgtxt[c], imgf[-1], imgf[-1]))

		print 'Processing image ' + numbers[i]
		
		imgout = 'tmp/imgout{}.png'.format(numbers[i])
		os.system(""" convert -pointsize 80 -fill white +append -draw 'text 50,80 "Frame: {}" ' {} {} """.format(i, ' '.join(imgf), imgout))
		
	os.system('ffmpeg -threads 2 -start_number {} -framerate {} -i tmp/imgout%0{}d.png {}'.format(startfrom, frameRate, len(digits), output))
	os.system('rm -r tmp')

make_movie(dir='movie2', imgs = ['compress'])

#make_movie(dir='movie2', imgs = ['nofriction', 'friction', 'stages'], \
#	imgtxt = ['no friction', 'friction', 'friction (3 stages)'], imgtxtloc = ['50,1970', '50, 1970', '50, 1970'])

# classical example 
# mixed up friction vs nofriction images with ovito ~ oops!
# make_movie(dir='movie', imgs = ['friction', 'nofriction', 'stages'], \
#	imgtxt = ['no friction', 'friction', 'friction (3 stages)'], imgtxtloc = ['150,990', '150, 990', '150, 990'])

# e.g. on using 'startfrom' to specify where image sequence starts from
# make_movie(dir='movie', imgs = ['beverloo-bottom', 'beverloo-side', 'beverloo-core'], \
#	imgtxt = ['bottom view', 'side view', 'core view'], imgtxtloc = ['150,990', '150,990', '150,990'], startfrom = 80)