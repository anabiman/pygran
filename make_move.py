import os, glob

def make_movie(dir, img, output='out.mp4', frameRate=10):

	os.system('mkdir tmp')
	os.system('cp {}/*.png tmp/'.format(dir))

	files = glob.glob("{}/*.png".format(dir))
	files = sorted(files, key=lambda item: int( (item.split(dir+'/'+img)[1]).split('.png')[0] )) 

	numbers = []
	digits = str(len(files))

	for i in range(len(files)):
		num = str(i)
		numbers.append('0' * ( len(digits) - len(num) ) + num)

	count = 0
	for file in files:

		imgf = 'img' + numbers[count] + '.png'
		print file
		os.system(""" convert -pointsize 40 -fill white -draw 'text 5,100 "Frame: {}" ' {} tmp/{} """.format(count, file, imgf))
		count += 1

	os.system('ffmpeg -framerate {} -i tmp/img%0{}d.png -c:v libx264 -pix_fmt yuv420p {}'.format(frameRate, len(digits), output))
	os.system('rm -r tmp')


make_movie(dir='Cohesion-movie', img='cohes')
