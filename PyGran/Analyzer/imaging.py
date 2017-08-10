from numbers import Number
from PIL import Image
from numpy import random, array, linspace, sqrt, fabs
import numpy as np

try:
	import cv2
except:
	pass

def readExcel(fname):
	"""
	reads an excel sheet(s) and appends each data column to 
	a dictionary
	"""
	from xlrd import open_workbook

	book = open_workbook(fname, on_demand=True)
	data = {}

	for name in book.sheet_names():
		sheet = book.sheet_by_name(name)

		for coln in range(sheet.ncols):
			for celln, cell in enumerate(sheet.col(coln)):
				if celln == 0:
					dname = cell.value
					data[dname] = []
				else:
					if isinstance(cell.value, Number):
						data[dname].append(cell.value)
					else:
						print cell.value, ' ignored'

	for key in data.keys():
		data[key] = array(data[key])

	return data

def _mapPositions(Particles, axis, resol=None):
	""" Maps particle positions to pixels """
	if axis == 'x':
		h = Particles.x

		if resol:
			x = array(Particles.y * resol, 'int')
			y = array(Particles.z * resol, 'int')
		else:
			x = Particles.y
			y = Particles.z

	elif axis == 'y':
		h = Particles.y

		if resol:
			x = array(Particles.x * resol, 'int')
			y = array(Particles.z * resol, 'int')
		else:
			x = Particles.x
			y = Particles.z

	elif axis == 'z':
		h = Particles.z

		if resol:
			x = array(Particles.x * resol, 'int')
			y = array(Particles.y * resol, 'int')
		else:
			x = Particles.x
			y = Particles.y
	else:
		raise ValueError('axis can be x, y, or z only.')

	return x,y,h

def slice(Particles, zmin, zmax, axis, resol=None, output=None, imgShow=False):
	"""
	Generates a 2D image from a slice (limited by 'zmin/zmax' and of reslution '1/dz') 
	of a 3D config in the Particles class. The resol is in distance per pixel.
	"""

	Particles = Particles.copy()
	maxRad = Particles.radius.max()

	if resol:
		resol = 1.0 / resol

	x,y,h = _mapPositions(Particles, axis, resol)

	if resol:
		length, width = max(x), max(y)
		img = Image.new('RGB', (length, width), "black") # create a new black image
		pixels = img.load() # create the pixel map

	Particles = Particles[h >= zmin - maxRad]

	x,y,h = _mapPositions(Particles, axis, resol)
	Particles = Particles[h <= zmax + maxRad]

	h = eval('Particles.{}'.format(axis))
	Particles = Particles[fabs(h - zmax) <= Particles.radius]

	h = eval('Particles.{}'.format(axis))
	Particles = Particles[fabs(h - zmin) <= Particles.radius]

	x,y,h = _mapPositions(Particles, axis, resol)
	N = Particles.natoms
	print 'natoms per slice = ', N
	if N > 0:

		# uncommnet the line below for alpha channel estimation
		# trans = array(z * 255 / z.max(), 'int')

		zmean = (zmax + zmin) / 2

		r = sqrt(Particles.radius**2.0 - (h - zmean)**2.0)

		if not resol:
			return Particles

		else:
			r = array(r * resol, 'int')

			for n in range(N):

				i, j = x[n], y[n]
				radius = r[n]
				
				for ic in range(-radius, radius+1):
					for jc in range(-radius, radius+1):
						if (ic)**2 + (jc)**2 <= radius**2:
							if ( (i + ic < length) and (i + ic >= 0) and (j + jc < width) and (j + jc >= 0) ):
								pixels[i+ic,j+jc] = (255, 255, 255) #  add trans[n] for transparency (e.g. png files) and then set the colour accordingly
							else:
								pass
			if imgShow:
				img.show()

			if output:
				img.save(output)

def reconstruct(fimg, imgShow=False):

	img = cv2.imread(fimg)
	gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
	ret, thresh = cv2.threshold(gray,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)

	# noise removal
	kernel = np.ones((3,3),np.uint8)
	opening = cv2.morphologyEx(thresh,cv2.MORPH_OPEN,kernel, iterations = 2)

	# sure background area
	sure_bg = cv2.dilate(opening,kernel,iterations=3)

	# Finding sure foreground area
	if int(cv2.__version__.split('.')[0]) < 3:
    		dist_transform = cv2.distanceTransform(opening, cv2.cv.CV_DIST_L2, 5)
	else:
		dist_transform = cv2.distanceTransform(opening, cv2.DIST_L2, 5)

   	ret, sure_fg = cv2.threshold(dist_transform, 0.7*dist_transform.max(), 255, 0)

   	# Finding unknown region
   	sure_fg = np.uint8(sure_fg)
   	unknown = cv2.subtract(sure_bg,sure_fg)

	if int(cv2.__version__.split('.')[0]) == 3:
		# Marker labelling
		ret, markers = cv2.connectedComponents(sure_fg)

    		# Add one to all labels to make sure background is not 0, but 1
    		markers = markers+1

    		# Now, mark the region of unknown with zero
		markers[unknown==255] = 0

		# Now our marker is ready. It is time for final step, apply watershed. 
		# Then marker image will be modified. The boundary region will be marked with -1.
		markers = cv2.watershed(img,markers)
		img[markers == -1] = [255,0,0]
	else:

		img = cv2.add(sure_bg,sure_fg)

	if imgShow:
		cv2.imshow('image', img)
		#cv2.imwrite('image.png',img)
		cv2.waitKey(0)
                cv2.destroyAllWindows()

	return img
