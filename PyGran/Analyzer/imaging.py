'''
Created on February 24, 2017
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------


from numbers import Number
from PIL import Image
import glob, os
from numpy import random, array, linspace, sqrt, fabs
import numpy as np

try:
	import cv2
except:
	pass

def readExcel(fname):
	"""
	Reads an excel sheet(s) and appends each data column to 
	a dictionary

	@fname: filename to read
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

def slice(Particles, zmin, zmax, axis, resol=None, output=None, size=None, imgShow=False):
	"""
	Generates a 2D image from a slice (limited by 'zmin/zmax' and of reslution '1/dz') 
	of a 3D config in the Particles class. The resol is in distance per pixel.
	@[size]: tuple (length, width) for generated image slice
	"""

	Particles = Particles.copy()
	maxRad = Particles.radius.max()

	if resol:
		resol = 1.0 / resol

	x,y,h = _mapPositions(Particles, axis, resol)

	if resol:
		if size:
			length, width = size
		else:
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

def readImg(file, order=False, processes=None, fillHoles=False):
	""" Loads image file(s) and returns an array 

	@file: a list of image file names, or a string containing the image filename(s). In
	the latter case, if the string ends in '*' (e.g. img*), then all image files starting
	with 'img' are read (in a chronological order if order is set to True).

	@[order]: read a list of image files chronologically
	@[processes]: run in parallel???
	@[fillHoles]: fill holes in 3D image

	TODO: fix processes
	"""

	if type(file) is list:
		pass
	elif type(file) is str:
		if file.endswith('*'):
			file = glob.glob("{}*".format(file))
			file.sort(key=os.path.getmtime)
		else:
			 # Read single image file
			 pic = Image.open(file)

			 if len(np.array(pic.getdata()).shape) > 1:
			 	data = np.array(pic.getdata()).reshape(pic.size[0], pic.size[1], np.array(pic.getdata()).shape[-1])
			 else:
			 	data = np.array(pic.getdata()).reshape(pic.size[0], pic.size[1])

			 return data

		def func(file):
			
			for i, img in enumerate(file):
				pic = Image.open(img)

				if i == 0:
					if len(np.array(pic.getdata()).shape) > 1:
				 		data = np.zeros((pic.size[0], pic.size[1], np.array(pic.getdata()).shape[-1], len(file)))
				 	else:
				 		data = np.zeros((pic.size[0], pic.size[1], len(file)))

				if len(np.array(pic.getdata()).shape) > 1:
					data[:,:,:,i] = np.array(pic.getdata()).reshape(pic.size[0], pic.size[1], np.array(pic.getdata()).shape[-1])
				else:
					data[:,:,i] = np.array(pic.getdata()).reshape(pic.size[0], pic.size[1])

			if fillHoles:
				from scipy import ndimage
				return ndimage.morphology.binary_fill_holes(data).astype(int)
			else:
				return data

		if processes:
			pool = multiprocessing.Pool(processes=processes)
			func = partial(file)
			pool.map(func, range(nFiles))
			pool.close()
			pool.join()
		else:
			return func(file)

def coarseDiscretize(images, binsize, order=False, fillHoles=False):
	""" Discretizes a 3D image into a coarse grid
	@images: list of image file strings
	@binsize: length of each discrete grid cell in pixels
	@[order]: read images in a chronological order if set to True
	@[fillHoles]: fill holes in 3D image

	Returns a list of volume fraction 3D arrays, vol fraction mean, 
	and vol fraction variance
	"""

	dataList = []

	# Construct a 3D representation of the system
	if isinstance(images, list):
		for imgs in images:
			im = readImg(imgs, order)
			if len(im.shape) > 3:
				dataList.append(im[:,:,0,:])
			else:
				dataList.append(im)
	else:
		im = readImg(images, order)
		dataList.append(im)

	# Discretize system into cells of size 'binsize'
	data = dataList[0]

	xmin, xmax = 0, data.shape[0]
	ymin, ymax = 0, data.shape[1]
	zmin, zmax = 0, data.shape[2]

	x = np.array(xrange(xmin, xmax))
	y = np.array(xrange(ymin, ymax))
	z = np.array(xrange(zmin, zmax))

	indi = array(x / binsize, dtype='int32')
	indj = array(y / binsize, dtype='int32')
	indk = array(z / binsize, dtype='int32')

	# Compute variance in volume fraction for each grid cell
	volFrac = []

	for data in dataList:

		Grid = np.zeros((max(indi) + 1, max(indj) + 1,max(indk) + 1))

		for i, ig in enumerate(indi):
			for j, jg in enumerate(indj):
				for k, kg in enumerate(indk):
					Grid[ig,jg,kg] += data[i,j,k]

		volFrac.append(Grid)

	fracTotal = volFrac[0] * 0
	dataVar = []
	dataMean = []

	if len(dataList) > 1:
		for frac in volFrac:
			fracTotal += frac

		# Ignore all voxels not contaning any data
		for frac in volFrac:
			frac[fracTotal > 0] /= fracTotal[fracTotal > 0]
			dataMean.append(frac[fracTotal > 0].mean())
			dataVar.append(frac[fracTotal > 0].std()**2.0)
	else:
		data = dataList[0]
		c1, c2 = data[data > 0].min(), data.max()
		volFrac = [data * 0, data * 0]
		volFrac[0][data == c1], volFrac[1][data == c2] = data[data == c1] / c1, data[data == c2] / c2

		v1, v2 = volFrac[0], volFrac[1]
		dataMean = [v1[data > 0].mean(), v2[data > 0].mean()]
		dataVar = [v1[data > 0].std()**2, v2[data > 0].std()**2]

	return volFrac, dataMean, dataVar

def intensitySegregation(images, binsize, order=False):
	""" Computes the intensity of segregation from a set of image files
	@images: list of image file strings
	@binsize: length of each discrete grid cell in pixels
	@[order]: read images in a chronological order if set to True

	Returns the mean volume fraction, variance, and intensity

	TODO: support multi-component systems, not just binary systems
	"""

	_, dataMean, dataVar = coarseDiscretize(images, binsize, order)

	# Assuming only a binary system .. must be somehow fixed/extended for tertiary systems, etc.
	return dataMean, dataVar[0], dataVar[0] / (dataMean[0] * dataMean[1])

def scaleSegregation(images, binsize, samplesize, resol, maxDist=None, order=False, fillHoles=False):
	""" Computes (through Monte Carlo sim) the linear scale of segregation from a set of image files
	@images: list of image file strings
	@binsize: length of each discrete grid cell in pixels
	@samplesize: number of successful Monte Carlo trials
	@resol: image resolution (distance/pixel)

	@[maxDist]: maximum distance (in pixels) to sample
	@[order]: read images in a chronological order if set to True
	@[fillHoles]: fill holes in 3D image

	Returns the coefficient of correlation R(r) and separation distance (r)

	TODO: support multi-component systems, not just binary systems
	"""

	volFrac, volMean, volVar = coarseDiscretize(images, binsize, order, fillHoles)

	a,b = volFrac[0], volFrac[1]

	volMean = volMean[0]
	volVar = volVar[0]

	# Begin Monte Carlo simulation
	maxDim = max(a.shape)
	
	if not maxDist:
		maxDist = int(np.sqrt(3 * maxDim**2)) + 1

	corrfunc = np.zeros(maxDist - 1)
	count = 0
	incr = 1

	for dist in np.arange(incr, maxDist, incr):
		nTrials = 0
		corr = 0

		while nTrials < samplesize:

			theta, phi = np.random.rand() * np.pi, np.random.rand() * 2 * np.pi
			i1 = np.random.randint(0, a.shape[0], size=1)
			j1 = np.random.randint(0, a.shape[1], size=1)
			k1 = np.random.randint(0, a.shape[2], size=1)

			i2 = i1 + np.int(dist * np.sin(theta) * np.sin(phi))
			j2 = j1 + np.int(dist * np.sin(theta) * np.cos(phi))
			k2 = k1 + np.int(dist * np.cos(theta))

			# Check for boundary pts
			if i2 < a.shape[0] and j2 < a.shape[1] and k2 < a.shape[2]:

				# Make sure we are sampling non-void spatial points
				if a[i1,j1,k1] > 0 or b[i1,j1,k1] > 0:
					if a[i2,j2,k2] > 0 or b[i2,j2,k2] > 0:

						corr += ((a[i1,j1,k1] - volMean) * (a[i2,j2,k2] - volMean)) / volVar
						nTrials += 1

		corrfunc[count] = corr / samplesize
		count += 1

	return corrfunc, np.arange(incr, maxDist, incr) * resol * binsize
