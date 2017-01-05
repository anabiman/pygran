from PIL import Image
from numpy import random, array, linspace, abs

LENGTH, WIDTH = 1024, 1024
img = Image.new( 'RGB', (LENGTH, WIDTH), "black") # create a new black image
pixels = img.load() # create the pixel map

N = 10
x = linspace(0, 0.99, N) #random.rand(255)
x = array(x * img.size[0], 'int')

y = linspace(0, 0.99, N) #random.rand(255)
y = array(y * img.size[1], 'int')

r = random.randint(10, 50, N*N)
count = 0

for i in x:
    for j in y:
	radius = r[count]
	count += 1

	if (i + radius < LENGTH and i - radius > 0) and (j + radius < WIDTH and j - radius > 0):
		for ic in range(-radius, radius+1):
			for jc in range(-radius, radius+1):
				if (ic)**2 + (jc)**2 <= radius**2:
		        		pixels[i+ic,j+jc] = (255, 255, 255) # set the colour accordingly
	else:
		pass

img.show()
