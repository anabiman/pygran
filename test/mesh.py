import numpy
from stl import mesh

# Create 3 faces of a cube
scale = 30.0
z = 60.0

tris2D1 = numpy.loadtxt('postris1.dat', delimiter=',') * scale
tris2D2 = numpy.loadtxt('postris2.dat', delimiter=',') * scale
tris2D3 = numpy.loadtxt('postris3.dat', delimiter=',') * scale

tris2D1[:,1] *= 2.0
tris2D2[:,1] *= 2.0
tris2D3[:,1] *= 2.0

cog = (tris2D1.sum(axis=0) + tris2D2.sum(axis=0) + tris2D3.sum(axis=0)) / len(tris2D1)
print cog

tris2D1 -= cog
tris2D2 -= cog
tris2D3 -= cog


data = numpy.zeros(6 * len(tris2D1), dtype=mesh.Mesh.dtype)

for i in range(len(tris2D1)):

	x1, y1 = tris2D1[i]
	x2, y2 = tris2D2[i]
	x3, y3 = tris2D3[i]

	data['vectors'][i] = numpy.array([[x1, y1, z],
					 [x2, y2, z],
					 [x3, y3, z]])
# Top of the cube
#data['vectors'][0] = numpy.array([[-0.5, 0.5, z/scale],
#                                  [0.5, -0.5, z/scale],
#                                  [-0.5, -0.5, z/scale]]) * scale

#data['vectors'][1] = numpy.array([[0.5, -0.5, z/scale],
#                                  [-0.5, 0.5, z/scale],
#                                  [0.5, 0.5, z/scale]]) * scale

meshes = mesh.Mesh(data)

meshes.save('example.stl')

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# Load the STL files and add the vectors to the plot
your_mesh = mesh.Mesh.from_file('hopper-2cm-6cm.stl')

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(meshes.vectors))

# Auto scale to the mesh size
scale = your_mesh.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()
