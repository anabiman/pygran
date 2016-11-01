import numpy
from stl import mesh, Mode
from stl.base import RemoveDuplicates

tris2D1 = numpy.loadtxt('trisCoords1.dat', delimiter=',')
tris2D2 = numpy.loadtxt('trisCoords2.dat', delimiter=',')
tris2D3 = numpy.loadtxt('trisCoords3.dat', delimiter=',')


data = numpy.zeros(6 * len(tris2D1), dtype=mesh.Mesh.dtype)

for i in range(len(tris2D1)):

	x1, y1, z1 = tris2D1[i]
	x2, y2, z2 = tris2D2[i]
	x3, y3, z3 = tris2D3[i]

	data['vectors'][i] = numpy.array([[x1, y1, z1],
					 [x2, y2, z2],
					 [x3, y3, z3]])

meshes = mesh.Mesh(data, remove_duplicate_polygons=True)

meshes.save('hopper-6-20.stl', mode=Mode.ASCII)

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# Load the STL files and add the vectors to the plot
your_mesh = mesh.Mesh.from_file('hopper-6-20.stl')

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(meshes.vectors))

# Auto scale to the mesh size
scale = your_mesh.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()
