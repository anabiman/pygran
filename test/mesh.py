import numpy
from stl import mesh

# Create 3 faces of a cube
data = numpy.zeros(6, dtype=mesh.Mesh.dtype)
scale = 20.0
z = 30.0

# Top of the cube
data['vectors'][0] = numpy.array([[-0.5, 0.5, z/scale],
                                  [0.5, -0.5, z/scale],
                                  [-0.5, -0.5, z/scale]]) * scale
data['vectors'][1] = numpy.array([[0.5, -0.5, z/scale],
                                  [-0.5, 0.5, z/scale],
                                  [0.5, 0.5, z/scale]]) * scale

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
