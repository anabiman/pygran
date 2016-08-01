from PyDEM import Visualizer
from numpy import array
from numpy.random import rand

N, scale = 1000, 1.0

v = Visualizer.Visualizer()
v.loadStl('hopper-2cm-6cm.stl', scale=(1e-3, 1e-3, 1e-3))

pos = array([rand(N), rand(N), rand(N)]).T * scale
pos -= pos.mean(axis=0)
pos[2] += 1

vel = 1.0 - 2.0 * array([rand(N), rand(N), rand(N)]).T
print vel
rad = rand(N) * 1e-2

v.visSpheres(pos, vel, rad)
v.addScalarBar()
v.render()
