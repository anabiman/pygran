from PyDEM import Visualizer
from numpy import array
from numpy.random import rand

N, scale = 1000, 10.0

pos = array([rand(N), rand(N), rand(N)]).T * scale

vel = 1.0 - 2.0 * array([rand(N), rand(N), rand(N)]).T
rad = rand(N) * 0.5

v = Visualizer.Visualizer()
v.attach_stl('hopper-2cm-6cm.stl', scale=(1e-1, 1e-1, 1e-1))

v.attach_pos(pos, rad)
v.attach_vel(vel, rad)
v.render()
