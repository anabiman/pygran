'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

from PyDEM import Visualizer, Analyzer
from numpy import array
from numpy.random import rand
import wx

N, scale = 1000, 10.0

g = Analyzer.Granular('traj.dump')
g.next()

app = wx.App(redirect=False)
frame = Visualizer.Visualizer(None, g.Particles, "Lights, Cameras, Action")
frame.Show()
app.MainLoop()
