
# !/usr/bin/python
# -*- coding: utf8 -*- 

'''
Created on July 10, 2016
@author: Andrew Abi-Mansour

  This is the 
   __________         ________                     
  ██████╗ ██╗   ██╗ ██████╗ ██████╗  █████╗ ███╗   ██╗
  ██╔══██╗╚██╗ ██╔╝██╔════╝ ██╔══██╗██╔══██╗████╗  ██║
  ██████╔╝ ╚████╔╝ ██║  ███╗██████╔╝███████║██╔██╗ ██║
  ██╔═══╝   ╚██╔╝  ██║   ██║██╔══██╗██╔══██║██║╚██╗██║
  ██║        ██║   ╚██████╔╝██║  ██║██║  ██║██║ ╚████║
  ╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝
													  
  DEM simulation and analysis toolkit
  http://www.pygran.org, support@pygran.org

  Core developer and main author:
  Andrew Abi-Mansour, andrew.abi.mansour@pygran.org

  PyGran is open-source, distributed under the terms of the GNU Public
  License, version 2 or later. It is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
  received a copy of the GNU General Public License along with PyGran.
  If not, see http://www.gnu.org/licenses . See also top-level README
  and LICENSE files.

-------------------------------------------------------------------------
	Python module for visualing DEM simulations with VTK library
-------------------------------------------------------------------------

'''

import vtk
from vtk.util import numpy_support
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

import numpy as np
import sys

from numpy.linalg import norm
import wx

class Panel(wx.Panel):
	""" A general purpose class for visualizing data using vtk """

	def __init__(self, parent):

		wx.Panel.__init__(self, parent)
		 
		#to interact with the scene using the mouse use an instance of vtkRenderWindowInteractor. 
		self.widget = wxVTKRenderWindowInteractor(self, -1)
		self.widget.Enable(1)
		self.widget.AddObserver("ExitEvent", lambda o,e,f=self: f.Close())
		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.sizer.Add(self.widget, 1, wx.EXPAND)
		self.SetSizer(self.sizer)
		self.Layout()
		self.isploted = False
		
		self._init()

		self._pos = None
		self._vel = None
		self._rad = None

	def _init(self):

		self._axes = vtk.vtkAxesActor()
		self._marker = vtk.vtkOrientationMarkerWidget()
		self._marker.SetInteractor( self.widget._Iren )
		self._marker.SetOrientationMarker( self._axes )
		self._marker.SetViewport(0.75,0,1,0.25)
		self._marker.SetEnabled(1)

		self._ren = vtk.vtkRenderer() 
		self._camera = self._ren.GetActiveCamera()

		self._camera.SetFocalPoint(0, 0, 0.05)
		self._camera.SetPosition(0,0,0)
		self._ren.SetActiveCamera(self._camera)

		self._style = vtk.vtkInteractorStyleTrackballCamera()

		self._renWin = self.widget.GetRenderWindow()

		self._renWin.AddRenderer(self._ren)
		self._renWin.SetWindowName("DEM Visualizer")
		
		self.widget._Iren.SetRenderWindow(self._renWin)

		#self._ren.AddActor(self._axes)

	def _setupColorFunction(self, minV, maxV):

		# Create a color transfer function to be used for both the balls and arrows.
		self._colorTransferFunction = vtk.vtkColorTransferFunction()
		self._colorTransferFunction.AddRGBPoint(minV, 1.0, 0.0, 0.0)
		self._colorTransferFunction.AddRGBPoint(0.5*(minV + maxV), 0.0, 1.0, 0.0)
		self._colorTransferFunction.AddRGBPoint(maxV, 0.0, 0.0, 1.0)

	def render(self):
		""" initialize renderer """
		
		self._ren.SetBackground(0, 0, 0)

		#self._close_window()

	def load_parts(self, Parts):

		self._pos = np.array([Parts.x, Parts.y, Parts.z]).T
		self._vel = np.array([Parts.vx, Parts.vy, Parts.vz]).T
		self._rad = Parts.radius
		self._points = vtk.vtkPoints()

		for i, r in enumerate(self._pos):
			self._points.InsertPoint(i, r[0], r[1], r[2])

	def attach_stl(self, fname, scale=None):
		"""Load a given STL file into a vtkPolyData object"""

		reader = vtk.vtkSTLReader()
		reader.SetFileName(fname)
		reader.Update() # polydata

		if scale is not None:
			trans = vtk.vtkTransform()
			trans.Scale(scale)

			filt = vtk.vtkTransformFilter()

			if vtk.VTK_MAJOR_VERSION <= 5:
				filt.SetInputConnection(reader.GetOutputPort)
			else:
				filt.SetInputConnection(reader.GetOutputPort())
				   
			filt.SetTransform(trans)

		mapper = vtk.vtkPolyDataMapper()

		if vtk.VTK_MAJOR_VERSION <= 5:
		   mapper.SetInput(filt.GetOutput())
		else:
		   mapper.SetInputConnection(filt.GetOutputPort())

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)

		self._ren.AddActor(actor)

	def loadVtk(self, fname):
		""" Load a given VTK file into a vtkPolyData object """

		reader = vtk.vtkUnstructuredGridReader()
		reader.SetFileName(filename)
		reader.Update()
		self._vtk = reader.GetOutput()

	@property
	def axes(self):
		return self.axes

	def attach_pos(self):

		pos = self._pos
		rad = self._rad

		if pos is not None and rad is not None:
			positions_vtk = numpy_support.numpy_to_vtk(num_array=pos.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
			positions_vtk.SetName("positions")

			radius_vtk = numpy_support.numpy_to_vtk(num_array=rad.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
			radius_vtk.SetName("radius")

			sphere = vtk.vtkSphereSource()
			sphere.SetRadius(1.0)

			ballGlyph = vtk.vtkGlyph3D()

			if vtk.VTK_MAJOR_VERSION <= 5:
				ballGlyph.SetSource(sphere.GetOutput())
			else:
				ballGlyph.SetSourceConnection(sphere.GetOutputPort())

			polydata = vtk.vtkPolyData()
			polydata.SetPoints(self._points)
			polydata.GetPointData().AddArray(radius_vtk)
			polydata.GetPointData().SetActiveScalars("radius") # this scales the source (sphere) radius (1.0)
			
			ballGlyph.SetInputData(polydata)

			#ballGlyph.SetScaleModeToDataScalingOn() 
			mapper = vtk.vtkPolyDataMapper()

			if vtk.VTK_MAJOR_VERSION <= 5:
			   mapper.SetInput(ballGlyph.GetOutput())
			else:
			   mapper.SetInputConnection(ballGlyph.GetOutputPort())

			# Set colors depending on the color transfer functions
			# mapper.SetLookupTable(self.colorTransferFunction)

			# actor
			ballActor = vtk.vtkActor()
			ballActor.GetProperty().SetColor(0,0,1)
			ballActor.SetMapper(mapper)

			#self._ren.AddActor(ballActor)
			self._marker2 = vtk.vtkOrientationMarkerWidget()
			self._marker2.SetInteractor( self.widget._Iren )
			self._marker2.SetOrientationMarker( ballActor )
			self._marker2.SetViewport(0.75,0,1,0.25)
			self._marker2.SetEnabled(1)

		else:
			print("No particles found. Make sure the particles loaded have positions and radii.")

	def attach_vel(self):

		vel = self._vel * 1e-2
		rad = self._rad

		if vel is not None and rad is not None:

			velMag = norm(vel, axis=1)
			velMag_vtk = numpy_support.numpy_to_vtk(num_array=velMag.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
			velMag_vtk.SetName("veloMag")

			vecVel = vtk.vtkFloatArray()
			vecVel.SetNumberOfComponents(3)

			for i, v in enumerate(vel):
				vecVel.InsertTuple3(i, v[0], v[1], v[2])

			#Put an arrow (vector) at each ball
			arrow = vtk.vtkArrowSource()
			arrow.SetTipRadius(rad.mean() * 10)
			arrow.SetShaftRadius(rad.mean() * 10)

			poly = vtk.vtkPolyData()
			poly.SetPoints(self._points)

			poly.GetPointData().AddArray(velMag_vtk)
			poly.GetPointData().SetActiveScalars("veloMag")

			arrowGlyph = vtk.vtkGlyph3D()
			
			arrowGlyph.SetInputData(poly)
			arrowGlyph.SetSourceConnection(arrow.GetOutputPort())
			arrowGlyph.SetVectorModeToUseVector()

			poly.GetPointData().SetVectors(vecVel)

			# If we do not want the Arrow's size to depend on the Scalar
			# then arrowGlyph.SetScaleModeToDataScalingOff() must be called
			arrowMapper = vtk.vtkPolyDataMapper()
			arrowMapper.SetInputConnection(arrowGlyph.GetOutputPort())

			self._addScalarBar(velMag)
			arrowMapper.SetLookupTable(self._colorTransferFunction)

			arrowActor = vtk.vtkActor()
			arrowActor.SetMapper(arrowMapper)
			arrowActor.GetProperty().SetColor(1,1,0)

			self._ren.AddActor(arrowActor)
		else:
			print("No particles found. Make sure the particles loaded have velocities and radii.")

	def _addScalarBar(self, val):

		if len(val.shape) > 1:
			val = norm(val, axis=1)

		self._setupColorFunction(val.min(), val.max())
		self._scalarBar = vtk.vtkScalarBarActor()
		self._scalarBar.SetLookupTable(self._colorTransferFunction)
		self._scalarBar.SetTitle("Velocities (m/s)")
		self._ren.AddActor2D(self._scalarBar)

	def _close_window(self):
		""" kills any active renderer / windows """
		render_window = self._iren.GetRenderWindow()
		render_window.Finalize()
		self._iren.TerminateApp()

		del render_window, self._iren, self._ren, self._renWin

class Visualizer(wx.Frame):
	def __init__(self,parent, particles, title):

		self._Particles = particles
		wx.Frame.__init__(self,parent,title=title,size=(800,600), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU| wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)
		self.sp = wx.SplitterWindow(self)
		self.p1 = Panel(self.sp)
		self.p2 = wx.Panel(self.sp,style=wx.SUNKEN_BORDER)
		 
		self.sp.SplitHorizontally(self.p1,self.p2,470)
 
		self.statusbar = self.CreateStatusBar()
		self.statusbar.SetStatusText("Click on the Plot Button")
		 
		self.plotbut = wx.Button(self.p2,-1,"plot", size=(80,40),pos=(10,10))
		self.plotbut.Bind(wx.EVT_BUTTON, self.plot)

		self.clearbut = wx.Button(self.p2,-1,"clear", size=(80,40),pos=(140,10))
		self.clearbut.Bind(wx.EVT_BUTTON, self.plot)

		self.loadbut = wx.Button(self.p2,-1,"load", size=(80,40),pos=(270,10))
		self.loadbut.Bind(wx.EVT_BUTTON, self.plot)
		 
		self.slidebut = wx.Slider(self.p2, value=200, minValue=150, maxValue=500, pos=(400, 10), size=(250, -1), style=wx.SL_HORIZONTAL)

		self.slidetxt = wx.StaticText(self.p2, label='resolution', pos=(475, 20)) 
 
	def plot(self,event):

		if not self.p1.isploted:

			self.p1.load_parts(self._Particles)
			#self.p1.attach_pos()
			#self.p1.attach_vel()
			self.p1.attach_stl('hopper-20-6.stl', scale=(1, 1, 1))

			self.statusbar.SetStatusText("Use your mouse to interact with the model")
			self.p1.render()
