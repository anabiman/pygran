'''
Created on Sep 8, 2016
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
 '''

import vtk
import numpy as np
import sys
 
def loadStl(fname):
	"""Load the given STL file, and return a vtkPolyData object for it."""

	reader = vtk.vtkSTLReader()
	reader.SetFileName(fname)
	reader.Update()
	polydata = reader.GetOutput()
	return polydata

def polyDataToActor(polydata):
	"""Wrap the provided vtkPolyData object in a mapper and an actor, returning the actor."""

	mapper = vtk.vtkPolyDataMapper()
	
	if vtk.VTK_MAJOR_VERSION <= 5:
		mapper.SetInput(polydata)
	else:
		mapper.SetInputData(polydata)
	
	actor = vtk.vtkActor()
	actor.SetMapper(mapper)

	actor.GetProperty().SetColor(0.5, 0.5, 1.0)
	return actor

def plotSpheres(ren, *args):

	source = vtk.vtkSphereSource()
	source.SetCenter(args[0], args[1], args[2])
	source.SetRadius(args[3])

	# mapper
	mapper = vtk.vtkPolyDataMapper()

	if vtk.VTK_MAJOR_VERSION <= 5:
	   mapper.SetInput(source.GetOutput())
	else:
	   mapper.SetInputConnection(source.GetOutputPort())
 
	# actor
	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
 
	# assign actor to the renderer
	ren.AddActor(actor)

def initialize():

	# Create a rendering window and renderer
	transform = vtk.vtkTransform()
	transform.Scale(10.0, 10.0, 10.0)
	 
	axes = vtk.vtkAxesActor()
	axes.SetUserTransform(transform)

	transform.Translate(3.0, -2.0, 0.0)
	axes.SetUserTransform(transform)

	ren = vtk.vtkRenderer()
	ren.AddActor(axes)
	ren.ResetCamera()

	renWin = vtk.vtkRenderWindow()
	renWin.AddRenderer(ren)

	# Create a RenderWindowInteractor to permit manipulating the camera
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetRenderWindow(renWin)
	style = vtk.vtkInteractorStyleTrackballCamera()
	iren.SetInteractorStyle(style)

	# enable user interface interactor
	iren.Initialize()
	renWin.Render()
	iren.Start()

	return ren, iren, renWin

def visualize(meshFname = None, dumpFname = None, nFrames = 100):

	# Create a rendering window and renderer
	transform = vtk.vtkTransform()
	transform.Scale(10.0, 10.0, 10.0)
	 
	axes = vtk.vtkAxesActor()
	axes.SetUserTransform(transform)

	transform.Translate(2.0, -1.0, 0.0)
	axes.SetUserTransform(transform)

	renderer = vtk.vtkRenderer()
	renderer.AddActor(axes)
	
	camera = vtk.vtkCamera()
	camera.SetFocalPoint(0, 0, 0);
	renderer.SetActiveCamera(camera);

	renWin = vtk.vtkRenderWindow()
	renWin.AddRenderer(renderer)

	 # Create a RenderWindowInteractor to permit manipulating the camera
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetRenderWindow(renWin)
	style = vtk.vtkInteractorStyleTrackballCamera()
	iren.SetInteractorStyle(style)

	scale = 0.1

	if dumpFname is not None:
		with open(dumFname) as fp:
			lines = (line for line in fp if not line.strip()[:3].isalpha())
			sphereData = np.array(lines)

		Natoms = len(sphereData) - nFrames

		for i in range(nFrames):
			for data in sphereData[Natoms * nFrames: Natoms * (nFrames+1)]:
				plotSpheres(ren, data[1] * scale, data[2] * scale, data[3] * scale, 0.001)

	if meshFname is not None:
		polydata = loadStl(meshFname)
		renderer.AddActor(polyDataToActor(polydata))
		renderer.SetBackground(0.1, 0.1, 0.1)

	style = vtk.vtkInteractorStyleTrackballCamera()
	iren.SetInteractorStyle(style)
	renderer.ResetCamera()

	iren.Initialize()
	renWin.Render()
	iren.Start()
