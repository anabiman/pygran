'''
Created on July 10, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
#
# ----------------------------------------------------------------------
#
#   Python module for visualing DEM simulations with VTK library
#
# ----------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

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

def plotSphereMesh(meshFname = None, dumpFname = None):

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
        with open(dumpFname) as fp:
            lines = (line for line in fp if not line.strip()[:3].isalpha())

            for line in lines:
                data = line.split()
                
                if len(data) == 1:
                    Natoms = int(line)
                    sphereData = np.zeros(Natoms)
                else:
                    sphereData[count] = float(line)

            print sphereData.shape()

            for data in sphereData:
                plotSpheres(ren, data[1] * scale, data[2] * scale, data[3] * scale, 0.00224)

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

def visualize_Vtk(file_name):
    # Read the source file.
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()
    scalar_range = output.GetScalarRange()

    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(scalar_range)

    # Create the Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create the Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1) # Set background to white

    # Create the RendererWindow
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    # Create the RendererWindowInteractor and display the vtk_file
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)
    interactor.Initialize()
    interactor.Start()
