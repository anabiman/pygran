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
 
def surfMesh(positions):
    """ Creates a Delaunay traingulation from a set of particle positions.
    """

    points = vtk.vtkPoints()

    for i, r in enumerate(positions):
        points.InsertPoint(i, r[0], r[1], r[2])
        profile = vtk.vtkPolyData()
        profile.SetPoints(points)

    delny = vtk.vtkDelaunay3D()
    delny.SetInputData(profile)
    delny.SetTolerance(0.01)
    delny.SetAlpha(0.2)
    delny.BoundingTriangulationOff()

    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(delny.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1, 0, 0)

    # Create the Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1) # Set background to white

    # Create the RendererWindow
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)
    renderer_window.SetSize(250, 250)

    # Create the RendererWindowInteractor and display the vtk_file
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    interactor.Initialize()
    renderer_window.Render()
    interactor.Start()

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

def plotSpheres(ren, positions, radii):

    source = vtk.vtkSphereSource()
    source.SetCenter(positions[:,0], positions[:,1], positions[:,2])
    source.SetRadius(radii)

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

def visSpheres(positions, radii):

    axes = vtk.vtkAxesActor()

    ren = vtk.vtkRenderer()
    ren.AddActor(axes)
    ren.ResetCamera()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)

     # Create a RenderWindowInteractor to permit manipulating the camera
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)

    plotSpheres(ren, positions, radii)

    camera = vtk.vtkCamera()
    camera.SetFocalPoint(0, 0, 0);
    renderer.SetActiveCamera(camera);

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
