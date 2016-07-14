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
from vtk.util import numpy_support

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

def plotSpheres(ren, points, radius):

    source = vtk.vtkSphereSource()
    source.SetRadius(0.01)

    ballGlyph = vtk.vtkGlyph3D()

    if vtk.VTK_MAJOR_VERSION <= 5:
        ballGlyph.SetSource(source.GetOutput())
    else:
        ballGlyph.SetSourceConnection(source.GetOutputPort())

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.GetPointData().AddArray(radius)
    grid.GetPointData().SetActiveScalars("radius")

    ballGlyph.SetInputData(grid)

    #ballGlyph.SetScaleModeToDataScalingOn() 
    mapper = vtk.vtkPolyDataMapper()

    if vtk.VTK_MAJOR_VERSION <= 5:
	   mapper.SetInput(ballGlyph.GetOutput())
    else:
	   mapper.SetInputConnection(ballGlyph.GetOutputPort())
 
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
 
    # assign actor to the renderer
    ren.AddActor(actor)

def visSpheres(positions, radius):

    positions_vtk = numpy_support.numpy_to_vtk(num_array=positions.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    positions_vtk.SetName("positions")

    radius_vtk = numpy_support.numpy_to_vtk(num_array=radius.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    radius_vtk.SetName("radius")

    points = vtk.vtkPoints()
    
    for i, r in enumerate(positions):
        points.InsertPoint(i, r[0], r[1], r[2])

    xc, yc, zc = np.mean(positions, axis=0)
    
    axes = vtk.vtkAxesActor()
    axes.SetPosition(xc, yc, zc)
    axes.VisibilityOff()

    ren = vtk.vtkRenderer() 
    ren.AddActor(axes)

     # Create a RenderWindowInteractor to permit manipulating the camera
    plotSpheres(ren, points, radius_vtk)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)

    camera = vtk.vtkCamera()
    camera.SetFocalPoint(0, 0, 0)
    camera.SetPosition(xc, yc, zc)
    ren.SetActiveCamera(camera)
    
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
