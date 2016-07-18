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

class Visualizer:
    """ A general purpose class for visualizing data using vtk """
    
    def __init__(self):
        #self._ren.AddActor(self._axes)
        
        self._axes = vtk.vtkAxesActor()
        self._ren = vtk.vtkRenderer() 
        
        self._style = vtk.vtkInteractorStyleTrackballCamera()

        self._renWin = vtk.vtkRenderWindow()
        self._renWin.AddRenderer(self._ren)
        self._renWin.SetWindowName("DEM Visualizer")

        self._iren = vtk.vtkRenderWindowInteractor()
        self._iren.SetRenderWindow(self._renWin)

    def render(self):
        """ initialize renderer """
        
        self._ren.SetBackground(0, 0, 0)
        #ren.AddActor(arrowActor)

        self._iren.Initialize()
        self._iren.SetInteractorStyle(self._style)
        self._iren.Start()

        self._close_window()

    def loadStl(self, fname):
        """Load the given STL file into a vtkPolyData object"""

        reader = vtk.vtkSTLReader()
        reader.SetFileName(fname)
        reader.Update()
        self._stl = reader.GetOutput() # polydata

    @property
    def axes(self):
        return self.axes

    def visSpheres(self, positions, velocities, radius):

        self._camera = vtk.vtkCamera()
        self._camera.SetFocalPoint(0, 0, 0)
        self._camera.SetPosition(np.mean(positions, axis=0))
        self._ren.SetActiveCamera(self._camera)

        self.positions_vtk = numpy_support.numpy_to_vtk(num_array=positions.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        self.positions_vtk.SetName("positions")

        self.radius_vtk = numpy_support.numpy_to_vtk(num_array=radius.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        self.radius_vtk.SetName("radius")

        self.points = vtk.vtkPoints()
        self.velPoints = vtk.vtkPoints()
        
        for i, r in enumerate(positions):
            self.points.InsertPoint(i, r[0], r[1], r[2])

        for i, v in enumerate(velocities):
            self.velPoints.InsertPoint(i, v[0], v[1], v[2])

        self.source = vtk.vtkSphereSource()
        self.source.SetRadius(1.0)

        self.ballGlyph = vtk.vtkGlyph3D()

        if vtk.VTK_MAJOR_VERSION <= 5:
            self.ballGlyph.SetSource(self.source.GetOutput())
        else:
            self.ballGlyph.SetSourceConnection(self.source.GetOutputPort())

        self.grid = vtk.vtkUnstructuredGrid()
        self.grid.SetPoints(self.points)
        self.grid.GetPointData().AddArray(self.radius_vtk)
        self.grid.GetPointData().SetActiveScalars("radius") # this scales the source (sphere) radius (1.0)

        self.ballGlyph.SetInputData(self.grid)

        # Create a color transfer function to be used for both the balls and arrows.
        colorTransferFunction = vtk.vtkColorTransferFunction()
        colorTransferFunction.AddRGBPoint(5.0 , 0.0, 0.0, 1.0)
        colorTransferFunction.AddRGBPoint(10.0, 0.0, 1.0, 1.0)
        colorTransferFunction.AddRGBPoint(15.0, 0.0, 1.0, 0.0)
        colorTransferFunction.AddRGBPoint(20.0, 1.0, 1.0, 0.0)
        colorTransferFunction.AddRGBPoint(25.0, 1.0, 0.0, 0.0)
        colorTransferFunction.AddRGBPoint(30.0, 1.0, 0.0, 1.0)

        #ballGlyph.SetScaleModeToDataScalingOn() 
        self.mapper = vtk.vtkPolyDataMapper()

        if vtk.VTK_MAJOR_VERSION <= 5:
           self.mapper.SetInput(self.ballGlyph.GetOutput())
        else:
           self.mapper.SetInputConnection(self.ballGlyph.GetOutputPort())

        # Set colors depending on the color transfer functions
        self.mapper.SetLookupTable(colorTransferFunction)

        #Put an arrow (vector) at each ball
        arrow = vtk.vtkArrowSource()
        arrow.SetTipRadius(0.05)
        arrow.SetShaftRadius(0.05)

        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(self.points)

        arrowGlyph = vtk.vtkGlyph3D()
        arrowGlyph.SetInputData(self.poly)
        arrowGlyph.SetSourceConnection(arrow.GetOutputPort())
        arrowGlyph.SetScaleFactor(1e-3)

        # We do not want the Arrow's size to depend on the Scalar
        #arrowGlyph.SetScaleModeToDataScalingOff()
        arrowMapper = vtk.vtkPolyDataMapper()
        arrowMapper.SetInputConnection(arrowGlyph.GetOutputPort())

        # Set colors depending on the color transfer functions
        arrowMapper.SetLookupTable(colorTransferFunction)

        # actor
        self.ballActor = vtk.vtkActor()
        self.ballActor.SetMapper(self.mapper)
     
        arrowActor = vtk.vtkActor()
        arrowActor.SetMapper(arrowMapper)
        arrowActor.GetProperty().SetColor(0,1,1)

        self._ren.AddActor(self.ballActor)
        self._ren.AddActor(arrowActor)

    def _close_window(self):
        """ kills any active renderer / windows """
        render_window = self._iren.GetRenderWindow()
        render_window.Finalize()
        self._iren.TerminateApp()

        del render_window, self._iren, self._ren, self._renWin