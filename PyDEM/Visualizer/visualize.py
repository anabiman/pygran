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
        self._init()
        
    def _init(self):

        self._axes = vtk.vtkAxesActor()
        self._ren = vtk.vtkRenderer() 
        
        self._style = vtk.vtkInteractorStyleTrackballCamera()

        self._renWin = vtk.vtkRenderWindow()
        self._renWin.AddRenderer(self._ren)
        self._renWin.SetWindowName("DEM Visualizer")

        self._iren = vtk.vtkRenderWindowInteractor()
        self._iren.SetRenderWindow(self._renWin)

        # Create a color transfer function to be used for both the balls and arrows.
        self.colorTransferFunction = vtk.vtkColorTransferFunction()
        self.colorTransferFunction.AddRGBPoint(-1.0 , 0.0, 0.0, 1.0)
        self.colorTransferFunction.AddRGBPoint(-0.5, 0.0, 1.0, 1.0)
        self.colorTransferFunction.AddRGBPoint(0, 0.0, 1.0, 0.0)
        self.colorTransferFunction.AddRGBPoint(0.5, 1.0, 1.0, 0.0)
        self.colorTransferFunction.AddRGBPoint(1, 1.0, 0.0, 0.0)

    def render(self):
        """ initialize renderer """
        
        self._ren.SetBackground(0, 0, 0)
        #ren.AddActor(arrowActor)

        self._iren.Initialize()
        self._iren.SetInteractorStyle(self._style)
        self._iren.Start()

        self._close_window()

        self._init()

    def loadStl(self, fname, scale=None):
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

        vecVel = vtk.vtkFloatArray()
        vecVel.SetNumberOfComponents(3)

        for i, v in enumerate(velocities):
            vecVel.InsertTuple3(i, v[0], v[1], v[2])
        
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

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(self.points)
        polydata.GetPointData().AddArray(self.radius_vtk)
        polydata.GetPointData().SetActiveScalars("radius") # this scales the source (sphere) radius (1.0)
        polydata.GetPointData().SetVectors(vecVel)

        self.ballGlyph.SetInputData(polydata)

        #ballGlyph.SetScaleModeToDataScalingOn() 
        mapper = vtk.vtkPolyDataMapper()

        if vtk.VTK_MAJOR_VERSION <= 5:
           mapper.SetInput(self.ballGlyph.GetOutput())
        else:
           mapper.SetInputConnection(self.ballGlyph.GetOutputPort())

        # Set colors depending on the color transfer functions
        # mapper.SetLookupTable(self.colorTransferFunction)

        #Put an arrow (vector) at each ball
        arrow = vtk.vtkArrowSource()
        arrow.SetTipRadius(0.05)
        arrow.SetShaftRadius(0.01)

        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(self.velPoints)

        arrowGlyph = vtk.vtkGlyph3D()
        
        arrowGlyph.SetInputData(self.poly)
        arrowGlyph.SetSourceConnection(arrow.GetOutputPort())
        arrowGlyph.SetScaleFactor(1e-1)

        arrowGlyph.SetScaleModeToDataScalingOff()
        arrowGlyph.SetVectorModeToUseVector()

        # We do not want the Arrow's size to depend on the Scalar
        #arrowGlyph.SetScaleModeToDataScalingOff()
        arrowMapper = vtk.vtkPolyDataMapper()
        arrowMapper.SetInputConnection(arrowGlyph.GetOutputPort())

        # actor
        ballActor = vtk.vtkActor()
        ballActor.GetProperty().SetColor(0,0,1)
        ballActor.SetMapper(mapper)
     
        arrowActor = vtk.vtkActor()
        arrowActor.SetMapper(arrowMapper)
        arrowActor.GetProperty().SetColor(1,1,0)

        self._ren.AddActor(ballActor)
        self._ren.AddActor(arrowActor)

    def addScalarBar(self):

        self.scalarBar = vtk.vtkScalarBarActor()
        self.scalarBar.SetLookupTable(self.colorTransferFunction)
        self.scalarBar.SetTitle("Colorbar")
        self._ren.AddActor2D(self.scalarBar)

    def _close_window(self):
        """ kills any active renderer / windows """
        render_window = self._iren.GetRenderWindow()
        render_window.Finalize()
        self._iren.TerminateApp()

        del render_window, self._iren, self._ren, self._renWin