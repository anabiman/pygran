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
from numpy.linalg import norm

class Visualizer:
    """ A general purpose class for visualizing data using vtk """
    
    def __init__(self):
        self._init()
        
    def _init(self):

        self._axes = vtk.vtkAxesActor()

        self._ren = vtk.vtkRenderer() 
        self._camera = vtk.vtkCamera()

        self._camera.SetFocalPoint(0, 0, 0)
        self._camera.SetPosition(0,0,0)
        self._ren.SetActiveCamera(self._camera)

        self._style = vtk.vtkInteractorStyleTrackballCamera()

        self._renWin = vtk.vtkRenderWindow()
        self._renWin.AddRenderer(self._ren)
        self._renWin.SetWindowName("DEM Visualizer")

        self._iren = vtk.vtkRenderWindowInteractor()
        self._iren.SetRenderWindow(self._renWin)

        self._ren.AddActor(self._axes)

    def _setupColorFunction(self, minV, maxV):

        # Create a color transfer function to be used for both the balls and arrows.
        self._colorTransferFunction = vtk.vtkColorTransferFunction()
        self._colorTransferFunction.AddRGBPoint(minV, 1.0, 0.0, 0.0)
        self._colorTransferFunction.AddRGBPoint(0.5*(minV + maxV), 0.0, 1.0, 0.0)
        self._colorTransferFunction.AddRGBPoint(maxV, 0.0, 0.0, 1.0)

    def render(self):
        """ initialize renderer """
        
        transform = vtk.vtkTransform()
        radMean, posMean = self._rad.mean() * 10, self._pos.mean(axis=0)
        posMean *= 1.2

        transform.Scale(radMean, radMean, radMean)

        self._axes.SetUserTransform(transform)
        transform.Translate(posMean)
        self._axes.SetUserTransform(transform)

        self._ren.SetBackground(0, 0, 0)
        self._camera.SetPosition(self._pos.mean(axis=0))
        self._ren.SetActiveCamera(self._camera)

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

    def attach_vel(self, vel, rad):

        velMag = norm(vel, axis=1)
        velMag /= velMag.max()

        print velMag.max()

        velMag_vtk = numpy_support.numpy_to_vtk(num_array=velMag.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        velMag_vtk.SetName("veloMag")

        vecVel = vtk.vtkFloatArray()
        vecVel.SetNumberOfComponents(3)

        for i, v in enumerate(vel):
            vecVel.InsertTuple3(i, v[0], v[1], v[2])


        #Put an arrow (vector) at each ball
        arrow = vtk.vtkArrowSource()
        arrow.SetTipRadius(rad.mean() * 0.2)
        arrow.SetShaftRadius(rad.mean() * 0.2)

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

    def attach_pos(self, pos, rad):

        self._pos = pos
        self._rad = rad

        positions_vtk = numpy_support.numpy_to_vtk(num_array=pos.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        positions_vtk.SetName("positions")

        radius_vtk = numpy_support.numpy_to_vtk(num_array=rad.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        radius_vtk.SetName("radius")

        self._points = vtk.vtkPoints()
        
        for i, r in enumerate(pos):
            self._points.InsertPoint(i, r[0], r[1], r[2])

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

        self._ren.AddActor(ballActor)

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