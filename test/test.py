'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

import wx
import vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
 
class p1(wx.Panel):
    def __init__(self,parent):
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
            
    def renderthis(self):
            # open a window and create a renderer
            ren = vtk.vtkRenderer()
            self.widget.GetRenderWindow().AddRenderer(ren)
   
            # to generate polygonal data for a cone.
            cone = vtk.vtkConeSource()
            cone.SetResolution(25)
 
            # o take the polygonal data from the vtkConeSource and
            # create a rendering for the renderer.
            coneMapper = vtk.vtkPolyDataMapper()
            coneMapper.SetInputData(cone.GetOutput())
 
            # create an actor for our scene
            coneActor = vtk.vtkActor()
            coneActor.SetMapper(coneMapper)
            # Add actor
            ren.AddActor(coneActor)
 
            axes = vtk.vtkAxesActor()
            self.marker = vtk.vtkOrientationMarkerWidget()
            self.marker.SetInteractor( self.widget._Iren )
            self.marker.SetOrientationMarker( axes )
            self.marker.SetViewport(0.75,0,1,0.25)
            self.marker.SetEnabled(1)
 
            ren.ResetCamera()
            ren.ResetCameraClippingRange()
            cam = ren.GetActiveCamera()
            cam.Elevation(10)
            cam.Azimuth(70)
            self.isploted = True
 
class TestFrame(wx.Frame):
    def __init__(self,parent,title):
        wx.Frame.__init__(self,parent,title=title,size=(650,600), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|
                  wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)
        self.sp = wx.SplitterWindow(self)
        self.p1 = p1(self.sp)
        self.p2 = wx.Panel(self.sp,style=wx.SUNKEN_BORDER)
         
        self.sp.SplitHorizontally(self.p1,self.p2,470)
 
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText("Click on the Plot Button")
         
        self.plotbut = wx.Button(self.p2,-1,"plot", size=(40,20),pos=(10,10))
        self.plotbut.Bind(wx.EVT_BUTTON,self.plot)
         
 
    def plot(self,event):
        if not self.p1.isploted:
            self.p1.renderthis()
            self.statusbar.SetStatusText("Use your mouse to interact with the model")
 
         
app = wx.App(redirect=False)
frame = TestFrame(None,"PyGran Visualizer")
frame.Show()
app.MainLoop()
