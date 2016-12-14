# adapted from:
# http://wiki.wxpython.org/Getting%20Started
# http://www.cs.colorado.edu/~kena/classes/5448/s11/presentations/pearse.pdf

import wx
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import pylab as pl

class GUIPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.parent = parent

        # create some sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        h_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # A button
        self.button1 = wx.Button(self, label="Buttton1")
        self.CloseBtn = wx.Button(self, label="Close")
        
        self.Bind(wx.EVT_BUTTON, self.OnClick,self.button1)
        self.Bind(event=wx.EVT_BUTTON, handler=self.OnClose, source=self.CloseBtn)

        # put up a figure
        self.figure = matplotlib.figure.Figure()
        self.axes = self.drawplot(self.figure)
        self.canvas = FigureCanvas(self, -1, self.figure)

        h_sizer.Add(item=self.button1)
        h_sizer.Add(item=self.CloseBtn)
        
        sizer.Add(self.canvas, 0, wx.ALIGN_CENTER | wx.ALL | wx.EXPAND)
        sizer.Add(h_sizer, 0, wx.ALIGN_CENTER | wx.ALL | wx.EXPAND)

        self.SetSizerAndFit(sizer)
         
    def log(self, fmt, *args):
        print (fmt % args)
        
    def OnClick(self,event):
        self.log("button clicked, id#%d\n", event.GetId())
        
    def OnClose(self, event):
        self.parent.Destroy()
        
    def drawplot(self, fig):
        ax = fig.add_subplot(1,1,1)
        t = pl.arange(0,1,0.001)
        ax.plot(t,t*t)
        ax.grid()
        return ax
