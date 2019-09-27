'''
Created on Sep 8, 2017
@author: Andrew Abi-Mansour
'''

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

from SiDEM import MainWindow
import wx

class DerivedWindow(MainWindow):

    def onSetTraj(self, event):
        pass

    def onAddSS(self, event):
        """
        """
        if not self.component_frame:

            self.component_frame = wx.Frame(parent=self, id=-1, name='Simulation Setup', size=(self._WIDTH * 1.6, self._LENGTH * 0.6))
            frame = self.component_frame
            addEmptySpace = lambda sizer:  "h_sizer{}.Add(wx.StaticText(frame, label='', style=wx.TE_READONLY), proportion=1)".format(sizer)

            v_sizer = wx.BoxSizer(wx.VERTICAL) 
            h_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            h_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            h_sizer3 = wx.BoxSizer(wx.HORIZONTAL)

            saveBtn = wx.Button(parent=frame, id=wx.ID_ANY, label='Ok', size=(80,40))
            resetBtn = wx.Button(parent=frame, id=wx.ID_ANY, label='Cancel', size=(80,40))

            materials = ['Glass']
            eval(addEmptySpace(1))
            h_sizer1.Add(wx.StaticText(frame, label="Materials: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            model = wx.ComboBox(parent=frame, choices = materials)
            h_sizer1.Add(model, proportion=4)

            eval(addEmptySpace(1))
            h_sizer1.Add(wx.StaticText(frame, label="Id: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            ids = wx.TextCtrl(frame)
            h_sizer1.Add(ids, proportion=4)

            eval(addEmptySpace(1))
            h_sizer1.Add(wx.StaticText(frame, label="# Particles: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            simBox = wx.TextCtrl(frame)
            h_sizer1.Add(simBox, proportion=4)
            eval(addEmptySpace(1))

            eval(addEmptySpace(2))
            radiusType = ['Constant', 'Gaussian Number', 'Uniform Number', 'Uniform Mass']
            h_sizer2.Add(wx.StaticText(frame, label="Radius: ", style=wx.TE_READONLY), proportion=1)
            #eval(addEmptySpace(2))
            radType = wx.ComboBox(parent=frame, choices=radiusType)
            h_sizer2.Add(radType, proportion=1)
            radArg1 = wx.TextCtrl(frame)
            h_sizer2.Add(radArg1, proportion=1)
            radArg2 = wx.TextCtrl(frame)
            h_sizer2.Add(radArg2, proportion=1)
            eval(addEmptySpace(2))

            eval(addEmptySpace(3))
            h_sizer3.Add(item=saveBtn, proportion=1)
            eval(addEmptySpace(3))
            h_sizer3.Add(item=resetBtn, proportion=1)

            # Add some empty space to separate h_sizer1 from the upper border
            v_sizer.Add(wx.StaticText(frame, label="", style=wx.TE_READONLY), proportion=1)
            v_sizer.Add(item=h_sizer1, proportion=1)
            v_sizer.Add(item=h_sizer2, proportion=1)
            v_sizer.Add(item=h_sizer3, proportion=1)

            frame.SetSizer(v_sizer)

            frame.Centre()
            frame.Show()
        else:
            self.UpdateDisplayPanel('Simulation settings window already open')

    def onSimSetup(self, event):
        """
        """

        if not self.settings_frame:
            def settToolBar(frame):
                """
                """
                ToolBar = frame.CreateToolBar(style=wx.DEFAULT_DIALOG_STYLE)
                
                loadBtn = ToolBar.AddLabelTool(wx.ID_ANY, '&Trajectory', wx.Bitmap(self.iconsDir + 'Settings.png'))
                ToolBar.AddSeparator()

                addSSBtn = ToolBar.AddLabelTool(wx.ID_ANY, 'Add &Subsystem', wx.Bitmap(self.iconsDir + 'AddSS.png'))
                ToolBar.AddSeparator()

                addWall = ToolBar.AddLabelTool(wx.ID_ANY, 'Add &Wall', wx.Bitmap(self.iconsDir + 'AddWall.png'))

                ToolBar.AddSeparator()
                
                ToolBar.Realize()
                frame.Bind(wx.EVT_TOOL, self.onSetTraj, loadBtn)
                frame.Bind(wx.EVT_TOOL, self.onAddSS, addSSBtn)

            self.settings_frame = wx.Frame(parent=self, id=-1, name='Simulation Setup', size=(self._WIDTH * 1.6, self._LENGTH * 1.1))
            frame = self.settings_frame
            addEmptySpace = lambda sizer:  "h_sizer{}.Add(wx.StaticText(frame, label='', style=wx.TE_READONLY), proportion=1)".format(sizer)

            settToolBar(frame)

            v_sizer = wx.BoxSizer(wx.VERTICAL) 
            h_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            h_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            h_sizer3 = wx.BoxSizer(wx.HORIZONTAL)
            h_sizer4 = wx.BoxSizer(wx.HORIZONTAL)

            saveBtn = wx.Button(parent=frame, id=wx.ID_ANY, label='Save', size=(80,40))
            resetBtn = wx.Button(parent=frame, id=wx.ID_ANY, label='Reset', size=(80,40))

            contactModels = ['Spring Dashpot', 'Hertz Mindlin', 'Hysteresis Thorn']
            eval(addEmptySpace(1))
            h_sizer1.Add(wx.StaticText(frame, label="Model: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            model = wx.ComboBox(parent=frame, choices = contactModels)
            h_sizer1.Add(model, proportion=4)

            eval(addEmptySpace(1))
            shapes = ['Spherical']
            h_sizer1.Add(wx.StaticText(frame, label="Shape: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            shape = wx.ComboBox(parent=frame, choices = shapes)
            h_sizer1.Add(shape, proportion=4)

            eval(addEmptySpace(1))
            h_sizer1.Add(wx.StaticText(frame, label="Box Size: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(1))
            simBox = wx.TextCtrl(frame)
            h_sizer1.Add(simBox, proportion=4)
            eval(addEmptySpace(1))

            eval(addEmptySpace(2))
            h_sizer2.Add(wx.StaticText(frame, label="Ext. Force: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(2))
            force = wx.TextCtrl(frame)
            h_sizer2.Add(force, proportion=4)

            eval(addEmptySpace(2))
            h_sizer2.Add(wx.StaticText(frame, label="Timestep: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(2))
            timeStep = wx.TextCtrl(frame)
            h_sizer2.Add(timeStep, proportion=4)

            eval(addEmptySpace(2))
            unitSys = ['SI', 'Metal', 'LJ']
            h_sizer2.Add(wx.StaticText(frame, label="Units: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(2))
            model = wx.ComboBox(parent=frame, choices = unitSys)
            h_sizer2.Add(model, proportion=4)
            eval(addEmptySpace(2))

            eval(addEmptySpace(3))
            h_sizer3.Add(wx.StaticText(frame, label="Newton 3rd Law: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(3))
            newton = wx.CheckBox(parent=frame, id=-1)
            wx.EVT_CHECKBOX(frame, newton.GetId(), self.CheckForHETATM)
            h_sizer3.Add(newton, proportion=1)

            eval(addEmptySpace(3))
            h_sizer3.Add(wx.StaticText(frame, label="Cohesion: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(3))
            cohesion = wx.CheckBox(parent=frame, id=-1)
            wx.EVT_CHECKBOX(frame, cohesion.GetId(), self.CheckForHETATM)
            h_sizer3.Add(cohesion, proportion=1)

            eval(addEmptySpace(3))
            h_sizer3.Add(wx.StaticText(frame, label="Rolling: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(3))
            rolling = wx.CheckBox(parent=frame, id=-1)
            wx.EVT_CHECKBOX(frame, rolling.GetId(), self.CheckForHETATM)
            h_sizer3.Add(rolling, proportion=1)

            eval(addEmptySpace(3))
            h_sizer3.Add(wx.StaticText(frame, label="Tangential: ", style=wx.TE_READONLY), proportion=2)
            eval(addEmptySpace(3))
            tangential = wx.CheckBox(parent=frame, id=-1)
            wx.EVT_CHECKBOX(frame, tangential.GetId(), self.CheckForHETATM)
            h_sizer3.Add(tangential, proportion=1)
            eval(addEmptySpace(3))

            eval(addEmptySpace(4))
            h_sizer4.Add(item=saveBtn, proportion=1)
            eval(addEmptySpace(4))
            h_sizer4.Add(item=resetBtn, proportion=1)
            
            eval(addEmptySpace(4))
            boundaryType = ['Periodic', 'Fixed', 'Shrink-wrapped']
            h_sizer4.Add(wx.StaticText(frame, label="Boundaries: ", style=wx.TE_READONLY), proportion=2)
            #eval(addEmptySpace(4))
            boundaryX = wx.ComboBox(parent=frame, choices = boundaryType)
            h_sizer4.Add(boundaryX, proportion=2)
            boundaryY = wx.ComboBox(parent=frame, choices = boundaryType)
            h_sizer4.Add(boundaryY, proportion=2)
            boundaryZ = wx.ComboBox(parent=frame, choices = boundaryType)
            h_sizer4.Add(boundaryZ, proportion=2)
            eval(addEmptySpace(4))

            # Add some empty space to separate h_sizer1 from the upper border
            v_sizer.Add(wx.StaticText(frame, label="", style=wx.TE_READONLY), proportion=1)
            v_sizer.Add(item=h_sizer1, proportion=1)
            v_sizer.Add(item=h_sizer2, proportion=1)
            v_sizer.Add(item=h_sizer3, proportion=1)
            v_sizer.Add(item=h_sizer4, proportion=1)

            frame.SetSizer(v_sizer)

            frame.Centre()
            frame.Show()
        else:
            self.UpdateDisplayPanel('Simulation settings window already open')