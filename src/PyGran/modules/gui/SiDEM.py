#!/usr/bin/python
# -*- coding: utf8 -*- 

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

'''
Created on Sep 8, 2017
@author: Andrew Abi-Mansour
'''

import wx, os, glob
import subprocess
import urllib2
import wx.lib.agw.multidirdialog as MDD
from PyGran import Simulator
from importlib import import_module

#import Plot
import sys
from mpi4py import MPI
import visualize
import Settings.Language as Lang

class MainWindow(wx.Frame):

 

    def __init__(self, parent, title, name):
        """
    	This is the "constructor". It initializes the main frame and (so far) the two panels: upper panel and display pannel. 
    	"""
        self._LENGTH = 650
        self._WIDTH = 800
        self.iconsDir = os.path.dirname(os.path.realpath(__file__)) + '/Icons/'

        __slots__ = ["display_panel", "display_box", "LoadedPanel", "LoadedBox", "panel_box", "contents_box", \
            "upper_box1", "upper_box2", "upper_box3", "ScriptFileTxt", "GenBtn", "ClearBtn", "SettBtn"\
            "Loaded_PDB", "Generated_PDB", "tmp_pdb_file", "command_txt", "read_pdb_options"]

        self.__engines__ = vars(Simulator.engines).keys()

        self.BGColor = wx.Colour(95,95,95)

        # Declare the children frames
        self.settings_frame = None
        self.component_frame = None

        # Loaded script files
        self.loadedScript = None
        self.loadedScriptPy = None

        # Initialize some GUI stuff
        super(MainWindow,self).__init__(parent=parent, name=name, title=title,  pos=wx.DefaultPosition, size=(self._WIDTH, self._LENGTH))

        # A panel for reading user input / commands
        self.InputPanel()

        # A panel for loaded scripts
        # self.LoadedPanel()

        # Setup language 
        self.MainMenu('English')

        self.ToolBar()
        self.NoteBook = wx.Notebook(parent=self, style=wx.NB_TOP, id=wx.ID_ANY)
        self.Page_Output = wx.NotebookPage(parent=self.NoteBook, id=wx.ID_ANY)
        self.NoteBook.AddPage(self.Page_Output, "Output")
        
        self.Page_AtomAnalysis = wx.NotebookPage(parent=self.NoteBook, id=wx.ID_ANY)
        self.NoteBook.AddPage(self.Page_AtomAnalysis, "Particle Analysis")
        self.Page_AtomAnalysis.SetBackgroundColour(self.BGColor)
        
        self.Page_CGAnalysis = wx.NotebookPage(parent=self.NoteBook, id=wx.ID_ANY)
        self.NoteBook.AddPage(self.Page_CGAnalysis, "Contact Force")
        self.Page_CGAnalysis.SetBackgroundColour(self.BGColor)
        
        self.panel_box = wx.BoxSizer(wx.VERTICAL)
        self.panel_box.Add(item=self.NoteBook, proportion=10, flag=wx.EXPAND)
        self.panel_box.Add(item=self.InputPanel, proportion=1, flag=wx.EXPAND)
        #self.panel_box.Add(item=self.LoadedPanel, proportion=1, flag=wx.EXPAND)

        self.DisplayPanel()
        
        self.SetAutoLayout(True)
        self.SetSizer(self.panel_box)
        
        self.Layout()
        self.Centre()
        
        # Initialize some I/O stuff
        self.loadedVars = {}
        self.CWD = os.getcwd()
        self.read_pdb_options = {'Check_for_HETATM':True, 'Debug':False, 'T_number':5} 
        # must read this from the user input
        
        #wx.FutureCall(1000, self.ShowMessage)
        self.Font = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas', encoding=wx.FONTENCODING_ISO8859_6)
        self.nProcs = 1

        self.UpdateDisplayPanel("Using a total of {} procs".format(self.nProcs))
        self.UpdateDisplayPanel('Awaiting input ...')
        
    def __del__(self):
        try:
            os.system('rm {}'.format(self.tmp_pdb_file))
            self.UpdateDisplayPanel('Cleaned up junk files.')
        except:
            pass
        
    def MainMenu(self, language):
        """
        Undocumented
        """
        
        # Create the main menu bar here. All menus will be appended to this bar.
        self.MainMenuBar = wx.MenuBar(style=wx.DEFAULT_STATUSBAR_STYLE)
        
        # start with the file menu: Open, Save, Import
        self.FileMenu = wx.Menu(style=wx.DEFAULT_STATUSBAR_STYLE)
        
        self.OpenBtn = self.FileMenu.Append(id=wx.ID_OPEN, text="&Open\tCtrl-O", help="Open a script file")
        self.SaveBtn = self.FileMenu.Append(id=wx.ID_CLOSE, text="&Save\tCtrl-S", help="Save output to file")
        self.CWDBtn = self.FileMenu.Append(id=wx.ID_CLOSE, text="&Change Dir \tCtrl-D", help="Change current working directory")
        
        self.Bind(event=wx.EVT_MENU, handler=self.OnOpen, source=self.OpenBtn)
        self.Bind(event=wx.EVT_MENU, handler=self.OnSave, source=self.SaveBtn)
        self.Bind(event=wx.EVT_MENU, handler=self.OnSelDir, source=self.CWDBtn)
        
        # Create import submenu
        self.FileMenu.AppendSeparator() # haven o clue what this does
        self.ImpSubMenu = wx.Menu()
        
        self.impMeshBtn = self.ImpSubMenu.Append(wx.ID_ANY, 'Mesh (STL)')
        self.Bind(event=wx.EVT_MENU, handler=self.onImpMeshBtn, source=self.impMeshBtn)
        
        self.impTrajBtn = self.ImpSubMenu.Append(wx.ID_ANY, 'Trajectory')
        self.Bind(event=wx.EVT_MENU, handler=self.onImpTrajBtn, source=self.impTrajBtn)

        self.ImpPDBBtn = self.ImpSubMenu.Append(wx.ID_ANY, 'PDB from RCSB')
        self.Bind(event=wx.EVT_MENU, handler=self.OnRCSB, source=self.ImpPDBBtn)
        
        self.FileMenu.AppendMenu(wx.ID_ANY, '&Import', self.ImpSubMenu)
        
        # Help menu goes here
        self.HelpMenu = wx.Menu(style=wx.DEFAULT_STATUSBAR_STYLE)
        
        self.AboutBtn = self.HelpMenu.Append(id=wx.ID_ABOUT, text="&About\tCtrl-A", help="About this program")
        self.Bind(event=wx.EVT_MENU, handler=self.onHelp, source=self.AboutBtn)
        
        # Configure menu goes here
        self.ConfMenu = wx.Menu(style=wx.DEFAULT_STATUSBAR_STYLE)
        
        self.ConfBtn = self.ConfMenu.Append(id=wx.ID_SETUP, text="&Simulation\tCtrl-S", help="Setup Simulation")
        self.Bind(event=wx.EVT_MENU, handler=self.onSimSetup, source=self.ConfBtn)
        
        self.setComputBtn = self.ConfMenu.Append(id=wx.ID_ANY, text='&Computation\tCtrl-C', help='')
        self.Bind(event=wx.EVT_MENU, handler=self.onSetProcs, source=self.setComputBtn)

        self.setProcsBtn = self.ConfMenu.Append(id=wx.ID_ANY, text='&Hardware\tCtrl-H', help='Set number of cores')
        self.Bind(event=wx.EVT_MENU, handler=self.onSetProcs, source=self.setProcsBtn)
        
        self.setEngBtn = self.ConfMenu.Append(id=wx.ID_ANY, text='&Engine\tCtrl-E', help='Set number of cores')
        self.Bind(event=wx.EVT_MENU, handler=self.onSetEng, source=self.setEngBtn)

        # Analysis menu goes here
        self.AnaMenu = wx.Menu(style=wx.DEFAULT_STATUSBAR_STYLE)
        
        self.AtomAnBtn = self.AnaMenu.Append(id=wx.ID_ANY, text="&Atomistic\tAlt-Ctrl-A", help="Begin atomistic analysis")
        self.CGAnBtn = self.AnaMenu.Append(id=wx.ID_ANY, text="&Coarse-Grained\tAlt-Ctrl-C", help="Begin macromolecular analysis")
        
        self.Bind(event=wx.EVT_MENU, handler=self.OnAtomAnBtn, source=self.AtomAnBtn)
        self.Bind(event=wx.EVT_MENU, handler=self.onSimSetup, source=self.CGAnBtn)
        
        
        # Setup language
        self.LanguageClass = eval('Lang.{}()'.format(language))

        self.MainMenuBar.Append(menu=self.FileMenu, title=self.LanguageClass.File)
        self.MainMenuBar.Append(menu=self.ConfMenu, title=self.LanguageClass.Configure)
        self.MainMenuBar.Append(menu=self.AnaMenu, title=self.LanguageClass.Analysis)
        self.MainMenuBar.Append(menu=self.HelpMenu, title=self.LanguageClass.Help)
        
        self.SetMenuBar(menubar=self.MainMenuBar)
        
    def ToolBar(self):
        """
        """
        iconsDir = self.iconsDir

        self.ToolBar = self.CreateToolBar(style=wx.DEFAULT_DIALOG_STYLE)
        
        self.OpenBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Open', wx.Bitmap(iconsDir+'Load.png'))
        self.SaveBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Save', wx.Bitmap(iconsDir+'Save.png'))
        
        self.ToolBar.AddSeparator()
        
        self.runBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Start', wx.Bitmap(iconsDir+'Generate.png'))
        self.visBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Plot', wx.Bitmap(iconsDir+'Plot.png'))
        
        self.ToolBar.AddSeparator()
        
        self.ClrBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Clear', wx.Bitmap(iconsDir+'Clear.png'))
        self.QuitBtn = self.ToolBar.AddLabelTool(wx.ID_ANY, '&Quit', wx.Bitmap(iconsDir+'Quit.png'))
        
        self.ToolBar.Realize()
        self.Bind(wx.EVT_TOOL, self.OnGen, self.runBtn)
        self.Bind(wx.EVT_TOOL, self.clearDisplayPanel, self.ClrBtn)
        self.Bind(wx.EVT_TOOL, self.OnOpen, self.OpenBtn)
        self.Bind(wx.EVT_TOOL, self.OnSave, self.SaveBtn)
        self.Bind(wx.EVT_TOOL, self.onVisualize, self.visBtn)
        self.Bind(wx.EVT_TOOL, self.OnQuit, self.QuitBtn)
        
    def UpdateDisplayPanel(self, txt, dtype=None):
        """
        This function updates the content text box in the display panel.
        """
        if dtype is None:
            self.contents_txt.AppendText('{}\n'.format(txt))
        else:
            self.contents_txt.AppendText('{} {}\n'.format(txt, dtype))

        self.Show()
        
    def clearDisplayPanel(self, event = None):
        """
        This function clears the display panel from all output.
        """
        self.contents_txt.Clear()
        self.Show()
            
    def DisplayPanel(self):
        """
        This function initiates the display panel, which is responsible for displaying any relevant info
        to the user. The actual mechanics of this function is yet to be decided.
        """
        #self.display_panel  = wx.Panel(self, -1, style=wx.TAB_TRAVERSAL|wx.BORDER, name='display_panel')
        #self.display_panel.SetBackgroundColour("GRAY")
        self.display_box = wx.BoxSizer(wx.VERTICAL)
        
        self.contents_txt = wx.TextCtrl(parent=self.Page_Output, size=(-1, -1),style = wx.TE_MULTILINE | wx.TE_READONLY)
        self.contents_txt.SetBackgroundColour(self.BGColor)
        #self.contents_txt.SetFont(self.Font)
        self.contents_txt.SetForegroundColour(wx.BLACK)
        
        self.display_box.Add(self.contents_txt, proportion=1, flag=wx.EXPAND)
        
        self.Page_Output.SetSizer(self.display_box)
    
    def OnTerminal_txt(self, event):

        process = subprocess.Popen([self.command_txt.Value], stdout=subprocess.PIPE, shell=True)
        invoked_command = process.communicate()[0].split('\n')

        for s in invoked_command:
            self.UpdateDisplayPanel(s)
        
    def LoadedPanel(self):

        self.LoadedPanel  = wx.Panel(self, -1, style=wx.TAB_TRAVERSAL|wx.BORDER, name='Loaded Panel')  
        self.LoadedBox = wx.BoxSizer(wx.HORIZONTAL)
   
        # upper_box1 text stuff
        self.ScriptFileTxt = wx.TextCtrl(self.LoadedPanel, value="Script file:", \
                                           style = wx.TE_READONLY | wx.TRANSPARENT_WINDOW | wx.ALIGN_CENTER_HORIZONTAL)
        
        self.LoadedScript = wx.TextCtrl(self.LoadedPanel, value = "None", style = wx.ALIGN_CENTER_HORIZONTAL)
        #self.loaded_file_txt.SetBackgroundColour("Gray")
        
       
        # Expansion set for the buttons here, i.e., they will expand both vertically and horizontally
        self.LoadedBox.Add(item=self.ScriptFileTxt, proportion=1, flag=wx.EXPAND)
        self.LoadedBox.Add(item=self.LoadedScript, proportion=5, flag=wx.EXPAND)
        
        
        self.LoadedPanel.SetSizer(self.LoadedBox)

    def InputPanel(self):

        self.InputPanel  = wx.Panel(self, -1, style=wx.TAB_TRAVERSAL|wx.BORDER, name='Input Panel')  
        self.InputBox = wx.BoxSizer(wx.HORIZONTAL)
   
        self.InputName = wx.TextCtrl(self.InputPanel, value="Input command:", \
                                           style =  wx.TE_READONLY | wx.TRANSPARENT_WINDOW | wx.ALIGN_CENTER_HORIZONTAL)

        self.InputTxt = wx.TextCtrl(self.InputPanel, value = "", style = wx.ALIGN_CENTER_HORIZONTAL | wx.TE_PROCESS_ENTER)

        self.InputTxt.SetBackgroundColour((15,15,15))
        self.InputTxt.SetForegroundColour(wx.WHITE)
       
        # Expansion set for the buttons here, i.e., they will expand both vertically and horizontally
        self.InputBox.Add(item=self.InputName, proportion=2, flag=wx.EXPAND)
        self.InputBox.Add(item=self.InputTxt, proportion=5, flag=wx.EXPAND)
        
        self.InputPanel.SetSizer(self.InputBox)

        self.Bind(wx.EVT_TEXT_ENTER, self.onReadCmd, self.InputTxt)

    def onReadCmd(self, event):
        """ reads text from InputBox, passes it asa command to Liggghts, updates the display panel,
        and finally clears InputBox"""

        command = self.InputTxt.GetValue()
        
        if len(command.split()) == 1:
            if command == 'clc':
                self.clearDisplayPanel()

            elif command == 'whos':
                for item in self.loadedVars:
                    self.UpdateDisplayPanel(item, type(self.loadedVars[item]))
            elif command == 'visualize':
                self.onVisualize()

            else:
                self.UpdateDisplayPanel('Unknown command')

        else:
            method, var = command.split(' ', 1)

            if method == 'plot':
                if var not in self.loadedVars:
                    self.UpdateDisplayPanel('Error: could not find {}'.format(var))
                else:
                    self.UpdateDisplayPanel('Generating plot ...')

                    try:
                        if var == 'mesh':
                            visualize.visualize(meshFname=self.loadedVars[var])
                        else:
                            if var == 'mesh and particles':
                                if var not in self.loadedVars:
                                    self.UpdateDisplayPanel('Error: no particles loaded')
                                else:
                                    try:
                                        visualize.visualize(meshFname=self.loadedVars['mesh'], dumpFname=self.loadedVars['mesh and particles'])
                                    except:
                                        self.UpdateDisplayPanel('Plot failed with unexpected error: {}'.format(sys.exc_info()[0]))
                            else:
                                pass # make 2D plot
                    except:
                        pass

            elif method == 'run':
                if self.__selEngine__:
                    # try to call a DEM-engine method based on user-supplied cmds
                    try:
                        self._module.command('{}'.format(var))
                    except:
                        self.UpdateDisplayPanel('Unexpected error: {}'.format(sys.exc_info()[0]))
                else:
                    self.UpdateDisplayPanel('No engine selected. Make sure an available DEM engine is installed.')
            elif method == 'unix':
                try:
                    output = os.popen(var).read()
                    self.UpdateDisplayPanel(output)
                except:
                    self.UpdateDisplayPanel('Unexpected error: {}'.format(sys.exc_info()[0]))
            else:   
                self.UpdateDisplayPanel('Unknown command')

        self.InputTxt.Clear()

    def onHelp(self, event):
        """
    	This is not yet documented.
    	"""
        
        description = self.LanguageClass.description
        
        license = """SiDEM is free software; you can redistribute 
it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation; either version 2 of the License, 
or (at your option) any later version.

SiDEM is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details. You should have 
received a copy of the GNU General Public License along with SiDEM; 
if not, write to the Free Software Foundation, Inc., 59 Temple Place, 
Suite 330, Boston, MA  02111-1307  USA"""

        info = wx.AboutDialogInfo()
        info.SetIcon(wx.Icon(self.iconsDir + 'CAAM.png', wx.BITMAP_TYPE_PNG))
        info.SetName('SiDEM')
        info.SetVersion('1.0')
        info.SetDescription(description)
        info.SetCopyright('(C) 2016 Andrew Abi Mansour')
        info.SetWebSite('http://www.chemical-computing.com')
        info.SetLicence(license)
        info.AddDeveloper('Andrew Abi Mansour')
        info.AddDocWriter('Andrew Abi Mansour')
        info.AddArtist('https://www.iconfinder.com/')
        info.AddTranslator('Andrew Abi Mansour')
        
        wx.AboutBox(info)
            
    def OnQuit(self, event):
        """
        Exit application (main loop killed)
        """
        self.Close()
        
    def onImpMeshBtn(self, event):
        dlg = wx.FileDialog(self, message="Open a mesh STL file...", defaultDir=self.CWD, defaultFile="*.stl", style=wx.OPEN)
        dlg.CentreOnParent()

        if dlg.ShowModal() == wx.ID_OK:
            self.loadedVars['mesh'] = dlg.GetPath()
            self.UpdateDisplayPanel('Loaded mesh from {}'.format(self.loadedVars['mesh']))

        dlg.Destroy()

    def onImpTrajBtn(self, event):
        dlg = wx.FileDialog(self, message="Open a mesh STL file...", defaultDir=self.CWD, defaultFile="*.stl", style=wx.OPEN)
        dlg.CentreOnParent()

        if dlg.ShowModal() == wx.ID_OK:
            self.loadedVars['mesh and particles'] = dlg.GetPath()
            self.UpdateDisplayPanel('Loaded trajectory from {}'.format(self.loadedVars['mesh and particles']))

        dlg.Destroy()

    def OnOpen(self, event):
        """
        This function Launches an instance of a new file dialog that allows the user
		to select a specific input file.
        """
        wildcard = "Python (*.py)|*.py|"     \
           "LIGGGHTS (*.in)|*.in|" \
           "All files (*.*)|*.*"
           
        dlg = wx.FileDialog(self, message="Open a file...", defaultDir=self.CWD, defaultFile="", \
                            style=wx.FD_FILE_MUST_EXIST | wx.FD_OPEN | wx.FD_PREVIEW, wildcard=wildcard)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.loadedScript = dlg.GetPath()
            
            if(self.loadedScript.split('.')[1] == 'py'):
                
                self.UpdateDisplayPanel('Loaded input script file: {}'.format(self.loadedScript))
                
            else:
                self.UpdateDisplayPanel('Input must be a valid python/LIGGGHTS file.')
        dlg.Destroy()
        
    def OnSave(self, event):
        """
        This function Launches an instance of a new file dialog that allows the user
                to save the (currently) generated pdb file.
        """
        dlg = wx.FileDialog(self, message="Save a file...", defaultDir=self.CWD, defaultFile="", \
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            try:
                cdir = dlg.GetDirectory()
                pdb_file_name = dlg.GetFilename()
                
                fp = open(pdb_file_name,'w')
                
                self.UpdateDisplayPanel('Attempting to save file {} ...'.format(pdb_file_name))
                fp.write(self.Generated_PDB)
                
                self.UpdateDisplayPanel('Success! Coordinates file {} written to {}'.format(pdb_file_name,cdir))
                fp.close()
                
            except:
                raise ValueError('Could not save output file {} to {}'.format(pdb_file_name,cdir))
        dlg.Destroy()

    def onVisualize(self, event):

        self.UpdateDisplayPanel('Launching ovito')
        cmd = 'ovito'

        if self.loadedScriptPy:
            pDict = self.loadedScriptPy.pDict
            traj = pDict['traj']['output'] + '/' + pDict['traj']['dir'] + '/' + pDict['traj']['file']
            cmd = cmd + ' ' + traj
            self.UpdateDisplayPanel('Loading trajectory file {}'.format(traj))
        try:
            os.popen(cmd).read()
        except:
            self.UpdateDisplayPanel('No visualization software found. Make sure ovito is properly installed on your system.')

    def OnAA(self, event):
        Analysis = event.GetSelection()
        diag = wx.Dialog(self, size=(800,600), title=Analysis, style=wx.DIALOG_NO_PARENT)
        GUIPanel = Plot.GUIPanel(parent=diag)
        
        diag.ShowModal()
        
    def OnGen(self, event):
        """
        This is the main computational engine that runs a DEM simulation.
        Only works in python 2.x
        For python 3.x must use: for line in popen.stdout: print(line.decode(), end='')
        """
        
        if self.loadedScript:
            try:
                split = self.loadedScript.split('/')
                wdir, script = '/'.join(split[:-1]), split[-1]

                cmd = 'python {}'.format(script)
                self.UpdateDisplayPanel('Changing current working dir to ' + wdir)
                os.chdir(wdir)
                self.CWD = wdir
                self.UpdateDisplayPanel('Running simulation: ' + cmd)

                output = os.popen(cmd).read()
                #for output in self.execute(cmd):
                self.UpdateDisplayPanel(output)
                self.UpdateDisplayPanel('Simulation done. Importing DEM params ...')

                self.loadedScriptPy = __import__(script.split('.py')[0])
                # ver hackish!!!
                self.loadedScriptPy.pDict['traj']['output'] = wdir + '/' + sorted(glob.glob('out-*'), key=os.path.getmtime)[-1]

            except:
                raise
        else:
            self.UpdateDisplayPanel('Load input script file first.')

    def execute(self, cmd):
        """ A function that executed commands using 'yield' to output any line from the command-line
        to the screen """
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    def OnGetIndices(self, event):
        dlg = wx.Dialog(parent=self,title="I am modal, close me first to get to main frame")
        dlg.ShowModal()

    def OnHelpOlig(self, event):
        """
        """
        dlg = wx.MessageDialog(self, "This modules generates the oligmer indices based on the oligomer number. For instance, when this"
	                                 " number is set to 5, a pentamer is generated (5 repeating subunits), when it's set to 6, a hexamer is generated, and so on. \n", "Help", wx.OK)
        dlg.ShowModal()
        dlg.Destroy() #destroy dialog when finished
            
    def ProgressState(self, max_val):
        """
        Copied from the wxPython demo
        """
 
        dlg = wx.ProgressDialog("Progress dialog example",
                               "Generating full capsid .. ",
                               maximum = max_val,
                               parent=self,
                               style = wx.PD_CAN_ABORT
                                | wx.PD_APP_MODAL
                                | wx.PD_ELAPSED_TIME
                                | wx.PD_ESTIMATED_TIME
                                | wx.PD_REMAINING_TIME
                                )   
 
        return dlg
        
    def OnSelDir(self, event):
        """
        The effect of this function must be taken into account for all other variables
        """
        dlg = wx.DirDialog(parent=self, message='Choose working directory', defaultPath=self.CWD, \
                           style = wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON | wx.DD_DIR_MUST_EXIST)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.CWD = dlg.GetPath()
            self.UpdateDisplayPanel('Changed current working directory to {}'.format(self.CWD))
                                               
        dlg.Destroy()
               
    def OnAtomAnBtn(self, event):
        """
        """
        dlg = wx.Dialog(parent=self, title="Atomic Analysis", style=wx.DEFAULT_DIALOG_STYLE | wx.BORDER_SUNKEN)
        dlg.ShowModal()
      
    def onSetProcs(self, event):
        dlg = wx.TextEntryDialog(self, 'Enter number of CPU cores:','Set CPU config')

        if dlg.ShowModal() == wx.ID_OK:
            self.nProcs = eval(dlg.GetValue())

            if type(self.nProcs) is not int:
                self.UpdateDisplayPanel('Error. Number of processors must be an integer.')
                self.nProcs = 1
            else:
                self.UpdateDisplayPanel('Number of procs set to {}'.format(self.nProcs))

        dlg.Destroy()

    def onSetEng(self, event):

        dlg = wx.SingleChoiceDialog(self, 'Available Engines', 'Select DEM engine', self.__engines__, wx.CHOICEDLG_STYLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.__selEngine__ = dlg.GetStringSelection()
            self.UpdateDisplayPanel('Importing PyGran.Simulator.engine.' + self.__selEngine__)
            self._module = import_module('PyGran.Simulator.engine_' + self.__selEngine__)

        dlg.Destroy()

    def OnRCSB(self, event):
        dlg = wx.TextEntryDialog(self, 'Enter PDB ID:','Download PDB file from RCSB')
        
        if dlg.ShowModal() == wx.ID_OK:
            self.PDB_ID = dlg.GetValue()
            self.url='http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={}'.format(self.PDB_ID)
            self.UpdateDisplayPanel('Fetching file from http://www.RCSB.org ... ')
            
            try:
                webFile = urllib2.urlopen(url=self.url)
                Success = True
            except:
                self.UpdateDisplayPanel('Failed to download file. File does not exist or connection timed out.')
                Success = False
                
            # Report this to the developer in case of breakage
            if Success:
                try:
                    self.Loaded_PDB_file = webFile.read()
                    self.Loaded_PDB = '{} (RCSB)'.format(self.PDB_ID)
                    webFile.close()
                    self.UpdateDisplayPanel('Loaded PDB ID {}.'.format(self.PDB_ID))
                    self.loaded_file_txt.ChangeValue(self.Loaded_PDB)
                except:
                    raise
                
                try:
                    self.Natoms, self.Operator = Transform.ReadPDB(self.Loaded_PDB_file)
                    self.UpdateDisplayPanel('Read {} atoms.'.format(self.Natoms))
                except:
                    self.UpdateDisplayPanel('PDB file does not contain BIOMT matrices')
                
        dlg.Destroy()
                
    def CustomDiag(self, event):
        """
        """
        font = wx.Font(11, wx.MODERN, wx.NORMAL, wx.NORMAL)
        
        self.mesh_dlg = wx.Dialog(size=(400,80), parent=self,title="Import mesh", style= wx.DEFAULT_DIALOG_STYLE | wx.BORDER_SUNKEN)
        v_sizer = wx.BoxSizer(wx.VERTICAL) 
        h_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        PDB_ID_RE_txt = wx.TextCtrl(parent=self.RCSB_dlg, style=wx.TE_READONLY, value='PDB ID:')
        self.PDB_ID_txt = wx.TextCtrl(parent=self.RCSB_dlg)
        
        SubmitBtn = wx.Button(parent=self.RCSB_dlg, id=wx.ID_ANY, label='Submit')
        CancelBtn = wx.Button(parent=self.RCSB_dlg, id=wx.ID_ANY, label='Cancel')
        
        self.RCSB_dlg.Bind(event=wx.EVT_BUTTON, handler=self.OnRCSBSubmit, source=SubmitBtn)
        self.PDB_ID_txt.Bind(wx.EVT_KEY_DOWN, self.OnRCSBSubmit)
        
        self.RCSB_dlg.Bind(event=wx.EVT_BUTTON, handler=self.OnRCSBCancel, source=CancelBtn)
        self.PDB_ID_txt.Bind(wx.EVT_KEY_DOWN, self.OnRCSBCancel)
        
        h_sizer.Add(item=PDB_ID_RE_txt, proportion=2, flag=wx.EXPAND)
        h_sizer.Add(item=self.PDB_ID_txt, proportion=4, flag=wx.EXPAND)
        h_sizer.Add(item=SubmitBtn, proportion=2, flag=wx.EXPAND)
        h_sizer.Add(item=CancelBtn, proportion=2, flag=wx.EXPAND)
        
        v_sizer.Add(item=h_sizer, proportion=1, flag=wx.EXPAND)
        
        self.RCSB_dlg.SetSizer(v_sizer)
        self.PDB_ID_txt.SetFocus()
        
        self.PDB_ID_txt.SetFont(font)
        PDB_ID_RE_txt.SetFont(font)
        SubmitBtn.SetFont(font)
        CancelBtn.SetFont(font)
        
        self.RCSB_dlg.ShowModal()
        
    def OnRCSBCancel(self, event):
        try:
            if event.GetKeyCode() == wx.WXK_ESCAPE:
                self.RCSB_dlg.Destroy()
            else:
                event.Skip()
        except:
            event.Skip()
                
        if event.GetEventType() == wx.wxEVT_COMMAND_BUTTON_CLICKED:
            self.RCSB_dlg.Destroy()
        else:
            event.Skip()
        
    def OnRCSBSubmit(self, event):
        """
        This function is not used but was kept for demonstration purposes only
        """
        Invoke = False
        
        try:
            if event.GetKeyCode() == wx.WXK_RETURN:
                Invoke = True
        except:
            event.Skip()
        
        if event.GetEventType() == wx.wxEVT_COMMAND_BUTTON_CLICKED:
            Invoke = True
        else:
            event.Skip()
            
        if Invoke:
            self.PDB_ID = self.PDB_ID_txt.GetValue()
            self.RCSB_dlg.Destroy()
            
            self.url='http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={}'.format(self.PDB_ID)
            self.UpdateDisplayPanel('Fetching file from http://www.RCSB.org ... ')
            
            try:
                webFile = urllib2.urlopen(url=self.url)
                self.Loaded_PDB = '{} (RCSB)'.format(self.PDB_ID)
                Success = True
            except:
                self.UpdateDisplayPanel('Failed to download file. File does not exist or connection timed out.')
                Success = False
                
            # Report this to the developer in case of breakage
            if Success:
                try:
                    self.Loaded_PDB_file = webFile.read()
                    webFile.close()
                    
                    self.UpdateDisplayPanel('Loaded PDB ID {}.'.format(self.PDB_ID))
                except:
                    raise
                
                try:
                    self.Natoms, self.Operator = Transform.ReadPDB(self.Loaded_PDB_file)
                except:
                    self.UpdateDisplayPanel('PDB file does not contain BIOMT matrices.')
        
    def OnOlig(self, event):
        """
    	Undocumented
    	"""
        frame_tmp = wx.Frame(parent=self, id=-1, name='Oligomer indices generator')
        panel_tmp = wx.Panel(parent=frame_tmp, id=-1, style=wx.TAB_TRAVERSAL|wx.BORDER, name='Oligomer indices generator panel')
    
        IndicesBtn = wx.Button(parent=panel_tmp, id=-1, label='Get &Indices')
        SelectBtn = wx.Button(parent=panel_tmp, id=-1, label='&Select PDB')
        HelpBtn = wx.Button(parent=panel_tmp, id=-1, label='&Help')
    
        frame_tmp.Bind(wx.EVT_BUTTON, self.OnGetIndices, IndicesBtn)
        frame_tmp.Bind(wx.EVT_BUTTON, self.OnGetIndices, SelectBtn)
        frame_tmp.Bind(wx.EVT_BUTTON, self.OnHelpOlig, HelpBtn)
    
        oligo_txt = wx.TextCtrl(panel_tmp, value="Oligomer number:", \
                                      style = wx.TE_READONLY | wx.ALIGN_CENTER_HORIZONTAL)
        
        olig_input = wx.TextCtrl(panel_tmp, value = "", style = wx.ALIGN_CENTER_HORIZONTAL)
        
        VertBox = wx.BoxSizer(wx.VERTICAL)
        HorzBox1 = wx.BoxSizer(wx.HORIZONTAL)
        HorzBox2 = wx.BoxSizer(wx.HORIZONTAL)
    
        HorzBox1.Add(IndicesBtn, proportion=1, flag=wx.EXPAND)
        HorzBox1.Add(SelectBtn, proportion=1, flag=wx.EXPAND)
        HorzBox1.Add(HelpBtn, proportion=1, flag=wx.EXPAND)
    
        HorzBox2.Add(oligo_txt, proportion=1, flag=wx.EXPAND, border=wx.BORDER_DOUBLE)
        HorzBox2.Add(olig_input, proportion=1, flag=wx.EXPAND, border=wx.BORDER_DOUBLE)
        
        VertBox.Add(HorzBox1, flag=wx.EXPAND)
        VertBox.Add(HorzBox2, flag=wx.EXPAND)
    
        panel_tmp.SetSizer(VertBox)	
    
        frame_tmp.Centre()
        frame_tmp.Show()
        
    def CheckForHETATM(self, event):
        """
        Undocumented (yet)
        """
        if self.cb.GetValue():
            self.read_pdb_options['Check_for_HETATM'] = True
            self.UpdateDisplayPanel('Hetero atoms will now be taken into account when reading pdb files(s).') 
        else:
            self.read_pdb_options['Check_for_HETATM'] = False
            self.UpdateDisplayPanel('Hetero atoms will NOT be taken into account when reading pdb file(s).')

