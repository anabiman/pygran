import windows
import wx
import Settings.Language as Lang

if __name__ == '__main__':

        app = wx.App(redirect=False)

        win = windows.DerivedWindow(parent=None,name=Lang.__SofName__, title=Lang.__SofName__)
        win.Show()

        app.MainLoop()