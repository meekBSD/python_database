import wx


class MyMenu(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(680, 500))
        menubar = wx.MenuBar()

        first_f = wx.Menu()
        analysis = wx.Menu()
        setting  = wx.Menu()
        help = wx.Menu()

        first_f.Append(101, '&Open', 'Open a new document')
        first_f.Append(102, '&Save', 'Save the document')
        first_f.AppendSeparator()

        quit = wx.MenuItem(first_f, 105, '&Quit\tCtrl+Q', 'Quit the Application')
        quit_image = wx.Bitmap('exit-1.png').ConvertToImage()
        quit_image.Rescale(50,50)
        #quit.SetBitmap(wx.Image('exit-1.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        quit.SetBitmap(wx.Bitmap(quit_image))
        first_f.Append(quit)

        analysis.Append(201, 'Find motif', 'find subsequence motif', kind=wx.ITEM_NORMAL)
        analysis.Append(202, 'GC content', 'calculation GC percentage', kind=wx.ITEM_NORMAL)
        analysis.Append(203, 'Reverse Comp', 'reverse complement seq', kind=wx.ITEM_NORMAL)

        #analysis.Append(203, 'Reverse Comp', 'reverse complement seq', kind=wx.ITEM_CHECK)

        help.Append(301, 'Usage', '', kind=wx.ITEM_NORMAL)  
        help.Append(302, 'version', '', kind=wx.ITEM_NORMAL)       

        menubar.Append(first_f, '&File')
        menubar.Append(analysis, '&Analysis')
        menubar.Append(setting, '&Setting')
        menubar.Append(help, '&Help')
        self.SetMenuBar(menubar)
        self.Centre()
        self.Bind(wx.EVT_MENU, self.OnQuit, id=105)

    def OnQuit(self, event):
        self.Close()

class MyApp(wx.App):
    def OnInit(self):
        frame = MyMenu(None, -1, 'menu_test.py')
        frame.Show(True)
        return True


app = MyApp(0)
app.MainLoop()



