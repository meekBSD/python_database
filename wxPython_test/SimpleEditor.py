import wx

app=wx.App()
win=wx.Frame(None,title="Simple Editor",size=(410,335))


bkg =wx.Panel(win)

def openfile(evt):
    file=filename.GetValue()
    fopen=open(file)
    contents.SetValue(fopen.read())
    fopen.close()

def savefile(evt):
    filepath=filename.GetValue()
    fopen=open(filepath,'w')
    fopen.write(contents.GetValue())
    fopen.close()

loadButton=wx.Button(bkg,label='Open')
saveButton=wx.Button(bkg,label='Save')

loadButton.Bind(wx.EVT_BUTTON,openfile)
saveButton.Bind(wx.EVT_BUTTON,savefile)

filename=wx.TextCtrl(bkg)
contents=wx.TextCtrl(bkg,style=wx.TE_MULTILINE|wx.HSCROLL)

hbox = wx.BoxSizer()
hbox.Add(filename,proportion=1,flag=wx.EXPAND)
hbox.Add(loadButton,proportion=0,flag=wx.LEFT,border=5)
hbox.Add(saveButton,proportion=0,flag=wx.LEFT,border=5)

vbox=wx.BoxSizer(wx.VERTICAL)
vbox.Add(hbox,proportion=0,flag=wx.EXPAND|wx.ALL,border=5)
vbox.Add(contents,proportion=1,flag=wx.EXPAND|wx.LEFT|wx.BOTTOM|wx.RIGHT,border=5)

bkg.SetSizer(vbox)

win.Show()
app.MainLoop()
