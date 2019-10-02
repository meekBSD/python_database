#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
motif finder for DNA sequence - wxPython tutorial
colored found motif in new panel.
"""
import os
import wx

APP_EXIT = 1

def find_mot_func(st, m, k):

    kmers = []
    for i in range(0,len(m)):
        if i + k <= len(m):
            m1 = m[i:i+k]
            kmers.append(m1) 

    match_profiles = [ st.index(a) for a in kmers if a in st ]
    if match_profiles == []:
        return False, None

    t0 = match_profiles[0]
    first_mismatch_t = st[t0-1 : t0 + k]
    mis_index = {}    

    n0 = m.find(first_mismatch_t[1:])
    mut_seq = first_mismatch_t[1:]

    for i in range(n0-1, 0, -1):
        mut_seq  = m[i] + mut_seq
        if mut_seq not in st:
            mut_seq = st[st.find(mut_seq[1:]) -1] + mut_seq[1:]
            mis_index[i] = st[st.find(mut_seq[1:]) -1]
        else:
            continue

    motif_result = ""
    for i in range(len(m)):
        if i in mis_index:
            motif_result += mis_index[i]
        else:
            motif_result += m[i]

        if motif_result not in st:
            break

    return motif_result, st.index(motif_result)

class Example(wx.Frame):

    def __init__(self, *args, **kwargs):
        super(Example, self).__init__(*args, **kwargs)

        self.InitUI()

    def InitUI(self):

        self.SetSize((720, 600))
        self.SetTitle('Test Seq Tools')

        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        analysis = wx.Menu()
        setting  = wx.Menu()
        help = wx.Menu()

        fileMenu.Append(101, '&Open', 'Open a new document')
        fileMenu.Append(102, '&Save', 'Save the document')
        fileMenu.AppendSeparator()

        quit= wx.MenuItem(fileMenu, APP_EXIT, '&Quit\tCtrl+Q')
        fileMenu.Append(quit)

        self.Bind(wx.EVT_MENU, self.OnQuit, id=APP_EXIT)

        analysis.Append(201, 'Find motif', 'find subsequence motif', kind=wx.ITEM_NORMAL)
        analysis.Append(202, 'GC content', 'calculation GC percentage', kind=wx.ITEM_NORMAL)
        analysis.Append(203, 'Reverse Comp', 'reverse complement seq', kind=wx.ITEM_NORMAL)

        #analysis.Append(203, 'Reverse Comp', 'reverse complement seq', kind=wx.ITEM_CHECK)

        help.Append(301, 'Usage', '', kind=wx.ITEM_NORMAL)  
        help.Append(302, 'version', '', kind=wx.ITEM_NORMAL)       

        menubar.Append(fileMenu, '&File')
        menubar.Append(analysis, '&Analysis')
        menubar.Append(setting, '&Setting')
        menubar.Append(help, '&Help')
        
        self.SetMenuBar(menubar)
        self.Centre()

        self.seq_names=wx.TextCtrl(self,style=wx.TE_MULTILINE, size = (200,500))
        self.contents=wx.TextCtrl(self,style=wx.TE_MULTILINE|wx.HSCROLL)

        #contents=wx.StaticText(frame, label="Colored text")
        #self.contents.SetForegroundColour((25,70,255)) # set text color
        self.contents.SetBackgroundColour((119,136,153)) # set text back color

        hbox = wx.BoxSizer()
        hbox.Add(self.seq_names, proportion=1,flag=wx.EXPAND)

        h1box=wx.BoxSizer(wx.HORIZONTAL)   ## VERTICAL is another option
        h1box.Add(hbox,proportion=0,flag=wx.EXPAND|wx.BOTTOM|wx.LEFT,border=5)
        h1box.Add(self.contents,proportion=1,flag=wx.EXPAND|wx.BOTTOM|wx.RIGHT|wx.TE_RICH2,border=5)

        self.SetSizer(h1box)
        
        self.Bind(wx.EVT_MENU, self.openfile, id=101)
        ## if add Button and bind Event to Button, use follow procedures
        # self.openButton=wx.Button(self, 2, label="open",pos=(225,290),size=(80,25))
        # self.openButton.Bind(wx.EVT_BUTTON, self.openfile, id = 2)  

        self.Bind(wx.EVT_MENU, self.Analyze_Motif, id=201)

    def openfile(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", os.getcwd(), "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath()
            fa_file=os.path.basename(filepath)
            fopen=open(fa_file, 'r')

            seqIds = []
            for line in fopen:
                if line.startswith(">"):
                    seqIds.append(line.rstrip()[1:])

            fopen.close()
            self.seq_names.SetValue('\n'.join(seqIds))
            self.contents.SetValue(open(fa_file).read())

    def Analyze_Motif(self, evt):
        dlg = wx.TextEntryDialog(self, 'Enter Motif Sequence', 'Motif Entry')
        dlg.SetValue('AAAGGTT')
        if dlg.ShowModal() == wx.ID_OK:
            self.findMot(dlg.GetValue())

    def findMot(self, m):
        fa_content = self.contents.GetValue()
        seqs_Dict = {}
        for i in fa_content.split('\n'):
            if i[0] == ">":
                Name = i.rstrip()[1:]
                seqs_Dict[Name] = ""
            else:
                try:
                    seqs_Dict[Name] += i.rstrip()
                except KeyError:
                    seqs_Dict[Name] = i.rstrip()
                except UnboundLocalError:
                    continue

        #self.font = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL)
        self.font1 = wx.Font(10, wx.MODERN, wx.NORMAL, wx.BOLD)

        outStr = ""
        starts_Of_Motif = []
        Desc_len = 0
        for i in seqs_Dict:
            motif_s , motif_index = find_mot_func(seqs_Dict[i], m, 9)
            if motif_s:
                outStr += i + '\n' + seqs_Dict[i] + '\n'
                starts_Of_Motif.append(Desc_len + len(i + '\n') + motif_index)
                Desc_len += len(i + '\n' + motif_s + '\n')
                
        self.contents.SetFont(self.font1)
        #seq_id = wx.StaticText(rM, 901 , i, (25,80), style = wx.ALIGN_LEFT)

        motif_len = len(m)
        self.contents.SetValue(outStr)
        for s in starts_Of_Motif:
            print(s, motif_len)
            self.contents.SetInsertionPoint(0)
            self.contents.SetStyle(s , s + motif_len + 1, wx.TextAttr("red", "blue"))

        
    def OnQuit(self, e):
        self.Close()


def main():

    app = wx.App()
    ex = Example(None)
    ex.Show()
    app.MainLoop()

    #app=wx.App()
    #frame=wx.Frame(None, title = "Seq Tools", size=(960, 750))
    #frame.Show(True)
    #app.MainLoop()

if __name__ == '__main__':
    main()

