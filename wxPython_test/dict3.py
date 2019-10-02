#!/usr/bin/python
#-*- coding: utf8 -*-

import wx
from urllib import request
from bs4 import BeautifulSoup


app=wx.App()
win=wx.Frame(None,title="YouDao",size=(410,335))


bkg =wx.Panel(win)

def main(evt):
        q = search_item.GetValue()
        page = request.urlopen('http://dict.youdao.com/search?le=eng&q='+ q)
        trans_soup = BeautifulSoup(page, features="html.parser")
        try:
            basic_list = trans_soup.find('div', attrs={'class': 'trans-container'}).findNext('ul').findAll('li')
            #translation=[unicode(item)[4:-5] for item in basic_list]
            translation=[str(item)[4:-5] for item in basic_list]
            result.SetValue("\n".join(translation))
        except:
            result.SetValue('未找到相关翻译。')
                  

loadButton=wx.Button(bkg,label='Search')
saveButton=wx.Button(bkg,label='Save')

loadButton.Bind(wx.EVT_BUTTON,main)


search_item=wx.TextCtrl(bkg)
result=wx.TextCtrl(bkg,style=wx.TE_MULTILINE|wx.HSCROLL)

hbox = wx.BoxSizer()
hbox.Add(search_item,proportion=1,flag=wx.EXPAND)
hbox.Add(loadButton,proportion=0,flag=wx.LEFT,border=5)
hbox.Add(saveButton,proportion=0,flag=wx.LEFT,border=5)

vbox=wx.BoxSizer(wx.VERTICAL)
vbox.Add(hbox,proportion=0,flag=wx.EXPAND|wx.ALL,border=5)
vbox.Add(result,proportion=1,flag=wx.EXPAND|wx.LEFT|wx.BOTTOM|wx.RIGHT,border=5)

bkg.SetSizer(vbox)

win.Show()
app.MainLoop()
