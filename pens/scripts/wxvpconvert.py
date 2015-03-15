#!/usr/bin/env python
'''WX GUI for vpconvert'''

##   Copyright (C) 2012 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys, os

try:
    import wx
except:
    sys.stderr.write('Please install wxPython!\n\n')
    sys.exit(1)
      
import rsf.vpconvert as vpconvert

def showerror(message):
    dlg = wx.MessageDialog(None,message,
                           wx.ICON_ERROR | wx.OK)
    dlg.ShowModal()

def showinfo(message):
    dlg = wx.MessageDialog(None,message,
                           wx.ICON_INFORMATION | wx.OK)
    dlg.ShowModal()

def convert(vpl,format,opts,pars):
    '''Convert a VPL file'''
    if not os.path.isfile(vpl):
        showerror("Can't find " + vpl)
        return

    new = '.'.join([os.path.splitext(vpl)[0],format.lower()])
    opts = ' '.join(map(lambda x: '='.join([x, str(pars[x].get())]), 
                        pars.keys())) + ' ' + opts
    
    fail = vpconvert.convert(vpl,new,format,None,opts,False)
    
    run = "%s to %s using \"%s\"" % (vpl,new,opts)
    if fail:
        showerror("Could not convert " + run)
    else:
        showinfo("Converted " + run)

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self,None,title='Vplot Converter',size=(600,300))

        sizer = wx.BoxSizer(wx.VERTICAL)

        panel = self.fileopen()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.format()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)
        
        panel = self.fat()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.bgcolor()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def fileopen(self):
        # Vplot File: _________________ [...]

        sizer = wx.BoxSizer(wx.HORIZONTAL)

        filetext  = wx.StaticText(self,-1,'Vplot File:')
        self.fileentry = wx.TextCtrl(self,size=(200,-1))
        self.Bind(wx.EVT_TEXT, self.file_entry,self.fileentry)

        filebutt  = wx.Button(self,-1,'Select File')
        self.Bind(wx.EVT_BUTTON,self.select_file,filebutt)

        sizer.Add(filetext)
        sizer.Add(self.fileentry,wx.EXPAND)
        sizer.Add(filebutt)

        return sizer

    def file_entry(self,event):
        '''Change file name using text box'''
        self.vpl = self.fileentry.GetValue()
        print self.vpl

    def select_file(self,event):
        '''Select a file into entry'''
        wildcard ="Vplot files (*.vpl)|*.vpl| All files (*.*)|*.*"
        
        dialog = wx.FileDialog(None,"Choose a file",os.getcwd(),'',
                               wildcard,wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            vpl = dialog.GetPath()
            self.fileentry.SetValue(vpl)
        dialog.Destroy()

    def format(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Format:')
        sizer.Add(text)

        formats = vpconvert.pens.keys()
        formats.sort()
        nf = len(formats)
        fsizer = wx.FlexGridSizer(rows=2,cols=(nf+1)/2,hgap=5,vgap=5)
        for i in range(nf):
            fmt = formats[i]
            rb = wx.RadioButton(self,-1,fmt.upper())
            self.Bind(wx.EVT_RADIOBUTTON,self.set_format,rb)
            fsizer.Add(rb)
        self.format = formats[0]
        sizer.Add(fsizer)

        return sizer

    def set_format(self,event):
        print self.format
        self.format = event.GetEventObject().GetLabel().lower()

    def fat(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Fat:')
        sizer.Add(text)

        slider = wx.Slider(self,-1,1,1,10,size=(250,-1),
                           style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.Bind(wx.EVT_SLIDER,self.set_fat,slider)
        self.fat = 1
        
        sizer.Add(slider)

        return sizer

    def set_fat(self,event):
        self.fat = event.GetEventObject().GetValue()

    def bgcolor(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Background Color:')
        sizer.Add(text)

        self.bgcolor = 'black'
        for color in ('light','dark','black','white'):
            rb = wx.RadioButton(self,-1,color.capitalize())
            if color==self.bgcolor:
                rb.SetValue(True)
            self.Bind(wx.EVT_RADIOBUTTON,self.set_bgcolor,rb)
            sizer.Add(rb)

        return sizer

    def set_bgcolor(self,event):
        self.bgcolor = event.GetEventObject().GetLabel().lower()

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    frame.Show(True)
    app.MainLoop()


