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

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self,None,title='Vplot Converter',size=(600,300))

        self.pars = {}

        sizer = wx.BoxSizer(wx.VERTICAL)

        panel = self.fileopen()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.format()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)
        
        panel = self.fat()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.bgcolor()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.other()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        panel = self.action()
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def fileopen(self):
        # Vplot File: _________________ [...]

        self.vpl = ''

        sizer = wx.BoxSizer(wx.HORIZONTAL)

        filetext  = wx.StaticText(self,-1,'Vplot File:')
        sizer.Add(filetext)

        self.fileentry = wx.TextCtrl(self,size=(200,-1))
        self.Bind(wx.EVT_TEXT, self.file_entry,self.fileentry)
        sizer.Add(self.fileentry,wx.EXPAND)

        butt  = wx.Button(self,-1,'Select File')
        self.Bind(wx.EVT_BUTTON,self.select_file,butt)
        sizer.Add(butt)

        return sizer

    def file_entry(self,event):
        '''Change file name using text box'''
        self.vpl = self.fileentry.GetValue()

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

        formats = list(vpconvert.pens.keys())
        formats.sort()
        nf = len(formats)
        fsizer = wx.FlexGridSizer(rows=2,cols=(nf+1)/2,hgap=5,vgap=5)
        for i in range(nf):
            fmt = formats[i]
            rb = wx.RadioButton(self,-1,fmt.upper())
            if 'png'==fmt:
                self.format = 'png'
                self.pars['format'] = self.format
                rb.SetValue(True)
            self.Bind(wx.EVT_RADIOBUTTON,self.set_format,rb)
            fsizer.Add(rb)
        sizer.Add(fsizer)

        return sizer

    def set_format(self,event):
        self.format = event.GetEventObject().GetLabel().lower()
        self.pars['format'] = self.format

    def fat(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Fat:')
        sizer.Add(text)

        slider = wx.Slider(self,-1,1,1,10,size=(250,-1),
                           style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.Bind(wx.EVT_SLIDER,self.set_fat,slider)
        self.pars['fat'] = '1'
        
        sizer.Add(slider,wx.EXPAND)

        return sizer

    def set_fat(self,event):
        self.pars['fat'] = str(event.GetEventObject().GetValue())

    def bgcolor(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Background Color:')
        sizer.Add(text)

        self.pars['bgcolor'] = 'black'
        for color in ('black','white','light','dark'):
            if 'black'==color:
                rb = wx.RadioButton(self,-1,color.capitalize(),style=wx.RB_GROUP)
                rb.SetValue(True)
            else:
                rb = wx.RadioButton(self,-1,color.capitalize())
            self.Bind(wx.EVT_RADIOBUTTON,self.set_bgcolor,rb)
            sizer.Add(rb)

        sizer.Add((0, 0), 1, wx.EXPAND)

        check = wx.CheckBox(self,-1,'Serifs')
        check.SetValue(wx.CHK_CHECKED)
        self.pars['serifs'] = 'y'
        self.Bind(wx.EVT_CHECKBOX, self.check_serifs, check)
        sizer.Add(check)

        return sizer

    def set_bgcolor(self,event):
        self.pars['bgcolor'] = event.GetEventObject().GetLabel().lower()

    def check_serifs(self,event):
        if event.IsChecked():
            self.pars['serifs'] = 'y'
        else:
            self.pars['serifs'] = 'n'

    def other(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        text  = wx.StaticText(self,-1,'Other Options:')
        sizer.Add(text)

        self.opts = ''
        entry = wx.TextCtrl(self,size=(200,-1))
        self.Bind(wx.EVT_TEXT, self.options,entry)
        sizer.Add(entry,wx.EXPAND)

        return sizer

    def options(self,event):
        self.opts = event.GetEventObject().GetValue()

    def action(self):
        # [Convert] [Quit]
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        run = wx.Button(self,-1,'Convert')
        self.Bind(wx.EVT_BUTTON,self.convert,run)
        run.SetBackgroundColour('yellow') 
        sizer.Add(run)

        sizer.Add((0, 0), 1, wx.EXPAND)

        end = wx.Button(self,-1,'Quit')
        self.Bind(wx.EVT_BUTTON,self.quit,end)
        sizer.Add(end)

        return sizer

    def convert(self,event):
        '''Convert a VPL file'''
        if not os.path.isfile(self.vpl):
            self.showerror("Can't find input file" + self.vpl)
            return

        new = '.'.join([os.path.splitext(self.vpl)[0],self.format.lower()])
        opts = ' '.join(['='.join([x, self.pars[x]]) for x in list(self.pars.keys())]) + ' ' + self.opts
    
        fail = vpconvert.convert(self.vpl,new,self.format,None,opts,False)
    
        run = "%s to %s using \"%s\"" % (self.vpl,new,opts)
        if fail:
            self.showerror("Could not convert " + run)
        else:
            self.showinfo("Converted " + run)

    def showerror(self,message):
        dlg = wx.MessageDialog(None,message,'Error',
                               wx.ICON_ERROR | wx.OK)
        dlg.ShowModal()

    def showinfo(self,message):
        dlg = wx.MessageDialog(None,message,'Message',
                               wx.ICON_INFORMATION | wx.OK)
        dlg.ShowModal()

    def quit(self,event):
        sys.exit()

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    frame.Show(True)
    app.MainLoop()


