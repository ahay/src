#!/usr/bin/env python
'Explore project results'

##   Copyright (C) 2010 University of Texas at Austin
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
import sys

try:
    import wx
except:
    sys.stderr.write('Please install wx!\n\n')
    sys.exit(1)
    
from signal import *
import os, commands

class TestFrame(wx.Frame):
    def __init__(self):
        cwd = os.getcwd()
        proj = os.path.join(
            os.path.basename(os.path.dirname(os.path.dirname(cwd))),
            os.path.basename(os.path.dirname(cwd)),
            os.path.basename(cwd))
        wx.Frame.__init__(self,None,-1,proj)
        self.pid = 0


        signal(SIGINT,self.handler)

        panel = wx.Panel(self,-1)
        buttons = self.set_buttons(panel)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttons)

        self.results = commands.getoutput("scons -s results").split()
        c = self.results[-1:][0]
        if (c < 'A' or c > 'z'): 
            self.results.pop() # remove scons junk
        length = max(map(len,self.results))

        self.flip = {}

        for fig in self.results:
            self.flip[fig] = False
            row = wx.BoxSizer(wx.HORIZONTAL)

            cb = wx.CheckBox(panel,-1,fig)
            cb.Bind(wx.EVT_CHECKBOX, self.set_flip)
            row.Add(cb,0,wx.ALL|wx.EXPAND,5)

            b = wx.Button(panel,-1,'show')
            self.Bind(wx.EVT_BUTTON,self.show(fig),b)
            row.Add(b,0,wx.ALL|wx.EXPAND,5)

            sizer.Add(row,0,wx.ALL|wx.EXPAND,5)

        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)
        self.SetSizer(sizer)
        sizer.Fit(self)

    def set_flip(self, event):        
        sender = event.GetEventObject()
        fig = sender.GetLabel()
        self.flip[fig] = sender.GetValue()
 
    def set_buttons(self,panel):
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        bcycl = wx.Button(panel,-1,'Cycle')
        bcycl.SetBackgroundColour('light yellow') 
        self.Bind(wx.EVT_BUTTON,self.showall,bcycl)
        sizer.Add(bcycl,0,wx.ALL|wx.EXPAND,5)

        bflip = wx.Button(panel,-1,'Flip')
        bflip.SetBackgroundColour('light green') 
        self.Bind(wx.EVT_BUTTON,self.flipit,bflip)
        sizer.Add(bflip,0,wx.ALL|wx.EXPAND,5)

        bquit = wx.Button(panel,-1,'Quit')
        bquit.SetBackgroundColour('pink') 
        self.Bind(wx.EVT_BUTTON,self.quit,bquit)
        sizer.Add(bquit,0,wx.ALL|wx.EXPAND,5)

        return sizer

    def quit(self,event):
        sys.exit(0)

    def handler(self,signum,frame):
        sys.stdout.write("\n%s: aborting...\n" % sys.argv[0])
        if self.pid:
            os.kill(self.pid,SIGINT)
        sys.exit(1)

    def showall(self,event):
        self.pid = os.fork()
        if not self.pid:
            for fig in self.results:
                os.system ("scons %s.view" % fig)
        sys.exit()

    def flipit(self,event):
        figs = map(lambda x: 'Fig/%s.vpl' % x,
                   filter(lambda y: self.flip[y],self.results))
        if figs:
            self.pid = os.fork()
            if not self.pid:
                command = "sfpen %s" % ' '.join(figs)
                sys.stderr.write(command+'\n')
                os.system (command)
                sys.exit()

    def show(self,fig):
        def showfig(event):
            os.system("scons %s.view" % fig)
        return showfig

if __name__=='__main__':
    app = wx.App()
    frame = TestFrame()
    frame.Show()

    app.MainLoop()
