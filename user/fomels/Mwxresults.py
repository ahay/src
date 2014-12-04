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
        panel = self.set_buttons()

        self.results = commands.getoutput("scons -s results").split()
        c = self.results[-1:][0]
        if (c < 'A' or c > 'z'): 
            self.results.pop() # remove scons junk
        length = max(map(len,self.results))

        self.flip = {}
        r=50

        for fig in self.results:
            self.flip[fig] = False
            cb = wx.CheckBox(panel,-1,fig,(35,r))
            cb.Bind(wx.EVT_CHECKBOX, self.set_flip)
            r += 30

#    Button(frame,text=fig,cursor='hand2',command=show(fig),width=length).grid(row=r,column=1)

    def set_flip(self, event):        
        sender = event.GetEventObject()
        fig = sender.GetLabel()
        self.flip[fig] = sender.GetValue()
        print self.flip

    def set_buttons(self):
        panel = wx.Panel(self,-1)

        bquit = wx.Button(panel,-1,'Quit',pos=(10,10))
        bquit.SetBackgroundColour('pink') 
        self.Bind(wx.EVT_BUTTON,self.quit,bquit)

        bflip = wx.Button(panel,-1,'Flip',pos=(100,10))
        bflip.SetBackgroundColour('light green') 
        self.Bind(wx.EVT_BUTTON,self.flipit,bflip)

        bcycl = wx.Button(panel,-1,'Cycle',pos=(190,10))
        bcycl.SetBackgroundColour('light yellow') 
        self.Bind(wx.EVT_BUTTON,self.showall,bcycl)

        return panel

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
                   filter(lambda y: self.flip[y].get(),self.results))
        if figs:
            self.pid = os.fork()
            if not self.pid:
                command = "sfpen %s" % ' '.join(figs)
                sys.stderr.write(command+'\n')
                os.system (command)
                sys.exit()

    def show(self,fig):
        def showfig():
            os.system("scons %s.view" % fig)
        return showfig

if __name__=='__main__':
    app = wx.App()
    frame = TestFrame()
    frame.Show()

    app.MainLoop()
