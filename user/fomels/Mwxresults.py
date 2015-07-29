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

        sizer = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self,-1)
        buttons = self.set_buttons(panel)
        panel.SetSizer(buttons)
        buttons.Fit(self)

        sizer.Add(panel,0,wx.ALL|wx.EXPAND,5)

        self.results = commands.getoutput("scons -s results").split()
        c = self.results[-1:][0]
        if (c < 'A' or c > 'z'): 
            self.results.pop() # remove scons junk

        self.flip = {}

        panel = wx.Panel(self,-1)
        rsizer = wx.FlexGridSizer(rows=len(self.results),cols=2,vgap=5)

        for fig in self.results:
            self.flip[fig] = False

            title = self.get_title(fig)
            if title:
                title = '%s [%s]' % (fig,title)
            else:
                title = fig
            cb = wx.CheckBox(panel,-1,title)
            cb.Bind(wx.EVT_CHECKBOX, self.set_flip)
            rsizer.Add(cb,0,0)

            b = wx.Button(panel,-1,'show')
            self.Bind(wx.EVT_BUTTON,self.show(fig),b)
            rsizer.Add(b,0,0)

        panel.SetSizer(rsizer)
        rsizer.Fit(self)
        
        sizer.Add(panel,0,wx.ALL|wx.EXPAND,15)
 
        self.SetSizer(sizer)
        sizer.Fit(self)

    def set_flip(self, event):        
        sender = event.GetEventObject()
        fig = sender.GetLabel()
        self.flip[fig] = sender.GetValue()
 
    def set_buttons(self,panel):
        bsizer = wx.BoxSizer(wx.HORIZONTAL)

        bcycl = wx.Button(panel,-1,'Cycle')
        bcycl.SetBackgroundColour('light yellow') 
        self.Bind(wx.EVT_BUTTON,self.showall,bcycl)
        bsizer.Add(bcycl,0,wx.ALL,5)

        bflip = wx.Button(panel,-1,'Flip')
        bflip.SetBackgroundColour('light green') 
        self.Bind(wx.EVT_BUTTON,self.flipit,bflip)
        bsizer.Add(bflip,0,wx.ALL,5)

        bquit = wx.Button(panel,-1,'Quit')
        bquit.SetBackgroundColour('pink') 
        self.Bind(wx.EVT_BUTTON,self.quit,bquit)
        bsizer.Add(bquit,0,wx.ALL,5)

        return bsizer

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

    def get_title(self,fig):
        vpl = os.path.join('Fig',fig+'.vpl')
        if not os.path.isfile(vpl):
            return None
        txt = os.path.join('Fig',fig+'.txt')
        os.system('sfpldb < %s > %s' % (vpl,txt))
        titles = []
        try:
            plot = open(txt,'r')
            for line in plot:
                if line[0] != '[':
                    continue
                line2 = plot.next()
                if line2[:5] != 'title':
                    continue
                while line2:
                    line2 = plot.next()
                    if line2[0] == 'G':
                        title = plot.next()
                        titles.append(title.rstrip())
                        break
            title = ','.join(titles)
        except:
            title = None
        
        plot.close()
        os.unlink(txt) 
        return title
  
if __name__=='__main__':
    app = wx.App()
    frame = TestFrame()
    frame.Show()

    app.MainLoop()
