#!/usr/bin/env python
'Show data with zoom'

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
import sys, os, subprocess, atexit, tempfile
import rsf.prog

try:
    import wx
except:
    sys.stderr.write('Please install wxPython!\n\n')
    sys.exit(1)

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s file.rsf <sfgrey options>\n\n' % sys.argv[0])
    sys.exit(2)

inp = sys.argv[1]

bindir   = os.path.join(rsf.prog.RSFROOT,'bin')
sfgrey   = os.path.join(bindir,'sfgrey')
sfget    = os.path.join(bindir,'sfget')
sfwindow = os.path.join(bindir,'sfwindow')
ppmpen   = os.path.join(bindir,'ppmpen')
sfrm     = os.path.join(bindir,'sfrm')

ppmfiles = []
rsffiles = []

width = 1024
height = 768

x0 = 170
x1 = 937

y0 = 96
y1 = 671

def rsf2image(rsf):
    global ppmfiles
    ppm = os.path.splitext(rsf)[0]+'.ppm'
    command = '< %s %s %s | %s > %s' % (rsf,sfgrey,' '.join(sys.argv[2:]),ppmpen,ppm)
    if os.system(command) or not os.path.isfile(ppm):
        sys.stderr.write('Failed to execute "%s"\n\n' % command)
        sys.exit(3)
    img = PhotoImage(file=ppm)
    if not ppm in ppmfiles:
        ppmfiles.append(ppm)
    return img

class Canvas(wx.Window):
    def __init__(self,parent,ID):
        wx.Window.__init__(self,parent,ID,size=(width,height))
        self.SetBackgroundColour('black')
        self.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

        image = self.rsf2image(inp)
        self.image = image.ConvertToBitmap()
        self.color = 'Yellow'
        self.thickness = 1
        self.pen = wx.Pen(self.color,self.thickness,wx.SOLID)
        self.brush = wx.Brush(self.color,wx.TRANSPARENT)
        self.Bind(wx.EVT_PAINT,self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN,self.OnLeftDown)
        self.Bind(wx.EVT_LEFT_UP,self.OnLeftUp)
        self.Bind(wx.EVT_MOTION,self.OnMotion)
    def rsf2image(self,rsf):
        global ppmfiles
        ppm = os.path.splitext(rsf)[0]+'.ppm'
        command = '< %s %s %s | %s > %s' % (rsf,sfgrey,' '.join(sys.argv[2:]),ppmpen,ppm)
        if os.system(command) or not os.path.isfile(ppm):
            sys.stderr.write('Failed to execute "%s"\n\n' % command)
            sys.exit(3)
        if not ppm in ppmfiles:
            ppmfiles.append(ppm)
        img = wx.Image(ppm)
        return img
    def OnPaint(self,event):
        self.dc = wx.PaintDC(self)
        brush = wx.Brush('black')
        self.dc.SetBackground(brush)
        self.dc.Clear() # clear with background brush
        self.dc.DrawBitmap(self.image,0,0,True)
    def OnLeftDown(self,event):
        self.start = event.GetPositionTuple()
        self.width = 0
        self.height = 0
        self.CaptureMouse()
    def OnLeftUp(self,event):
        if self.HasCapture():
            self.ReleaseMouse()
    def OnMotion(self,event):
        if event.Dragging() and event.LeftIsDown():
            self.drawSquare(event)
        event.Skip()
    def drawSquare(self,event):
        pos = event.GetPositionTuple()
        self.width = pos[0]-self.start[0]
        self.height = pos[1]-self.start[1]
        self.OnPaint(event)
        self.dc.SetBrush(self.brush)
        self.dc.SetPen(self.pen)
        self.dc.DrawRectangle(self.start[0], self.start[1], 
                              self.width, self.height)

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self,None,title=inp,size=(width+1,height+26))

        self.sketch = Canvas(self,-1)
        self.sketch.Bind(wx.EVT_MOTION,self.OnSketchMotion)

        self.status = self.CreateStatusBar()
        self.status.SetFieldsCount(1)

    def OnSketchMotion(self,event):
        self.status.SetStatusText('Pos: ' + str(event.GetPositionTuple()),0)
        event.Skip()

@atexit.register
def cleanup():
    global ppmfiles,rsffiles,sfrm
    for ppm in filter(os.path.isfile,ppmfiles):
        os.unlink(ppm)
    for rsf in filter(os.path.isfile,rsffiles):
        os.system(' '.join([sfrm,rsf]))

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    frame.Show(True)
    app.MainLoop()


