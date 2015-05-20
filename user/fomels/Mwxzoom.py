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

def hist(inp,func,var,default=None):
    command = '< %s %s %s parform=n' % (inp,sfget,var)
    devnull = open(os.devnull,"w")
    pipe = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=devnull,shell=True)
    val = pipe.stdout.read().rstrip()
    if val:
        return func(val)
    else:
        return default

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

        self.SetFrame(inp)

        self.label1 = hist(inp,str,'label1','Y')
        self.unit1  = hist(inp,str,'unit1','')

        self.label2 = hist(inp,str,'label2','X')
        self.unit2  = hist(inp,str,'unit2','')

    def SetFrame(self,inp):
        n1 = hist(inp,int,'n1',1)
        d1 = hist(inp,float,'d1',1.0)
        self.o1 = hist(inp,float,'o1',0.0)
        self.yscale = (n1-1)*d1/(y1-y0)
        
        n2 = hist(inp,int,'n2',1)
        d2 = hist(inp,float,'d2',1.0)
        self.o2 = hist(inp,float,'o2',0.0)
        self.xscale = (n2-1)*d2/(x1-x0)
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
        x,y = event.GetPositionTuple()
        self.x = min(max(x,x0),x1)
        self.y = min(max(y,y0),y1)
        self.width = 0
        self.height = 0
        self.CaptureMouse()
    def OnLeftUp(self,event):
        if self.HasCapture():
            self.ReleaseMouse()
            x,y = event.GetPositionTuple()
            x = min(max(x,x0),x1)
            y = min(max(y,y0),y1)
            self.Redraw(x,y)
            self.OnPaint(event)
    def Redraw(self,x,y):
        x1 = self.o2+(self.x-x0)*self.xscale
        x2 = self.o2+(     x-x0)*self.xscale
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        y1 = self.o1+(self.y-y0)*self.yscale
        y2 = self.o1+(     y-y0)*self.yscale
        ymin = min(y1,y2)
        ymax = max(y1,y2)
        inp2 = tempfile.mktemp(suffix='.rsf')
        command = '< %s %s min1=%g max1=%g min2=%g max2=%g > %s' % (inp,sfwindow,ymin,ymax,xmin,xmax,inp2)
        if os.system(command) or not os.path.isfile(inp2):
            sys.stderr.write('Failed to execute "%s"\n\n' % command)
            sys.exit(4)
        rsffiles.append(inp2)
        self.SetFrame(inp2)
        image = self.rsf2image(inp2)
        self.image = image.ConvertToBitmap()
    def OnMotion(self,event):
        if event.Dragging() and event.LeftIsDown():
            self.drawSquare(event)
        event.Skip()
    def drawSquare(self,event):
        x,y = event.GetPositionTuple()
        x = min(max(x,x0),x1)
        y = min(max(y,y0),y1)
        self.width  = x-self.x
        self.height = y-self.y
        self.OnPaint(event)
        self.dc.SetBrush(self.brush)
        self.dc.SetPen(self.pen)
        self.dc.DrawRectangle(self.x, self.y, 
                              self.width, self.height)
    def Position(self,x,y):
        if x >= x0 and y >= y0 and x <= x1 and y <= y1:
            x = self.o2+(x-x0)*self.xscale
            y = self.o1+(y-y0)*self.yscale
            return '(%s = %g %s, %s = %g %s)' % (self.label1,y,self.unit1,
                                                 self.label2,x,self.unit2)
        else:
            return ''
        
class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self,None,title=inp,size=(width+1,height+26))

        panel = wx.Panel(self, wx.ID_ANY)
        panel.Bind(wx.EVT_CHAR, self.OnKeyDown)
        panel.SetFocus()

        self.sketch = Canvas(panel,-1)
        self.sketch.Bind(wx.EVT_MOTION,self.OnSketchMotion)

        self.status = self.CreateStatusBar()
        self.status.SetFieldsCount(1)
    def OnSketchMotion(self,event):
        x,y = event.GetPositionTuple()
        position = self.sketch.Position(x,y)
        self.status.SetStatusText(position,0)
        event.Skip()
    def OnKeyDown(self,event):
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_ESCAPE:
            print 'Escape'
        elif chr(keycode) == 'q':
            sys.exit(0)
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


