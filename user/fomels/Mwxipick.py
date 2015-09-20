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
from rsf.prog import RSFROOT

try:
    import wx
except:
    sys.stderr.write('Please install wxPython!\n\n')
    sys.exit(1)

if os.isatty(sys.stdin.fileno()):
    sys.stderr.write('Usage: %s < file.rsf [sfgrey/ppmpen options] > picks.txt\n\n' % sys.argv[0])
    sys.exit(2)

byte = tempfile.mktemp(suffix='.rsf')
name = os.path.splitext(byte)[0]

args = ' '.join(sys.argv[1:])

ppms = []

bindir   = os.path.join(RSFROOT,'bin')
sfwindow = os.path.join(bindir,'sfwindow')
sfbyte   = os.path.join(bindir,'sfbyte')
sfgrey   = os.path.join(bindir,'sfgrey')
sfget    = os.path.join(bindir,'sfget')
sfrm     = os.path.join(bindir,'sfrm')
ppmpen   = os.path.join(bindir,'ppmpen')

command = '%s %s > %s' % (sfbyte,args,byte)
p = subprocess.Popen(command,shell=True)
res = p.communicate()

if not os.path.isfile(byte):
    sys.stderr.write('Failed to execute "%s"\n\n' % command)
    sys.exit(2)

def hist(func,var,default=None):
    global byte,sfget
    command = '< %s %s %s parform=n' % (byte,sfget,var)
    devnull = open(os.devnull,"w")
    pipe = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=devnull,shell=True)
    val = pipe.stdout.read().rstrip()
    if val:
        return func(val)
    else:
        return default

width = 1024
height = 768

x0 = 170
x1 = 937

y0 = 96
y1 = 671

r = 5 # circle radius

n3 = hist(int,'n3',1)
i3 = 0

picks = []
npick = 0

for i in range(n3):
    picks.append({})

class Canvas(wx.Window):
    def __init__(self,parent,ID):
        wx.Window.__init__(self,parent,ID,size=(width,height))
        self.SetBackgroundColour('black')
        self.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

        image = self.rsf2image(0)
        self.image = image.ConvertToBitmap()
        self.thickness = 1
        pencolor='White'
        color='Yellow'
        self.pen = wx.Pen(pencolor,self.thickness,wx.SOLID)
        self.brush = wx.Brush(color)
        self.black = wx.Brush('black')
        self.Bind(wx.EVT_PAINT,self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN,self.OnClick)
        self.Bind(wx.EVT_RIGHT_DOWN,self.Delete)
        
        self.SetFrame()

        self.label1 = hist(str,'label1','Y')
        self.unit1  = hist(str,'unit1','')

        self.label2 = hist(str,'label2','X')
        self.unit2  = hist(str,'unit2','')
    def rsf2image(self,i3):
        global byte
        ppm = '%s%d.ppm' % (name,i3)
        if not os.path.isfile(ppm):
            command = '< %s %s n3=1 f3=%d | %s %s | %s %s > %s' % \
                (byte,sfwindow,i3,sfgrey,args,ppmpen,args,ppm)
            if os.system(command) or not os.path.isfile(ppm):
                sys.stderr.write('Failed to execute "%s"\n\n' % command)
                sys.exit(3)
            ppms.append(ppm)
        img = wx.Image(ppm)
        return img
    def OnPaint(self,event):
        self.dc = wx.PaintDC(self)
        self.dc.SetBackground(self.black)
        self.dc.Clear() # clear with background brush
        self.dc.DrawBitmap(self.image,0,0,True)
    def Redraw(self,i3):
        image = self.rsf2image(i3)
        self.image = image.ConvertToBitmap()
        self.dc.SetBackground(self.black)
        self.dc.DrawBitmap(self.image,0,0,True)
        self.DrawPicks(i3)
    def SetFrame(self):
        n1 = hist(int,'n1',1)
        d1 = hist(float,'d1',1.0)
        self.o1 = hist(float,'o1',0.0)
        self.yscale = (n1-1)*d1/(y1-y0)
        
        n2 = hist(int,'n2',1)
        d2 = hist(float,'d2',1.0)
        self.o2 = hist(float,'o2',0.0)
        self.xscale = (n2-1)*d2/(x1-x0)
    def Position(self,x,y):
        if x >= x0 and y >= y0 and x <= x1 and y <= y1:
            x = self.o2+(x-x0)*self.xscale
            y = self.o1+(y-y0)*self.yscale
            return '(%s = %g %s, %s = %g %s)' % (self.label1,y,self.unit1,
                                                 self.label2,x,self.unit2)
        else:
            return ''
    def ScalePick(self,x,y):
        xs = self.o2+(x-x0)*self.xscale
        ys = self.o1+(y-y0)*self.yscale
        return (ys,xs,i3)
    def UnscalePick(self,pick):
        ys = pick[0]
        xs = pick[1]
        x = x0+(xs-self.o2)/self.xscale
        y = y0+(ys-self.o1)/self.yscale
        return (x,y)
    def DrawPicks(self,i3):
        for pick in picks[i3].values():
            (x,y) = self.UnscalePick(pick[0])
            self.dc.SetPen(self.pen)
            self.dc.SetBrush(pick[1])
            self.dc.DrawCircle(x,y,r)
    def SetColor(self,color):
        self.brush = wx.Brush(color)
    def OnClick(self,event):
        global npick, r
        self.dc.SetPen(self.pen)
        self.dc.SetBrush(self.brush)
        x,y = event.GetPositionTuple()
        if x >= x0 and y >= y0 and x <= x1 and y <= y1:
            npick += 1
            tag = 'pick%d' % npick
            self.dc.DrawCircle(x,y,r)
            #canvas.tag_bind(tag,'<ButtonPress-2>',selectpick)
            #canvas.tag_bind(tag,'<B2-Motion>',movepick)
            #canvas.tag_bind(tag,'<ButtonRelease-2>',movedpick)
            picks[i3][tag]=[self.ScalePick(x,y),self.brush]
        event.Skip()
    def Delete(self,event):
        x,y = event.GetPositionTuple()
        for p in picks[i3].keys():
            xp,yp = self.UnscalePick(picks[i3][p][0])
            if (x-xp)*(x-xp)+(y-yp)*(y-yp) < 4:
                del picks[i3][p]
                self.Redraw(i3)
                break
        event.Skip()
        
class MainFrame(wx.Frame):
    def __init__(self):
        global n3
        
        wx.Frame.__init__(self,None,title='Interactive Picking',size=(width+1,height+26))

        sizer = wx.BoxSizer(wx.VERTICAL)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        if n3 > 1:
            self.bnext = wx.Button(self,-1,'Next >')
            self.Bind(wx.EVT_BUTTON,self.NextFrame,self.bnext)
            buttons.Add(self.bnext,0,wx.ALL,5)
            self.bprev = wx.Button(self,-1,'< Prev')
            self.bprev.Disable()
            self.Bind(wx.EVT_BUTTON,self.PrevFrame,self.bprev)
            buttons.Add(self.bprev,0,wx.ALL,5)
            self.label = wx.StaticText(self,-1,'%d of %d' % (i3+1,n3))
            buttons.Add(self.label,0,wx.CENTER,5)
        buttons.AddStretchSpacer()
        colorpick = wx.Button(self,-1,'Change Color')
        colorpick.SetBackgroundColour('light yellow') 
        self.Bind(wx.EVT_BUTTON,self.PickColor,colorpick)
        buttons.Add(colorpick,0,wx.ALL,5)
        bquit = wx.Button(self,-1,'Quit')
        bquit.SetBackgroundColour('pink') 
        self.Bind(wx.EVT_BUTTON,self.Quit,bquit)
        buttons.Add(bquit,0,wx.ALL,5)
        sizer.Add(buttons,flag=wx.EXPAND)

        panel = wx.Panel(self, wx.ID_ANY)
        panel.Bind(wx.EVT_CHAR, self.OnKeyDown)
        panel.SetFocus()

        self.sketch = Canvas(panel,-1)
        self.sketch.Bind(wx.EVT_MOTION,self.OnSketchMotion)

        sizer.Add(panel,wx.EXPAND)
        self.SetSizer(sizer)

        self.status = self.CreateStatusBar()
        self.status.SetFieldsCount(1)
    def NextFrame(self,event):
        global i3,n3,canvas,next,prev,picks
        i3 += 1
        if i3 == n3-1:
            self.bnext.Disable()
        self.bprev.Enable()
        self.label.SetLabel('%d of %d' % (i3+1,n3))
        self.sketch.Redraw(i3)
    def PrevFrame(self,event):
        global i3,n3,canvas,next,prev,picks
        i3 -= 1
        if i3 == 0:
            self.bprev.Disable()
        self.bnext.Enable()
        self.label.SetLabel('%d of %d' % (i3+1,n3))
        self.sketch.Redraw(i3)
    def OnSketchMotion(self,event):
        x,y = event.GetPositionTuple()
        position = self.sketch.Position(x,y)
        self.status.SetStatusText(position,0)
        event.Skip()
    def OnKeyDown(self,event):
        keycode = event.GetKeyCode()
        if chr(keycode) == 'q':
            self.Quit(event)
        event.Skip()
    def PickColor(self,event):
        dlg = wx.ColourDialog(self)
        dlg.GetColourData().SetChooseFull(True)
        if dlg.ShowModal() == wx.ID_OK:
            color = dlg.GetColourData().GetColour()
            self.sketch.SetColor(color)
        dlg.Destroy()
    def Quit(self,event):
        sys.exit(0)

@atexit.register
def cleanup():
    global ppms, picks
    for i in range(n3):
        for pick in picks[i].values():
            sys.stdout.write('%g\t%g\t%d\n' % pick[0])
    for ppm in ppms:
        if os.path.isfile(ppm):
            os.unlink(ppm)
    if os.path.isfile(byte):
        os.system(' '.join([sfrm,byte]))

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    frame.Show(True)
    app.MainLoop()


