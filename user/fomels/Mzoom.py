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
import sys, os, subprocess, atexit
import rsf.prog

try:
    from Tkinter import *
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s file.rsf <sfgrey options>\n\n' % sys.argv[0])
    sys.exit(2)

inp = sys.argv[1]
ppm = os.path.splitext(inp)[0]+'.ppm'

bindir = os.path.join(rsf.prog.RSFROOT,'bin')
sfgrey = os.path.join(bindir,'sfgrey')
sfget  = os.path.join(bindir,'sfget')
ppmpen  = os.path.join(bindir,'ppmpen')

command = '< %s %s %s | %s > %s' % (inp,sfgrey,' '.join(sys.argv[2:]),ppmpen,ppm)
if os.system(command) or not os.path.isfile(ppm):
    sys.stderr.write('Failed to execute "%s"\n\n' % command)
    sys.exit(3)

root = Tk()
root.title(inp)

coords = StringVar()

label = Label(root,textvariable=coords,relief=RIDGE,borderwidth=3)
label.pack(side=BOTTOM,fill=X)

width = 1024
height = 768

x0 = 170
x1 = 937

y0 = 96
y1 = 671

canvas = Canvas(root,cursor='crosshair',
                width=width,height=height,
                background='black')

def hist(func,var,default=None):
    command = '< %s %s %s parform=n' % (inp,sfget,var)
    pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
    val = pipe.stdout.read().rstrip()
    if val:
        return func(val)
    else:
        return default

n1 = hist(int,'n1',1)
d1 = hist(float,'d1',1.0)
o1 = hist(float,'o1',0.0)
label1 = hist(str,'label1','Y')
unit1  = hist(str,'unit1','')
yscale = (n1-1)*d1/(y1-y0)

n2 = hist(int,'n2',1)
d2 = hist(float,'d2',1.0)
o2 = hist(float,'o2',0.0)
label2 = hist(str,'label2','X')
unit2  = hist(str,'unit2','')
xscale = (n2-1)*d2/(x1-x0)

def display(event):
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)    
    if x >= x0 and y >= y0 and x <= x1 and y <= y1:
        x = o2+(x-x0)*xscale
        y = o1+(y-y0)*yscale
	coords.set("(%s = %g %s, %s = %g %s)" % (label1,y,unit1,
                                                 label2,x,unit2))
    else:
	coords.set("")

def startwin(event):
    global xstart, ystart
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)
    xstart = min(max(x,x0),x1)
    ystart = min(max(y,y0),y1)
    canvas.delete('window')
    canvas.create_rectangle(xstart,ystart,xstart,ystart,
                            outline='yellow',tags='window')

def drawwin(event):
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)
    x = min(max(x,x0),x1)
    y = min(max(y,y0),y1)
    canvas.coords('window',xstart,ystart,x,y)

def endwin(event):
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)
    xend = min(max(x,x0),x1)
    yend = min(max(y,y0),y1)
    print "end with (xend,yend)"

image = PhotoImage(file=ppm)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.bind("<Motion>",display)
canvas.bind("<ButtonPress-1>",startwin)
canvas.bind("<B1-Motion>",drawwin)
canvas.bind("<ButtonRelease-1>",endwin)
canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    if os.path.isfile(ppm):
        os.unlink(ppm)

def bye(event):
    sys.exit(0)

root.bind("q",bye)
root.mainloop()

