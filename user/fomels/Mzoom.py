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
    from tkinter import *
except:
    sys.stderr.write('Please install tkinter!\n\n')
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

def hist(inp,func,var,default=None):
    command = '< %s %s %s parform=n' % (inp,sfget,var)
    devnull = open(os.devnull,"w")
    pipe = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=devnull,shell=True)
    val = pipe.stdout.read().rstrip().decode('utf-8')
    if val:
        return func(val)
    else:
        return default

label1 = hist(inp,str,'label1','Y')
unit1  = hist(inp,str,'unit1','')

label2 = hist(inp,str,'label2','X')
unit2  = hist(inp,str,'unit2','')

def setframe(inp):
    global o1,d1,o2,d2,xscale,yscale
    n1 = hist(inp,int,'n1',1)
    d1 = hist(inp,float,'d1',1.0)
    o1 = hist(inp,float,'o1',0.0)
    yscale = (n1-1)*d1/(y1-y0)

    n2 = hist(inp,int,'n2',1)
    d2 = hist(inp,float,'d2',1.0)
    o2 = hist(inp,float,'o2',0.0)
    xscale = (n2-1)*d2/(x1-x0)

setframe(inp)

def redraw(x1,y1,x2,y2):
    global o1,d1,o2,d2,xscale,yscale
    global canvas, image, rsffiles
    # convert to true coordinates
    x1 = o2+(x1-x0)*xscale
    x2 = o2+(x2-x0)*xscale
    xmin = min(x1,x2)
    xmax = max(x1,x2)
    y1 = o1+(y1-y0)*yscale
    y2 = o1+(y2-y0)*yscale
    ymin = min(y1,y2)
    ymax = max(y1,y2)
    inp2 = tempfile.mktemp(suffix='.rsf')
    command = '< %s %s min1=%g max1=%g min2=%g max2=%g > %s' % (inp,sfwindow,ymin,ymax,xmin,xmax,inp2)
    if os.system(command) or not os.path.isfile(inp2):
        sys.stderr.write('Failed to execute "%s"\n\n' % command)
        sys.exit(4)
    rsffiles.append(inp2)
    setframe(inp2)
    image = rsf2image(inp2)
    canvas.itemconfigure('image', image=image)
    canvas.delete('window')

def restore(event):
    global inp, canvas, image
    setframe(inp)
    image = rsf2image(inp)
    canvas.itemconfigure('image', image=image)

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
    global xstart, ystart
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)
    xend = min(max(x,x0),x1)
    yend = min(max(y,y0),y1)
    redraw(xstart,ystart,xend,yend)


image = rsf2image(inp)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.bind("<Motion>",display)
canvas.bind("<ButtonPress-1>",startwin)
canvas.bind("<B1-Motion>",drawwin)
canvas.bind("<ButtonRelease-1>",endwin)
canvas.bind("<Button-3>",restore)
canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    global ppmfiles,rsffiles,sfrm
    for ppm in filter(os.path.isfile,ppmfiles):
        os.unlink(ppm)
    for rsf in filter(os.path.isfile,rsffiles):
        os.system(' '.join([sfrm,rsf]))

def bye(event):
    sys.exit(0)

root.bind("q",bye)
root.bind("<Escape>",restore)
root.mainloop()

