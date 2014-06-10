#!/usr/bin/env python
'Simple interactive picking'

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
    from Tkinter import *
    from tkColorChooser import askcolor
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s file.rsf pick=pick.txt <sfgrey options>\n\n' % sys.argv[0])
    sys.exit(2)

inp = sys.argv[1]

bindir   = os.path.join(rsf.prog.RSFROOT,'bin')
sfgrey   = os.path.join(bindir,'sfgrey')
sfget    = os.path.join(bindir,'sfget')
ppmpen   = os.path.join(bindir,'ppmpen')
ppm = None

def rsf2image(rsf):
    global ppm
    ppm = os.path.splitext(rsf)[0]+'.ppm'
    command = '< %s %s %s | %s > %s' % (rsf,sfgrey,' '.join(sys.argv[2:]),ppmpen,ppm)
    if os.system(command) or not os.path.isfile(ppm):
        sys.stderr.write('Failed to execute "%s"\n\n' % command)
        sys.exit(3)
    img = PhotoImage(file=ppm)
    return img

color='#ffff00'
def pickcolor():
    global color
    col = askcolor()
    color = col[1]
   
root = Tk()
root.title(inp)

coords = StringVar()

frame = Frame(root)

button = Button(frame,text='Quit',bg='red',command=sys.exit)
button.pack(side=RIGHT)

button = Button(frame,text='Set Color',command=pickcolor)
button.pack(side=RIGHT)

label = Label(frame,textvariable=coords)
label.pack()

frame.pack(side=BOTTOM,fill=X)

width = 1024
height = 768

x0 = 170
x1 = 937

y0 = 96
y1 = 671

r = 5 # circle radius

picks =[]

canvas = Canvas(root,cursor='crosshair',
                width=width,height=height,
                background='black')

def hist(inp,func,var,default=None):
    command = '< %s %s %s parform=n' % (inp,sfget,var)
    devnull = open(os.devnull,"w")
    pipe = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=devnull,shell=True)
    val = pipe.stdout.read().rstrip()
    if val:
        return func(val)
    else:
        return default

label1 = hist(inp,str,'label1','Y')
unit1  = hist(inp,str,'unit1','')

label2 = hist(inp,str,'label2','X')
unit2  = hist(inp,str,'unit2','')

n1 = hist(inp,int,'n1',1)
d1 = hist(inp,float,'d1',1.0)
o1 = hist(inp,float,'o1',0.0)
yscale = (n1-1)*d1/(y1-y0)

n2 = hist(inp,int,'n2',1)
d2 = hist(inp,float,'d2',1.0)
o2 = hist(inp,float,'o2',0.0)
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

def getpick(event):
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)   
    if x >= x0 and y >= y0 and x <= x1 and y <= y1:
        canvas.create_oval(x-r,y-r,x+r,y+r,fill=color)
        xs = o2+(x-x0)*xscale
        ys = o1+(y-y0)*yscale
        picks.append((ys,xs))

image = rsf2image(inp)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.bind("<Motion>",display)
canvas.bind("<Button-1>",getpick)
canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    global ppm, picks
    for pick in picks:
        sys.stdout.write('%g\t%g\n' % pick)
    if os.path.isfile(ppm):
        os.unlink(ppm)

def bye(event):
    sys.exit(0)
 
root.bind("q",bye)
root.mainloop()

# Add:
# 2. remove picks with Button-3
# 3. 3-D
