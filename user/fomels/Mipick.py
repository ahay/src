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

def hist(inp,func,var,default=None):
    command = '< %s %s %s parform=n' % (inp,sfget,var)
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

n3 = hist(inp,int,'n3',1)
d3 = hist(inp,float,'d3',1.0)
o3 = hist(inp,float,'o3',0.0)

i3 = 0

root = Tk()
root.title(inp)

coords = StringVar()

frame = Frame(root)

button = Button(frame,text='Quit',bg='red',command=sys.exit)
button.pack(side=RIGHT)

button = Button(frame,text='Set Color',command=pickcolor)
button.pack(side=RIGHT)

def nextframe():
    global i3
    i3 += 1
    if i3 == n3-1:
        next.config(state='disabled')
    prev.config(state='normal')
    i3frame.set('%d of %d' % (i3+1,n3))

def prevframe():
    global i3
    i3 -= 1
    if i3 == 0:
        prev.config(state='disabled')
    next.config(state='normal')
    i3frame.set('%d of %d' % (i3+1,n3))

if n3 > 1:
    i3frame = StringVar()
    i3frame.set('%d of %d' % (i3+1,n3))

    next = Button(frame,text='Next >',command=nextframe)
    next.pack(side=LEFT)

    prev = Button(frame,text='< Prev',command=prevframe)
    prev.config(state='disabled')
    prev.pack(side=LEFT)

    label = Label(frame,textvariable=i3frame,relief=RIDGE,borderwidth=3)
    label.pack(side=LEFT)

label = Label(frame,textvariable=coords)
label.pack()

frame.pack(side=BOTTOM,fill=X)

picks = {}
npick = 0

canvas = Canvas(root,cursor='crosshair',
                width=width,height=height,
                background='black')

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

current = None

def scalepick(x,y):
    global o1,o2,x0,y0,xscale,yscale,o3,i3,d3
    xs = o2+(x-x0)*xscale
    ys = o1+(y-y0)*yscale
    if n3 > 1:
        zs = o3+i3*d3
        return (ys,xs,zs)
    else:
        return (ys,xs)

def selectpick(event):
    global current
    current = canvas.find_closest(event.x, event.y)

def movepick(event):
    global current,r
    if current:
        x=event.x
        y=event.y
        canvas.coords(current,x-r,y-r,x+r,y+r)

def movedpick(event):
    global current
    tag = canvas.gettags(current)[0]
    picks[tag] = scalepick(event.x,event.y)
    current = None

def deletepick(event):
    pick = canvas.find_closest(event.x, event.y)
    tag = canvas.gettags(pick)[0]
    canvas.delete(pick)
    del picks[tag]

def getpick(event):
    global npick,r
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)   
    if x >= x0 and y >= y0 and x <= x1 and y <= y1:
        npick += 1
        tag = 'pick%d' % npick
        canvas.create_oval(x-r,y-r,x+r,y+r,fill=color,tags=tag)
        canvas.tag_bind(tag,'<ButtonPress-2>',selectpick)
        canvas.tag_bind(tag,'<B2-Motion>',movepick)
        canvas.tag_bind(tag,'<ButtonRelease-2>',movedpick)
        canvas.tag_bind(tag,'<Button-3>',deletepick)
        picks[tag]=scalepick(x,y)

image = rsf2image(inp)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.bind("<Motion>",display)
canvas.bind("<Button-1>",getpick)
canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    global ppm, picks
    for pick in picks.values():
        if n3 > 1:
            sys.stdout.write('%g\t%g\t%g\n' % pick)
        else:
            sys.stdout.write('%g\t%g\n' % pick)
    if os.path.isfile(ppm):
        os.unlink(ppm)

def bye(event):
    sys.exit(0)
 
root.bind("q",bye)
root.mainloop()

# 3-D: load image, dots
