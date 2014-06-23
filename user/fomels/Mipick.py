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
import sys, os, atexit, tempfile, subprocess
from rsf.prog import RSFROOT

try:
    from Tkinter import *
    from tkColorChooser import askcolor
except:
    sys.stderr.write('Please install Tkinter!\n\n')
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

def rsf2image(i3):
    global byte
    ppm = '%s%d.ppm' % (name,i3)
    if not os.path.isfile(ppm):
        command = '< %s %s n3=1 f3=%d | %s %s | %s %s > %s' % \
            (byte,sfwindow,i3,sfgrey,args,ppmpen,args,ppm)
        if os.system(command) or not os.path.isfile(ppm):
            sys.stderr.write('Failed to execute "%s"\n\n' % command)
            sys.exit(3)
        ppms.append(ppm)
    img = PhotoImage(file=ppm)
    return img

color='#ffff00'
def pickcolor():
    global color
    col = askcolor()
    color = col[1]

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

label1 = hist(str,'label1','Y')
unit1  = hist(str,'unit1','')

label2 = hist(str,'label2','X')
unit2  = hist(str,'unit2','')

n1 = hist(int,'n1',1)
d1 = hist(float,'d1',1.0)
o1 = hist(float,'o1',0.0)
yscale = (n1-1)*d1/(y1-y0)

n2 = hist(int,'n2',1)
d2 = hist(float,'d2',1.0)
o2 = hist(float,'o2',0.0)
xscale = (n2-1)*d2/(x1-x0)

n3 = hist(int,'n3',1)
i3 = 0

root = Tk()
root.title('Interactive Picking')

coords = StringVar()

frame = Frame(root)

button = Button(frame,text='Quit',bg='red',command=sys.exit)
button.pack(side=RIGHT)

button = Button(frame,text='Set Color',command=pickcolor)
button.pack(side=RIGHT)

if n3 > 1:
    def nextframe(event=None):
        global i3,n3,canvas,next,prev,picks
        for pick in picks[i3].keys():
            canvas.itemconfigure(pick,state=HIDDEN)
        i3 += 1
        if i3 == n3-1:
            next.config(state='disabled')
        prev.config(state='normal')
        for pick in picks[i3].keys():
            canvas.itemconfigure(pick,state=NORMAL)   
        i3frame.set('%d of %d' % (i3+1,n3))
        image = rsf2image(i3)
        canvas.itemconfigure('image', image=image)
        canvas.image=image

    def prevframe(event=None):
        global i3,n3,canvas,next,prev,picks
        for pick in picks[i3].keys():
            canvas.itemconfigure(pick,state=HIDDEN)
        i3 -= 1
        if i3 == 0:
            prev.config(state='disabled')
        next.config(state='normal')
        for pick in picks[i3].keys():
            canvas.itemconfigure(pick,state=NORMAL)   
        i3frame.set('%d of %d' % (i3+1,n3))
        image = rsf2image(i3)
        canvas.itemconfigure('image', image=image)
        canvas.image=image

    i3frame = StringVar()
    i3frame.set('%d of %d' % (i3+1,n3))

    next = Button(frame,text='Next >',command=nextframe)
    next.pack(side=LEFT)

    prev = Button(frame,text='< Prev',command=prevframe)
    prev.config(state='disabled')
    prev.pack(side=LEFT)

    label = Label(frame,textvariable=i3frame,relief=RIDGE,borderwidth=3)
    label.pack(side=LEFT)

    root.bind('n',nextframe)
    root.bind('p',prevframe)
    root.bind('m',prevframe)

label = Label(frame,textvariable=coords)
label.pack()

frame.pack(side=BOTTOM,fill=X)

picks = []
npick = 0

for i in range(n3):
    picks.append({})

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
    global o1,o2,x0,y0,xscale,yscale,i3
    xs = o2+(x-x0)*xscale
    ys = o1+(y-y0)*yscale
    return (ys,xs,i3)

def selectpick(event):
    global current
    current = canvas.find_closest(event.x, event.y)

def movepick(event):
    global current,r
    if current:
        x=event.x
        y=event.y
        if x < x0:
            x=x0
        elif x > x1:
            x=x1
        if y < y0:
            y=y0
        elif y > y1:
            y=y1
        canvas.coords(current,x-r,y-r,x+r,y+r)

def movedpick(event):
    global current
    tag = canvas.gettags(current)[0]
    x=event.x
    y=event.y
    if x < x0:
        x=x0
    elif x > x1:
        x=x1
    if y < y0:
        y=y0
    elif y > y1:
        y=y1
    picks[i3][tag] = scalepick(x,y)
    current = None

def deletepick(event):
    pick = canvas.find_closest(event.x, event.y)
    tag = canvas.gettags(pick)[0]
    del picks[i3][tag]
    canvas.delete(pick)

def addpick(event):
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
        picks[i3][tag]=scalepick(x,y)

image = rsf2image(0)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.image=image
canvas.bind("<Motion>",display)
canvas.bind("<Button-1>",addpick)
canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    global ppms, picks
    for i in range(n3):
        for pick in picks[i].values():
            sys.stdout.write('%g\t%g\t%d\n' % pick)
    for ppm in ppms:
        if os.path.isfile(ppm):
            os.unlink(ppm)
    if os.path.isfile(byte):
        os.system(' '.join([sfrm,byte]))

def bye(event):
    sys.exit(0)
 
root.bind("q",bye)
root.mainloop()

