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
gif = os.path.splitext(inp)[0]+'.gif'

bindir = os.path.join(rsf.prog.RSFROOT,'bin')
sfgrey = os.path.join(bindir,'sfgrey')
sfget  = os.path.join(bindir,'sfget')
gdpen  = os.path.join(bindir,'gdpen')

command = '< %s %s %s | %s bgcolor=b type=gif > %s' % \
          (inp,sfgrey,' '.join(sys.argv[2:]),gdpen,gif)
if os.system(command) or not os.path.isfile(gif):
    sys.stderr.write('Failed to execute "%s"\n\n' % command)
    sys.exit(3)

root = Tk()
root.title(inp)

coords = StringVar()

label = Label(root,textvariable=coords,relief=RIDGE,borderwidth=3)
label.pack(side=BOTTOM,fill=X)

width = 1024
height = 768

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

n2 = hist(int,'n2',1)
d2 = hist(float,'d2',1.0)
o2 = hist(float,'o2',0.0)
label2 = hist(str,'label2','X')
unit2  = hist(str,'unit2','')

def display(event):
    canvas = event.widget
    x = canvas.canvasx(event.x)
    y = canvas.canvasx(event.y)    
    if x >= 170 and y >= 96 and x <= 937 and y <= 671:
        x = o2+(x-170)*(n2-1)*d2/(937-170)
        y = o1+(y-96 )*(n1-1)*d1/(671-96 )
	coords.set("(%s = %g %s, %s = %g %s)" % (label1,y,unit1,
                                                 label2,x,unit2))
    else:
	coords.set("")

image = PhotoImage(file=gif)
canvas.create_image(0,0,image=image,anchor=NW,tags="image")
canvas.bind("<Motion>",display)

canvas.pack(side=BOTTOM)

@atexit.register
def cleanup():
    if os.path.isfile(gif):
        os.unlink(gif)

root.mainloop()

