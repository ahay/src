#!/usr/bin/env python
import os, sys

try:
    from Tkinter import *
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)

root = Tk()
root.title('Canny Demo')

rect_frame = Frame(root,relief=SUNKEN,borderwidth=2)
rect_frame.pack(side=TOP,fill=X)

rect = IntVar()
rect.set(10)

scale = Scale(rect_frame,from_=1,to=50,resolution=1,orient=HORIZONTAL,
              variable=rect,length=200)
scale.pack(side=RIGHT)
Label(rect_frame,text='Smoothing Radius').pack(side=RIGHT,anchor=SE)

frame = Frame(root)
frame.pack(side=TOP,fill=X)

quit = Button(frame,text='Quit',background='red',command=sys.exit)
quit.pack(side=RIGHT)

def scons():
    'Get parameters from GUI and pass them to SCons'
    os.system ('scons -Q rect=%d overlay.view' % rect.get())

cycle = Button(frame,text='Run',background='yellow',command=scons)
cycle.pack()

root.mainloop()

