#!/usr/bin/env python
import os, sys

try:
    from Tkinter import *
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)

root = Tk()

frame = Frame(root,relief=GROOVE,borderwidth=2)
frame.pack(side=TOP,fill=X)

quit = Button(frame,text="Quit",background="red",command=sys.exit)
quit.pack(side=LEFT)

def scons():
    os.system ("scons -Q view")

cycle = Button(frame,text="Run",background="yellow",command=scons)
cycle.pack(side=RIGHT)

root.mainloop()

