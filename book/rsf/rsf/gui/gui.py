#!/usr/bin/env python
import os, sys

try:
    from Tkinter import *
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)

root = Tk()
root.title('Wavelet Demo')

wtype = StringVar()
wtype.set('b')

type_frame = Frame(root,relief=SUNKEN,borderwidth=2)
type_frame.pack(side=TOP,fill=X)

Label(type_frame,text='Wavelet Type').pack(side=TOP)
types = {'h':'Haar',
         'l':'Linear',
         'b':'Bi-orthogonal'}
for t in 'hlb':
    rbut = Radiobutton(type_frame,text=types[t],value=t,variable=wtype)
    rbut.pack(side=LEFT)

pclip_frame = Frame(root,relief=SUNKEN,borderwidth=2)
pclip_frame.pack(side=TOP,fill=X)
  
pclip = IntVar()
pclip.set(50)

scale = Scale(pclip_frame,from_=1,to=99,resolution=1,orient=HORIZONTAL,
              variable=pclip,length=200)
scale.pack(side=RIGHT)
Label(pclip_frame,text='Threshold\nPercentile').pack(side=RIGHT,anchor=SE)

frame = Frame(root)
frame.pack(side=TOP,fill=X)

quit = Button(frame,text='Quit',background='red',command=sys.exit)
quit.pack(side=RIGHT)

def scons():
    'Get parameters from GUI and pass them to SCons'
    os.system ('scons -Q type=%s pclip=%d view' % (wtype.get(),pclip.get()))
    
cycle = Button(frame,text='Run',background='yellow',command=scons)
cycle.pack()

root.mainloop()

