#!/usr/bin/env python
'Explore project results'

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
import sys

try:
    from Tkinter import *
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)
    
from signal import *
import time, os, commands

pid = 0

def handler(signum,frame):
    global pid
    sys.stdout.write("\n%s: aborting...\n" % sys.argv[0])
    if pid:
        os.kill(pid,SIGINT)
    sys.exit(1)

signal(SIGINT,handler)

root = Tk()

date = time.asctime()
root.title(date[:19])

frame = Frame(root,relief=GROOVE,borderwidth=2)
frame.pack(side=TOP,fill=X)

def showall():
    global results, pid
    pid = os.fork()
    if not pid:
        for fig in results:
            os.system ("scons %s.view" % fig)
        sys.exit()

quit = Button(frame,text="Quit",background="red",command=sys.exit)
quit.pack(side=LEFT)

cycle = Button(frame,text="Cycle",background="yellow",command=showall)
cycle.pack(side=RIGHT)

results = commands.getoutput("scons -s results").split()

def show(fig):
    def showfig():
        os.system("scons %s.view" % fig)
    return showfig

for fig in results:
    button = Button(root,text=fig,cursor='hand2',command=show(fig))
    button.pack(fill=X)

root.mainloop()
