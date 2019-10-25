#!/usr/bin/env python
##   Copyright (C) 2014 University of Texas at Austin
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

# modified from vp_annotate by Martin Karrenbach
import os, sys, tempfile, signal
from rsf.prog import RSFROOT
        
def handler(signum, frame):
    'signal handler for abortion [Ctrl-C]'
    sys.stderr.write('\n[Ctrl-C] Aborting...\n')
    if child:
        os.kill (signal.SIGINT,child)
    sys.exit(-1)

signal.signal(signal.SIGINT,handler) # handle interrupt

child = None
def syswait(comm):
    'Interruptable system command'
    global child
    child = os.fork()
    if child:
        (pid,exit) = os.waitpid(child,0)
        child = 0
        return exit
    else:
        os.system(comm)
        os._exit(0)

def annotate(files,args,interactive,textfile):
    inp = files[0]

    message = '''
Move cursor to the place, where the balloon arrow 
should point to, then click left mouse button.
Fill out the label and eventually change the defaults,
then click on CONFIRM.  Repeat for more annotations.
To create the annotated file, QUIT out of xtpen.
    '''

    xtpen = os.path.join(RSFROOT,'bin','xtpen')
    box   = os.path.join(RSFROOT,'bin','sfbox')
    vppen = os.path.join(RSFROOT,'bin','vppen')

    if interactive:
        run = '%s message="%s" %s interact=%s boxy=y < %s' % (xtpen,message,args,textfile,inp)
        syswait(run)

    boxes = []

    try:
        tfile = open(textfile,'r')
        for line in tfile.readlines():
            boxvpl = tempfile.mktemp(suffix='.vpl')
    
            run = '%s %s %s > %s' % (box,line.rstrip(),args,boxvpl)
            syswait(run)

            boxes.append(boxvpl)
        tfile.close()
    except:
        pass

    boxvpl = ' '.join(boxes)

    result = '''
This is the annotated vplot figure.
You might play with vpstyle=y, if you only want to 
see the original portion.
    '''
   
    if interactive:
        run = '%s %s %s erase=once vpstyle=n %s | %s  message="%s" ' % (vppen,inp,boxvpl,args,xtpen,result)
    else:
        run = '%s %s %s erase=once vpstyle=n %s > %s' % (vppen,inp,boxvpl,args,files[1])
    syswait(run)

    # cleanup
    for boxvpl in boxes:
        os.unlink(boxvpl)

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)
    prog = sys.argv.pop(0)
    
    usage = '''
    Annotates a Vplot file with a box.

    Usage:
    %s [batch=0] [text=box.par] < file.vpl [annotated.vpl]
    ''' % prog

    if argc < 2:
        print(usage)
        sys.exit(2)

    interactive = 1
    textfile = "box.par" 

    args = []
    files = []
    for arg in sys.argv:
        if '=' in arg:
            if arg[:6] == 'batch=':
                if arg[6]=='y' or arg[6]=='1':
                    interactive = 0
                else:
                    interactive = 1
            elif arg[:5] == 'text=':
                textfile = arg[5:]
            else:
                args.append(arg)
        else:
            files.append(arg)
    args = ' '.join(args)

    if interactive:
        needfiles = 1
    else:
        needfiles = 2
    
    if len(files) < needfiles:
        print(usage)
        sys.exit(3)

    annotate(files,args,interactive,textfile)
