#!/usr/bin/env python
'''Converts 2D velocity files from handvel.txt to handvel.rsf
-The program converts time samples from ms to s

-The rsf output file will have traces equal to the number
of CMP locations in handvel.txt. You need to interploate
between traces for a denser grid e.g. using sfremap1

-CMP locations in handvel.txt are not used in the program.
o2 tells the program where to put the first trace and d2
tells the program of the locations of the remaining traces.

-This program uses sfinvbin1 default parameters. The program
could possibly be enhanced to work with additional sfinvbin1
parameters.
'''

##   Copyright (C) 2007 University of Texas at Austin
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
import sys, os, string , tempfile
try: 
    import subprocess
except:
    sys.stderr.write("subprocess module is needed")
    sys.exit(2) 

basename=os.path.basename(sys.argv[0])
usage= '''
Name
        %s
DESCRIPTION 
        Converts 2D velocity files from handvel.txt to handvel.rsf
SYNOPSIS
        %s < handvels.txt o1=0 d1=.001 n1=3000 o2=5391.88 d2=625 > handvel.rsf
PARAMETERS
        float   o1= origin of the first axis
        float   d1= sampling in the first axis
        int     n1= size of the first axis
        float   o2= origin of the second axis
        float   d2= sampling in the second axis
             
COMMENTS:
        -The program converts time samples from ms to s

        -The rsf output file will have traces equal to the number
         of CMP locations in handvel.txt. You need to interploate
         between traces for a denser grid e.g. using sfremap1

        -CMP locations in handvel.txt are not used in the program.
         o2 tells the program where to put the first trace and d2
         tells the program of the locations of the remaining traces.

        -This program uses sfinvbin1 default parameters. The program
         could possibly be enhanced to work with additional sfinvbin1
         parameters.
SOURCE
        %s
                  
''' %(basename,basename,sys.argv[0])

bindir   = os.path.join(rsf.prog.RSFROOT,'bin')
sfcat    = os.path.join(bindir,'sfcat')
sfrm     = os.path.join(bindir,'sfrm')
sfinvbin1  = os.path.join(bindir,'sfinvbin1')
sfdd       = os.path.join(bindir,'sfdd')
sfput      = os.path.join(bindir,'sfput')
datapath = rsf.path.datapath().rstrip('/')

vs=[]
def myfunction(l,its):
    global sfcat, sfrm, sfinvbin1, sfdd, sfput, datapath
    tcmd='echo '
    vcmd='echo '
    
    # open temp files 
    td,tpath = tempfile.mkstemp(suffix=".rsf",dir=datapath,text=True)
    vd,vpath = tempfile.mkstemp(suffix=".rsf",dir=datapath,text=True)
    tb,tpathb = tempfile.mkstemp(suffix=".rsf",dir=datapath)
    vb,vpathb = tempfile.mkstemp(suffix=".rsf",dir=datapath)
    tvd,tvpath = tempfile.mkstemp(suffix=".rsf",dir=datapath)

    # let me put vel for time o1
    if its[0] != o1 :
        t=o1
        v=its[1]
        # insert v then t
        its.insert(0,v)
        its.insert(0,t)

    # create cmds for time and velocity files
    for i in range(0,len(its),2):
        # time samples are converted to seconds
        tcmd=tcmd+ str(float(its[i])/1000.0)+" "
        vcmd=vcmd+ str(float(its[i+1]))+" "
        
    tcmd= tcmd + '''n1=%d data_format=ascii_float in=%s\
                 '''%((len(its)/2),tpath)
    #print tcmd 
    vcmd= vcmd + '''n1=%d data_format=ascii_float in=%s\
                 '''%((len(its)/2),vpath)
    #print vcmd
    nx=n1
    dx=d1
    x0=o1
    tvcmd=''' %s head=%s nx=%d dx=%f x0=%f \
          '''%(sfinvbin1,tpathb,nx,dx,x0)
    #print tvpath
    #print vpathb
    #print tpathb    
    # execute cmds
    subprocess.call(tcmd,stdout=td,shell=True)
    subprocess.call(vcmd,stdout=vd,shell=True)
    os.lseek(td,0,0)
    os.lseek(vd,0,0)
    subprocess.call('%s form=native'%(sfdd),stdin=td,stdout=tb,shell=True)
    subprocess.call('%s form=native'%(sfdd),stdin=vd,stdout=vb,shell=True)
    os.lseek(tb,0,0)
    os.lseek(vb,0,0)  
    subprocess.call(tvcmd,stdin=vb,stdout=tvd,shell=True)

    # maintain a list of interpolated traces
    vs.append(tvpath)

    # close files
    for k in [td,vd,tb,vb,tvd]:
       os.close(k)
     
    # remove files  
    for k in [tpathb,vpathb]:
       try:
           subprocess.call(sfrm + ' ' + k,shell=True)
       except:
           pass
    for k in [tpath,vpath]:
       os.remove(k)
    
    return
if __name__ == "__main__":
    par=salah.Par()
    n1=par.int("n1")   # size of the first axis
    o1=par.float("o1") # origin of the first axis
    d1=par.float("d1") # sampleing in the first axis
    o2=par.float("o2") # origin of the second axis
    d2=par.float("d2") # sampling in the second axis
    if not (n1 or o1 or d1 or o2 or d2):
       #sys.stderr.write(usage)
       rsf.prog.selfdoc()
       sys.exit(2)
    if sys.stdout.isatty() or sys.stdin.isatty():
       #sys.stderr.write(usage)
       rsf.prog.selfdoc() 
       sys.exit(2)
    items=[]
    for line in sys.stdin:
        if line.startswith("HANDVEL"):
           loc=line.split()[1]
           if items:
              myfunction(loc,items)
              items=[]
        else:
           items = items + line.split()
    myfunction(loc,items)
    
    # concatinate traces in the second axis
    cmd='%s axis=2 %s | %s d2=%f o2=%f '%(sfcat,' '.join(vs),sfput,d2,o2,)
    #print cmd
    subprocess.call(cmd,stdout=sys.stdout,shell=True)

    # removing temp files of individual traces
    for tmp in vs:
        try:
           subprocess.call(sfrm + ' ' + tmp,shell=True)
        except:
           pass
    sys.exit(0)




