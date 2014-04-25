#!/usr/bin/env python
'''Converts 2D/3D velocity files from handvel.txt to handvel.rsf

- sfhandvel2rsf < handvels.txt o1=0 d1=.001 n1=3000 > handvel.rsf

- The program converts time samples from ms to s

- The rsf output file will have traces equal to the number
  of CMP locations in handvel.txt. You need to interploate
  between traces for a denser grid e.g. using sfremap1

- This program uses sfspline for interpolation. 

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
import sys, os, string , tempfile, subprocess, collections, multiprocessing
import rsf.path, rsf.prog
import re
try: 
    import subprocess
except:
    sys.stderr.write("subprocess module is needed")
    sys.exit(2) 
import rsf.api as salah

basename=os.path.basename(sys.argv[0])
usage= '''
Name
        %s
DESCRIPTION 
        Converts 2D/3D velocity files from handvel.txt to handvel.rsf
SYNOPSIS
        %s < handvels.txt o1=0 d1=.001 n1=3000 > handvel.rsf
PARAMETERS
        float   o1= origin of the first axis
        float   d1= sampling in the first axis
        int     n1= size of the first axis
             
COMMENTS:
        -The program converts time samples from ms to s

        -The rsf output file will have traces equal to the number
         of CMP locations in handvel.txt. You need to interploate
         between traces for a denser grid e.g. using sfremap1

        -This program uses sfspline default parameters. The program
         could possibly be enhanced to work with additional sfinvbin1
         parameters.
         
        -Input file is not QCed, thus, check your input file
        
SOURCE
        %s
                  
''' %(basename,basename,sys.argv[0])

bindir   = os.path.join(rsf.prog.RSFROOT,'bin')
sfcat    = os.path.join(bindir,'sfcat')
sfrm     = os.path.join(bindir,'sfrm')
sfinvbin1  = os.path.join(bindir,'sfinvbin1')
sfdd       = os.path.join(bindir,'sfdd')
sfput      = os.path.join(bindir,'sfput')
sfspline = os.path.join(bindir,'sfspline')
datapath = rsf.path.datapath().rstrip('/')
sftransp = os.path.join(bindir,'sftransp')


def myfunction(l,vs,its,nx,x0,dx):
    global sfcat, sfrm, sfinvbin1, sfdd, sfput, datapath,sfspline
    vcmd='echo '

    # open temp files 
    vd,vpath = tempfile.mkstemp(suffix=".rsf",dir=datapath,text=True)
    tvd,tvpath = tempfile.mkstemp(suffix=".rsf",dir=datapath)
    
    # let me put vel for time o1
    if float(its[0]) != o1:
        t=o1
        v=its[1]
        # insert v then t
        its.insert(0,v)
        its.insert(0,t)

    # create cmds for velocity file
    for i in range(0,len(its),2):
        # time samples are converted to seconds
        vcmd=vcmd+ str(float(its[i])/1000.0)+" "+str(float(its[i+1]))+" "
        
    vcmd= vcmd + '''n1=2 n2=%d data_format=ascii_float in=%s\
                 '''%((len(its)/2),vpath)
    #print vcmd
    #nx=n1
    #dx=d1
    #x0=o1
    
    #time velocity command
    tvcmd='''%s form=native | %s n1=%d o1=%f d1=%f fp=0,0\
          '''%(sfdd,sfspline,nx,x0,dx)
    
    # execute cmds
    subprocess.call(vcmd,stdout=vd,shell=True)
    os.lseek(vd,0,0)
    
    #print tvcmd  
    subprocess.call(tvcmd,stdin=vd,stdout=tvd,shell=True)

    # maintain a list of interpolated traces
    
    #print vs
    # close files
    for k in [vd,tvd]:
       os.close(k)
     
    for k in [vpath]:
       os.remove(k)
    vs.append(tvpath)
    return

if __name__ == "__main__":
    mgr = multiprocessing.Manager()
    vs = mgr.list()
    par=salah.Par()
    n1=par.int("n1")   # size of the first axis
    o1=par.float("o1") # origin of the first axis
    d1=par.float("d1") # sampling in the first axis
    #
    if not (n1 or o1 or d1):
       #sys.stderr.write(usage)
       rsf.prog.selfdoc()
       sys.exit(2)
    if sys.stdout.isatty() or sys.stdin.isatty():
       #sys.stderr.write(usage)
       rsf.prog.selfdoc() 
       sys.exit(2)
    inline=collections.OrderedDict()
    loc=collections.OrderedDict()
    i=None
    x=None
    for line in sys.stdin:
        line=line.strip()
        if re.match(r"^\*", line):
           continue 
        if re.match(r"^\HANDVEL|^\VFUNC", line):
           if 3 != len(line.split()):
           	 sys.stderr.write("wrong input file format\n")
           	 sys.stderr.write("%s\n"%(line))
                 sys.exit(2)
           # get inline and xline
           i=line.split()[1]
           x=line.split()[2]
           if i in inline.keys():
           	inline[i].append(x)
           else:
           	inline[i]=[x]
	   if (i,x) in loc.keys():
		 sys.stderr.write("duplicate location %s,%s\n"%(i,x))
                 sys.exit(2)
           else:
		 loc[i,x]=[]
        else:
           loc[i,x]= loc[i,x] + line.split()
    
    #for y in inline.keys():
    	#for x,v in inline[y].items():
            #print y
            #print inline[y]
            
     # compute o2, d2, n2, o3, d3, and n3
    n2=len(inline.keys())
    d2=1. if n2==1 else float(inline.keys()[1])-float(inline.keys()[0])
    o2=inline.keys()[0]
    
    n3=len(inline[o2])
    d3=1. if n3==1 else float(inline[o2][1])-float(inline[o2][0])
    o3=inline[o2][0]
    lock=multiprocessing.Lock()
    #print "o2="+str(o2)+" d2="+str(d2)+" n2="+str(n2)+" o3="+str(o3)+" d3="+str(d3)+" n3="+str(n3)
    jobs=[]
    for y in loc.keys():
        #myfunction(loc[y],n1,o1,d1)
        p = multiprocessing.Process(target=myfunction, args=(lock,vs,loc[y],n1,o1,d1))
        jobs.append(p)
        p.start()
    #print len(jobs)
    for job in jobs:
	job.join()
    # concatinate traces in the second axis
    cmd='''
        %s axis=2 %s | %s  n3=%d o3=%f d3=%f n2=%d o2=%f d2=%f label1=time label2=xline label3=inline| %s plane=23
        '''%(sfcat,' '.join(vs),sfput,n2,float(o2),d2,n3,float(o3),d3,sftransp)
    
    #cmd='%s axis=2 %s'%(sfcat,' '.join(vs))
    #print cmd
    subprocess.call(cmd,stdout=sys.stdout,shell=True)

    # removing temp files of individual traces
    for tmp in vs:
        try:
           subprocess.call(sfrm + ' ' + tmp,shell=True)
        except:
           pass
    sys.exit(0)




