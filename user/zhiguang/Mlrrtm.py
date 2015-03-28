#!/usr/bin/env python
'Lowrank prestack RTM as a linear operator'

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
import os, sys, tempfile, subprocess 
import rsf.prog, rsf.path

# Madagascar bin directory
bindir = os.path.join(rsf.prog.RSFROOT,'bin')

# Madagascar DATAPATH
datapath = rsf.path.datapath().rstrip('/')

# Madagascar commands
cp = os.path.join(bindir,'sfcp')
rtm = os.path.join(bindir,'sfmpilrrtm')

# Random files for input and output
inpd,inpfile = tempfile.mkstemp(dir=datapath)
outd,outfile = tempfile.mkstemp(dir=datapath)

p = subprocess.Popen([cp],stdout=inpd,close_fds=True)
p.wait()

run = 'ibrun tacc_affinity %s input=%s output=%s %s' % (rtm,inpfile,outfile,' '.join(sys.argv[1:]))
print run

os.system(run)

p = subprocess.Popen([cp],stdin=outd)
p.wait()

