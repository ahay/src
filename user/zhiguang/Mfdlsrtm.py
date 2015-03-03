#!/usr/bin/python
'Finite difference RTM as a linear operator'

import os, sys, tempfile, subprocess
import rsf.prog, rsf.path

# Madagascar bin directory
bindir=os.path.join(rsf.prog.RSFROOT,'bin')

# Madagascar DATAPATH
datapath=rsf.path.datapath().rstrip('/')

# Madagascar commands
cp=os.path.join(bindir,'sfcp')
rtm=os.path.join(bindir,'sfmpifdlsrtm')

# Random files for input and output
inpd,inpfile=tempfile.mkstemp(dir=datapath)
outd,outfile=tempfile.mkstemp(dir=datapath)

p=subprocess.Popen([cp],stdout=inpd, close_fds=True)
p.wait()

run='ibrun tacc_affinity %s input=%s output=%s %s' %(rtm, inpfile, outfile,' '.join(sys.argv[1:]))
print run

os.system(run)

p=subprocess.Popen([cp],stdin=outd)
p.wait
