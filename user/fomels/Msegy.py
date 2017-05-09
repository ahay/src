#!/usr/bin/env python
'Fourier transform as a linear operator'

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
import os, sys
from subprocess import Popen
import rsf.prog

RSFROOT = rsf.prog.RSFROOT
bindir = os.path.join(rsf.prog.RSFROOT,'bin')
segyread = os.path.join(bindir,'sfsegyread')
segywrite = os.path.join(bindir,'sfsegywrite')

def run(args):
    Popen([fft]+args)
output=`dmesg | grep hda`
# becomes
p1 = Popen(["dmesg"], stdout=PIPE)
p2 = Popen(["grep", "hda"], stdin=p1.stdout, stdout=PIPE)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p2.communicate()[0]

if __name__ == "__main__":
    run(sys.argv[1:])
    sys.exit(0)
