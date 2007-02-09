#!/usr/bin/env python
 
# Author: G. Hennenfent
#         Seismic Laboratory for Imaging and Modeling
#         Department of Earch & Ocean Sciences
#         The University of British Columbia
#         
# Date  : January, 07

#  Copyright (C) 2006 The University of British Columbia at Vancouver
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy as np
import rsf as sf
try:
    import pyct as ct
except:
    import sys
    sys.stderr.write('''
    sffdct needs PyCurvelab and did not find it on your system. Check
    your PYTHONPATH or go to
    https://wave.eos.ubc.ca/Software/Licenced/ to download it.\n\n''')
 
par = sf.Par()
input  = sf.Input()
output = sf.Output()
assert 'float' == input.type,"sffdct needs float input"

### read FDCT parameters from cmd line
# number of scale for the decomposition
nbs = par.int("nbs") 
# number of angle at the 2nd coarsest scale
nba = par.int("nba")
# curvelets at finest scale
ac = par.bool("ac") 
if ac is None:
    ac = False
# adj/inv transform
adj = par.bool("adj",False) 

if adj:
    ### inverse transform
    assert input.size(1) == 1,"sffdct with adj=y needs a vector input"
    n1,n2 = input.ints("sizes",2)
    assert n1 and n2,"vector provided is not a curvelet vector (missing sizes header)"

    # load input in memory
    x = np.zeros(input.int("n1"),'f')
    input.read(x)

    # apply transform
    shot = np.float32(ct.fdct2((n2,n1),nbs,nba,ac).inv(x) )

    # write result to file
    output.put("n1",n1)
    output.put("n2",n2 )
    output.write(shot)
else:
    ### forward transform
    # read in size of input
    n1 = input.int("n1")
    n2 = input.int("n2")
    ni = input.size(2)
    assert ni == 1,"sffdct needs 2D input"
 
    # load input in memory
    shot = np.zeros((n2,n1),'f')
    input.read(shot)

    # apply transform
    x = np.float32(ct.fdct2(shot.shape,nbs,nba,ac).fwd(shot) )

    # write result to file
    output.put("n1",len(x) )
    output.put("n2",1 )
    output.put("sizes",[n1,n2] )
    output.write(x)
