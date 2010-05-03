#!/usr/bin/env python
'''
Return mask to remove random traces in 2D and 3D using a maximum gap
size constraint
'''

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

import rsf.api as sf
import numpy as np

try: # Give precedence to local version
    from hegilles import killtraces
except: # Use distributed version
    from rsf.user.hegilles import killtraces

par = sf.Par()

input  = sf.Input()
output = sf.Output()

n1 = input.int("n1")
n2 = input.int("n2")
ni = input.size(2)
if ni != 1:
    n3 = ni
    ni = input.size(3)
else:
    n3 = 1
assert ni == 1,"sfkilltraces needs 2D or 3D input"

perc = par.float("perc",.75) # percentage of traces to remove
assert (perc>0 and perc<1),"perc should be between 0 and 1"

maxfactor = par.float("maxfactor",1.) # maximum gap factor
assert maxfactor>=1,"maxfactor should be greater or equal to 1"

seed = par.int("seed",np.random.randn() ) # seed for random number generator

output.put("n1",n2)
output.put("n2",n3)
output.put("n3",1)

for i in range(n3):
    mask = np.float32(killtraces(n2,100*(1-perc),maxfactor,seed))
    output.write(mask)
    seed += 1

# $Id$
