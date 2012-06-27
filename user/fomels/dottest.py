#!/usr/bin/env python
'Dot product test'

##   Copyright (C) 2008 University of Texas at Austin
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

import rsf.api as rsf

random = rsf.noise(rep=1,type=0)

def dottest(oper,mod,dat):
    'Dot product test'
    mod = random[mod]
    dat = random[dat]
    print " L[m]*d=%g" % oper(adj=0)[mod].dot(dat)
    print "L'[d]*m=%g" % oper(adj=1)[dat].dot(mod)

if __name__ == "__main__":
    # Create random matrix
    matrix = rsf.spike(n1=10,n2=5,d1=1,d2=1).noise(rep=1,type=0)[0]
    # Matrix multiplication operator
    oper = rsf.matmult(mat=matrix)
    # Model space and data space vectors
    model  = rsf.window(n2=1)[matrix]
    data   = rsf.window(n1=1)[matrix]
    # Using sfdottest
    rsf.dottest(0,oper,mod=model,dat=data)
    # Using function above
    dottest(oper,model,data)


