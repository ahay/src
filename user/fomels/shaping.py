#!/usr/bin/env python
'Conjugate-gradient method'

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

def conjgrad(oper,shape,eps,d,p0,niter):
    'Conjugate-gradient algorithm for shaping regularization'
    p = p0
    x = shape(adj=0)[p]
    r = oper(adj=0)[x]-d      
    for iter in range(niter):
        gx = oper(adj=1)[r]-x*eps
        gp = shape(adj=1)[gx]+p*eps
        gx = shape(adj=0)[gp]
        gr = oper(adj=0)[gx]

        gn = gp.dot(gp)
        print "iter %d: %g" % (iter+1,gn)
        if 0==iter:
            sp = gp
            sx = gx
            sr = gr
        else:
            beta = gn/gnp
            sp = gp + sp*beta
            sx = gx + sx*beta
            sr = gr + sr*beta
        gnp = gn
        
        alpha = sr.dot(sr)+eps*(sp.dot(sp)-sx.dot(sx))
        alpha = -gn/alpha
        p = p + sp*alpha
        x = x + sx*alpha
        r = r + sr*alpha
    return x

if __name__ == "__main__":
    # test matrix and data
    matrix = rsf.File([[1,1,1,0],
                       [1,2,0,0],
                       [1,3,1,0],
                       [1,4,0,1],
                       [1,5,1,1]])
    y = rsf.File([3,3,5,7,9])
    x0 = rsf.File([0,0,0,0])
    # matrix multiplication operator
    matmult = rsf.matmult(mat=matrix)
    copy = rsf.cp(x=1)

    x = conjgrad(matmult,copy,1,y,x0,6)
    y2 = matmult[x]
    print x[:]
    print y2[:]

