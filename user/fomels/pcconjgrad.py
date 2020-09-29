#!/usr/bin/env python
'Complex Conjugate-gradient method with Preconditioner'

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

import m8r

def pcconjgrad(oper,prec,dat,p0,niter):
    'Conjugate-gradients for minimizing |oper prec p - dat|^2'
    p = p0
    x = prec(adj=0)[p0]
    R = oper(adj=0)[x]-dat   
    for iter in range(niter):
        f = oper(adj=1)[R]
        g = prec(adj=1)[f]
        F = prec(adj=0)[g]
        G = oper(adj=0)[F]
        gn = g.cdot2()
        print "iter %d: %g" % (iter+1,gn)
        if 0==iter:
            s = g
            S = G
        else:
            beta = gn/gnp
            s = g+s*beta
            S = G+S*beta
        gnp = gn
        alpha = -gn/S.cdot2()
        p = p+s*alpha
        R = R+S*alpha
    x = prec(adj=0)[p]
    return x

if __name__ == "__main__":
    # test matrix and data
    matrix = m8r.File([[1,1,1,0],
                       [1,2,0,0],
                       [1,3,1,0],
                       [1,4,0,1],
                       [1,5,1,1]])
    cmatrix = m8r.cmplx[matrix,matrix]
    cmatmult = m8r.cmatmult(mat=cmatrix)
    yr = m8r.File([3,2,5,6,9])
    yi = m8r.File([3,4,5,8,9])
    y = m8r.cmplx[yr,yi]
    x0 = m8r.File([0,0,0,0])
    x1 = m8r.cmplx[x0,x0]

    # Using sfcconjgrad
    x = m8r.cconjgrad(mod=x1,niter=6)[y,cmatmult]
    y1 = cmatmult[x]
    print x[:]
    print y1[:]

    # Using function above
    ccausint = m8r.ccausint
    x = pcconjgrad(cmatmult,ccausint,y,x1,6)
    y1 = cmatmult[x]
    print x[:]
    print y1[:]
    

