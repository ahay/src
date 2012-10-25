#!/usr/bin/env python
'Conjugate-directions method'

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

def conjdir(oper,dat,x0,niter):
    'Conjugate-directions algorithm for minimizing |A x - dat|^2'
    x = x0
    R = oper(adj=0)[x]-dat      
    s,S,Sn=[],[],[]
    for iter in range(niter):
        f = random[x0]
        F = oper(adj=0)[f]
        for k in range(iter):
            beta = - F.dot(S[k])/Sn[k]
            f = f+s[k]*beta
            F = F+S[k]*beta
        Fn = F.dot(F)

        s.append(f)
        S.append(F)
        Sn.append(Fn)

        print "iter %d: %g" % (iter+1,Fn)
 
        alpha = - R.dot(F)/Fn
        x = x+f*alpha
        R = R+F*alpha
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

    # Using function above
    x = conjdir(matmult,y,x0,6)
    y2 = matmult[x]
    print x[:]
    print y2[:]

