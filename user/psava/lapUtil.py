import numpy as np
from scipy.sparse import spdiags

def lapop1D(n):
    e = np.ones( n, dtype='float')    
    Lopr = spdiags([e, -2*e, e], [-1, 0, +1], n,n)    
    return Lopr

def lapop2D(nx,nz):
    n = nx*nz
    e = np.ones( n, dtype='float') 
    
    u = e.copy()
    u[::nz]=0
    
    l = e.copy()
    l[nz-1::nz]=0
    
    Lopr = spdiags([e, l, -4*e, u,e], [-nz, -1, 0, +1, +nz], n,n)
    return Lopr
