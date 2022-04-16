import numpy as np
from scipy.sparse import spdiags

def idnop1D(n):
    e = np.ones( n, dtype='float')    
    Lopr = spdiags([e], [0], n,n)    
    return Lopr

def idnop2D(nx,nz):
    n = nx*nz
    e = np.ones( n, dtype='float') 
    
    Lopr = spdiags([e], [0], n,n)
    return Lopr
