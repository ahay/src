from __future__ import print_function

import m8r
import numpy as np
import math as m
import sys

try:
    import numba
except:
    sys.stderr.write('Please install numba!\n\n')
    sys.exit(1)

# gradient functions

@numba.jit(parallel=True,nopython=True)
def grad3(b,dx=1.0):
    '''compute 3 dimensional derivative'''
    c = np.zeros(b.shape)
    c[:,:,1:-1] = (b[:,:,2:] - b[:,:,:-2])/2
    c[:,:,0] = b[:,:,1] - b[:,:,0]
    c[:,:,-1] = b[:,:,-1] - b[:,:,-2]
    return c/dx

@numba.jit(parallel=True,nopython=True)
def grad2(b,dx=1.0):
    '''compute 2 dimensional derivative'''
    c = np.zeros(b.shape)
    c[:,1:-1] = (b[:,2:] - b[:,:-2])/2
    c[:,0] = b[:,1] - b[:,0]
    c[:,-1] = b[:,-1] - b[:,-2]
    return c/dx

@numba.jit(nopython=True)
def gradn(a,axis,dx=1):
    '''compute 2 or 3 dimensional derivative'''
    if a.ndim == 3:
        if axis == 2:
            return grad3(a,dx)
        if axis == 1:
            transp_tuple = (0,2,1)
        if axis == 0:
            transp_tuple = (2,1,0)
        return np.transpose(grad3(np.transpose(a,transp_tuple),dx),transp_tuple)
    if a.ndim == 2:
        if axis == 1:
            return grad2(a,dx)
        if axis == 0:
            return np.transpose(grad2(np.transpose(a),dx))
    
    return
    
def numbgrad(a):
    '''compute gradient'''
    g = []
    for axis in range(a.ndim):
        g.append(gradn(a,axis))
    return g

# picking functions
@numba.jit(parallel=True,nopython=True)
def index_remainder(X, xo, dx):
    '''return the coodinate index and remainder'''
    inds = ((X-xo)//dx)
    rems = (X-(inds*dx + xo))/dx
    return[ inds,rems]

 
@numba.jit(parallel=True,nopython=True)
def numb_pick_element(ind,array):
    '''pick elements from an array over a surface'''
    arr2 = np.transpose(array,(1,0,2))
    out = np.zeros(ind.shape)
    for i in range(ind.shape[0]):
        for j in range(ind.shape[1]):
            out[i,j] = arr2[ind[i,j],i,j]
    return out

@numba.jit(parallel=True, nopython=True)
def numbmax(A,b):
    '''clip A so it is \geq b'''
    B = np.zeros(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i,j] > b:
                B[i,j] = A[i,j]
            else:
                B[i,j] = b
    return B

@numba.jit(parallel=True,nopython=True)
def numbmin(A,b):
    '''clip A lo it is \leq b'''
    B = np.zeros(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i,j] < b:
                B[i,j] = A[i,j]
            else:
                B[i,j] = b
    return B

@numba.jit(parallel=True,nopython=True)
def numb_interp(A,B,r):
    '''returns C = (1-r)*A + r*B'''
    C = np.zeros(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            C[i,j] = (1.0-r[i,j])*A[i,j]+r[i,j]*B[i,j]
    return C

@numba.jit(parallel=True,forceobj=True)
def numb_pick_surface(v, vo, dv, semb):
    '''pick a semblance surface using a velocity field'''
    nv = semb.shape[1]
    # get index and remainder
    inds, rems = index_remainder(v, vo, dv)
    inds0 = to_int(inds)
    inds1 = to_int(numbmin(inds0+1,nv-1))
    return numb_interp(numb_pick_element(inds0,semb),numb_pick_element(inds1,semb),rems)

# helpful numba functions 

@numba.jit(parallel=True,nopython=True)
def numba_multiply2(A,B):
    '''returns A*B'''
    return A*B

@numba.jit(parallel=True,nopython=True)
def numba_multiply3(A,B,C):
    '''returns A*B*C'''
    return A*B*C

@numba.jit(parallel=True,nopython=True)
def numb_add2(A,B):
    '''returns A+B'''
    return A+B

@numba.jit(parallel=True,nopython=True)
def numb_add3(A,B,C):
    '''returns A+B+C'''
    return A+B+C

@numba.jit(forceobj=True)
def numb_integrate(A):
    '''integrate over an array'''
    return np.sum(A)

@numba.jit(forceobj=True)
def numb_dot(A,B):
    '''compute the inner product of (A,B)'''
    return np.dot(A.flatten(),B.flatten())

@numba.jit(parallel=True,nopython=True)
def numb_divide(a,b):
    '''returns a/b pointwise'''
    return a/b

@numba.jit(nopython=True,parallel=True)
def numb_force_range(A,lo,hi):
    '''force array A to be between lo and hi'''
    C = A.flatten()
    B = np.zeros(C.size)
    for dex in range(C.size):
        B[dex]= max( lo, min( hi, C[dex]))
        
    return B.reshape(A.shape)

@numba.jit(parallel=True,forceobj=True)
def to_int(X):
    '''convert array X to int array'''
    return X.astype(int)


# helpful functions for computing cost, gradient
@numba.jit(parallel=True,nopython=True)
def lmbda_term(dv0, dv1, lmbda):
    '''create lambda term with added epsilon for stability'''
    lt = np.zeros(dv0.shape)
    for i in range(dv0.shape[0]):
        for j in range(dv0.shape[1]): 
            lt[i,j] = np.sqrt(dv0[i,j]**2+dv1[i,j]**2+lmbda**2)
    return lt

@numba.jit(parallel=True,nopython=True)
def epsilon_term(dv0, dv1, epsilon):
    '''create epsilon term for H1 regularization'''
    eps_magnitude = np.zeros(dv0.shape)
    if epsilon <= 0.0:
        return eps_magnitude
    for i in range(dv0.shape[0]):
        for j in range(dv0.shape[1]):
            eps_magnitude[i,j] = epsilon/2 * (dv0[i,j]**2 + dv1[i,j]**2)
    return eps_magnitude    

@numba.jit(parallel=True,nopython=True)
def numb_secondterm_interior(expsemb_v, lmbdatrm, dvi, epsilon):
    '''create interior of second'''
    return expsemb_v*(dvi/lmbdatrm+dvi*epsilon)

@numba.jit(parallel=True,forceobj=True)
def numb_secondterm(expsemb_v,lmbdatrm,dvi, epsilon, axis):
    '''create second term'''
    return gradn(numb_secondterm_interior(expsemb_v,lmbdatrm, dvi, epsilon),axis)

@numba.jit(parallel=True,forceobj=True)
def compute_cost(v, semb, dX, oX, lmbda, epsilon):
    '''compute cost for semb and v'''
    # get expsemb(v)
    expsemb_v = np.exp(-1*numb_pick_surface(v, oX[1], dX[1], semb))
    # get grad v
    dv = numbgrad(v)
    # get lambda term
    lmbdatrm = lmbda_term(dv[0]/dX[0], numb_divide(dv[1],dX[2]*v), lmbda)
    # and epsilon term 
    epsterm = epsilon_term(dv[0]/dX[0], numb_divide(dv[1],dX[2]*v), epsilon)
    # return result
    return numb_integrate(numba_multiply2(expsemb_v, numb_add2(lmbdatrm, epsterm)))*dX[0]*dX[2]


@numba.jit(forceobj=True)
def compute_grad(v, semb, dsembv, dX, oX, lmbda, epsilon):
    '''compute variational gradient'''
    # get expsemb(v)
    expsemb_v = np.exp(-1*numb_pick_surface(v, oX[1], dX[1], semb))
    # get dexpsemb/dv(v)
    dexpsemb_dv = numb_pick_surface(v, oX[1], dX[1], dsembv)
    # get gradv
    dv = numbgrad(v)
    # get lambda term
    lmbdatrm = lmbda_term(dv[0]/dX[0], numb_divide(dv[1],dX[2]*v), lmbda)
    # get epsilon term
    epsterm = epsilon_term(dv[0]/dX[0], numb_divide(dv[1],dX[2]*v), epsilon)
    # get first term
    first = numba_multiply3(expsemb_v,dexpsemb_dv,numb_add2(lmbdatrm, epsterm))
    # get second term by axis
    second_0 = numb_secondterm(expsemb_v,lmbdatrm, dv[0], epsilon, 0)/dX[0]
    second_1 = numb_divide(numb_secondterm(expsemb_v,lmbdatrm, dv[1], epsilon, 1),dX[2]*v)
 
    return -1*numb_add3(first,second_0,second_1)



def search(u, semb, dX, oX, nX,lmbda, epsilon, direction, a, b, tol=1e-5):
    '''golden search over cost'''
    
    def f(a):
        return compute_cost(numb_force_range(numb_add2(u, -a*direction),oX[1],oX[1]+dX[1]*(nX[1]-1)), 
                            semb, dX, oX,lmbda, epsilon)
    gr = (m.sqrt(5)+1)/2
    c = b-(b-a)/gr
    d = a+(b-a)/gr
    while abs(c-d)>tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b-(b-a)/gr
        d = a+(b-a)/gr
    return (b+a)/2, f(b+a/2)


def get_axis(File,axis):
    o = File.float("o%d"%axis)
    d = File.float("d%d"%axis)
    n = File.int("n%d"%axis)
    return o, d, n 

def put_axis(File,axis,o,d,n):
    File.put("o%d"%axis,o*1.0)
    File.put("d%d"%axis,d*1.0)
    File.put("n%d"%axis,int(n))
    return o, d, n 


# big kahuna
def variational_velocity(semb, dsembv, vo, rho, lmbda, epsilon, dX, oX, nX,
                               niter, mem, rhoscl=1e-1, eps=1e-8,
                               Fmovie=None, search_type='lbfgs',verb=False):




    '''compute the optimal velocity in variational manner'''
    if search_type not in ['lbfgs','line','grad']:
        search_type = 'lbfgs'
        print("Search type not found, changing to LBFGS",file=sys.stderr)
        print("Choices are lbfgs line or grad",file=sys.stderr)
   
    if search_type == 'lbfgs':
        lbfgs_bool = True
    else:
        lbfgs_bool = False
    if search_type == 'grad':
        grad_bool = True
    else:
        grad_bool = False
    if search_type == 'line' or lbfgs_bool:
        line_bool = True
    else:
        line_bool = False
    # compute demb/dv
    #dsembv = gradn(semb,1)/dX[1]
    
    v = numb_force_range(vo,oX[1],oX[1]+dX[1]*(nX[1]-1))

    # compute initial cost
    oldcost = compute_cost(v, semb, dX, oX, lmbda, epsilon)   

    costlst = [oldcost]

    # initialize step memory
    syr_lst = []

    # initialize frame tracker
    frames = 0 # no longer counting starting model
    it1 = 0
    rho0 = rho
    for it in range(niter):
        g = compute_grad(v, semb, dsembv, dX, oX, lmbda, epsilon)
        if lbfgs_bool:
            # nested loop LBFGS
            if it1 > 0:
                q = np.array(g)
                y = numb_add2(g,-g_old)
                s = numb_add2(v,-v_old)
                r = numb_dot(s,y)
                gamma = r / numb_dot(y,y)
                g_old = np.array(g)
                v_old = np.array(v)
                syr_lst.insert(0,[s,y,r])
                while len(syr_lst)>mem:
                    syr_lst.pop(-1)
                alpha = np.zeros(len(syr_lst))
                for jt in range(len(syr_lst)):
                    syr = syr_lst[jt]
                    a_i = numb_dot(q,syr[0])/syr[2]
                    q = q - a_i*syr[1]
                    alpha[jt] = a_i
                # z is newton search direction
                z = gamma*q
                for jt2 in range(len(syr_lst)):
                    jt = len(syr_lst)-1-jt2
                    syr = syr_lst[jt]
                    b_i = numb_dot(z,syr[1])/syr[2]
                    z = numb_add2(z,( alpha[jt]-b_i)*syr[0])
                    # LBFGS has near unit step size for search direction
                if it == 1:
                    eps = eps/rho
                    rho = 1
            else:
                z = np.array(g)
                v_old = np.array(v)
                g_old = np.array(g)
        else:
            z = np.array(g)
        # line search in newton search direction
        if line_bool:
            step, cost = search(v, semb, dX, oX, nX, lmbda, epsilon, z , rho*rhoscl, rho, rho*rhoscl)
            if lbfgs_bool:
               if cost+eps > oldcost:
                  print("Cost of %g Diverging on LBFGS, trying to reset search direction"%cost,file=sys.stderr)
                  # reset 
                  syr_lst = []
                  z = np.array(g)
                  g_old = np.array(g)
                  v_old = np.array(v)
                  it1 = 0
                  rho = rho0
                  step, cost = search(v, semb, dX, oX, nX, lmbda, epsilon, z , rho*rhoscl, rho, rho*rhoscl)
        else:
            step = rho
            cost = compute_cost(numb_add2(v,-rho*z), semb, dX, oX, lmbda, epsilon)   
        if verb:
            print('iteration %i'%(it+1)+', cost %g'%cost,file=sys.stderr)
        # are we descending?
        if oldcost > cost+eps:
            v = numb_force_range(numb_add2(v, -step*z),oX[1],oX[1]+dX[1]*(nX[1]-1))
            oldcost = cost
            costlst += [cost]
            rho = step   
            if Fmovie is not None:
                Fmovie.write(v)
                frames = frames+1
                
        else:
            # make sure we have a non-empty output for updates
            if (Fmovie is not None) and (frames == 0):
                Fmovie.write(v)
                frames = frames+1
            if verb:
                print("Terminating for cost convergence",file=sys.stderr)
            break
        # how are the steps looking
        if rho < eps :
            if verb:
                print("Terminating for step size",file=sys.stderr)
            break
        it1 = it1+1
    return v, costlst, frames#, movie_lst_step, movie_lst_vel


