#!/usr/bin/env python

'''
Performs automatic picking of a surface from a semblance-like volume.
Input is semblance volume, also needs vo= starting model, dsemb= partial derivative of semblance with respect to v
'''


# key imports
import m8r
import numpy as np

try:
    import varitools as vt
except:
    import rsf.user.varitools as vt

# the actual program...
par = m8r.Par()

# files
Fsemb = m8r.Input()
# semblance volume
Fvel  = m8r.Input("vo")
# starting model
Fdsemb  = m8r.Input("dsemb")
# partial derivative of semblance volume with respect to v

assert 'float' == Fsemb.type
assert 'float' == Fvel.type
assert 'float' == Fdsemb.type

# get axis sampling
o1,d1,n1 = vt.get_axis(Fsemb,1)
o2,d2,n2 = vt.get_axis(Fsemb,2)
o3,d3,n3 = vt.get_axis(Fsemb,3)

# get input parameters
lmbda = par.float('lambda')
# positive regularizaton parameter
assert lmbda
niter = par.int('niter')
# number of iterations
assert niter
rho = par.float('rho')
# step size limit
assert rho
epsilon = par.float('epsilon')
# positive regularization parameter for gradient term
assert epsilon
search_type = par.string('type')
#lbfgs, line, or grad
# put in helpful arrays
dX = np.array([d3,d2,d1])
oX = np.array([o3,o2,n1])
nX = np.array([n3,n2,n1])

# read semblance and velocity
s = np.zeros(n3*n2*n1,'f')
s = Fsemb.read(shape=(n3,n2,n1))
ds = np.zeros(n3*n2*n1,'f')
ds = Fdsemb.read(shape=(n3,n2,n1))
vo = np.zeros(n3*n1,'f')
vo = Fvel.read(shape=(n3,n1))

# make sure velocity is in range
vout = vt.numb_force_range(vo,oX[1],oX[1]+dX[1]*(nX[1]-1))

#lmbda = 1#10
#niter = 20
#rho = 0.001
rhoscl = 1e-2
eps = 1e-12
mem = 3
verb = True
#epsilon = 0.001

#Fmovie = m8r.Output("updates")


moviename = par.string('updates')
if moviename is not None:
       Fmovie = m8r.Output(moviename)
       vt.put_axis(Fmovie,2,oX[0],dX[0],nX[0])
       vt.put_axis(Fmovie,3,0,1,niter)
#       vt.put_axis(Fmovie,3,0,1,niter+1) no longer including starting model
#       Fmovie.write(vo) no longer including starting model in updates
else:
       Fmovie = None
#Fmovie.write(vo)


#, movie_lst_step, movie_lst_vel
vout, costlst, frames = vt.variational_velocity(s, ds, vo, rho, lmbda, epsilon,
                                            dX, oX, nX, niter, mem, 
                                            rhoscl, eps, Fmovie ,search_type,verb)

vout = vt.numb_force_range(vout,oX[1],oX[1]+dX[1]*(nX[1]-1))

if moviename is not None:
#    if frames != niter+1: no longer including staring model in updates
    if frames != niter:
        moviefile = open(moviename,'a')
        moviefile.write("\n\tn3=%i\n"%(frames))
        moviefile.close()

# write output
#vt.put_axis(Fmovie,3,0,1,niter)
Fout = m8r.Output()
# set output sampling
vt.put_axis(Fout,2,o3,d3,n3)
vt.put_axis(Fout,3,0,1,1)
Fout.write(vout)
Fout.close()

Fsemb.close()
Fvel.close()


                
