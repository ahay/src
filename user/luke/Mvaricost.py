#!/usr/bin/env python



'''
Determines costs associated with velocity models from varipick.
Input is velocity model, needs semb= semblance volume

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
Fvel  = m8r.Input()
Fsemb = m8r.Input("semb")



assert 'float' == Fsemb.type
assert 'float' == Fvel.type

# get axis sampling
o1,d1,n1 = vt.get_axis(Fsemb,1)
o2,d2,n2 = vt.get_axis(Fsemb,2)
o3,d3,n3 = vt.get_axis(Fsemb,3)
# get number of velocities
ov,dv,nv = vt.get_axis(Fvel,3)

# put in helpful arrays
dX = np.array([d3,d2,d1])
oX = np.array([o3,o2,n1])
nX = np.array([n3,n2,n1])

# read semblance and velocity
s = np.zeros(n3*n2*n1,'f')
s = Fsemb.read(shape=(n3,n2,n1))
v = np.zeros(n3*n1,'f')

# get input parameters
lmbda = par.float('lambda')
assert lmbda
epsilon = par.float('epsilon')
assert epsilon


costlst = []
for i in range(nv):
    v = Fvel.read(shape=(n3,n1))
    costlst.append(vt.compute_cost(vt.numb_force_range(v,oX[1],oX[1]+dX[1]*(nX[1]-1)),s, dX, oX,lmbda,epsilon))


costs = np.asarray(costlst)


# write output
Fout = m8r.Output()
# set output sampling
vt.put_axis(Fout,1,0,1,i+1)
vt.put_axis(Fout,2,0,1,1)
vt.put_axis(Fout,3,0,1,1)
Fout.write(costs)
Fout.close()

Fsemb.close()
Fvel.close()


                
