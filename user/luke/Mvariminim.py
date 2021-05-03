#!/usr/bin/env python
'''
outputs the minimum and maximum cost velocity models
input is series of models
semb= semblance volume
costs= output for cost
worst= worst model output
primary output is lowest cost model

'''
# key imports
import m8r
import numpy as np
try:
    import varitools as vt
except:
    import rsf.user.varitools as vt
import sys
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

# initialize best velocity and array
best = np.zeros(n3*n1,'f')
bestcost = -9999
bestindex = 0
# and worst index and cost
worstcost = -9999
worstindex = 0
worst = np.zeros(n3*n1,'f')


costlst = []
# and mean cost
for i in range(nv):
    v = Fvel.read(shape=(n3,n1))
    cost = vt.compute_cost(vt.numb_force_range(v,oX[1],oX[1]+dX[1]*(nX[1]-1)),s, dX, oX,lmbda,epsilon)
    costlst.append(cost)
    # if better cost
    if (cost < bestcost) or (i == 0):
        # set best cost
        bestcost = cost
        # set array to best
        best = np.array(v)
        # record best index
        bestindex = i
    # if worst cost
    if (cost > worstcost) or (i == 0):
        worstcost = cost
        worstindex = i
        worst = np.array(v)
print('Best model was %i with cost of %g'%(bestindex,bestcost),file=sys.stderr)
print('Worst model was %i with cost of %g'%(worstindex,worstcost),file=sys.stderr)
# write output
Fout = m8r.Output()
# set output sampling
vt.put_axis(Fout,3,0,1,1)
Fout.write(best)


# record all costs if desired
costname = par.string('costs')
if costname is not None:
    Fcost = m8r.Output(costname)
    vt.put_axis(Fcost,1,0,1,i+1)
    vt.put_axis(Fcost,2,0,1,1)
    vt.put_axis(Fcost,3,0,1,1)
    costs = np.asarray(costlst)
    Fcost.write(costs)
    Fcost.close()

worstname = par.string('worst')
if worstname is not None:
    Fworst = m8r.Output(worstname)
    # set output sampling
    vt.put_axis(Fworst,3,0,1,1)
    Fworst.write(worst)
    Fworst.close()

Fout.close()

Fsemb.close()
Fvel.close()


                
