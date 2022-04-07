#!/usr/bin/env python
'''
Write data datumed to windowed point cloud to original point cloud
Created by: Dylan Hickson, Colorado School of Mines
Created on: Mar 9, 2022
'''
import rsf.api as rsf
import numpy as np

par = rsf.Par()

# RSF input file names from parameters
verb = par.bool('verb',False) # verbosity flag
modTinf = par.string('modTin')
topoWf  = par.string('topoW')
topof   = par.string('topo')

# Local functions
def read_rsf2D(rsffile):
	Fin = rsf.Input(rsffile)

	n1 = Fin.int('n1') # number of coordinates (9)
	n2 = Fin.int('n2') # number of points

	data = np.zeros((n1,n2),dtype=float)

	# Read data to numpy
	j2 = 1
	dn1 = np.zeros(n1,'f')
	i=0
	while i < n2:
		Fin.read(dn1)
		data[:,i] = dn1[:]
		i+=j2

	Fin.close()
	return data

# Read data into numpy arrays
topo   = read_rsf2D(topof)
topoW  = read_rsf2D(topoWf)
modTin = read_rsf2D(modTinf)

# Get header parameters from input data
modTin  = read_rsf2D(modTinf)
Fin = rsf.Input(modTinf)
ng = topo.shape[1]
nt = Fin.int('n1')
dt = Fin.float('d1')
ot = Fin.float('o1')

# Define header parameters for output data
Fou = rsf.Output()
Fou.put('n1',nt)
Fou.put('d1',dt)
Fou.put('o1',ot)
Fou.put('n2',ng)
Fou.put('d2',1)
Fou.put('o2',0)

for i in range(topo.shape[1]):
	if topo[0,i] in topoW[0,:] and topo[1,i] in topoW[1,:] and topo[2,i] in topoW[2,:]:
		xind = set(np.nonzero(topoW[0,:]==topo[0,i])[0])
		yind = set(np.nonzero(topoW[1,:]==topo[1,i])[0])
		zind = set(np.nonzero(topoW[2,:]==topo[2,i])[0])

		intersection = list(xind.intersection(yind, zind))

		if len(intersection) > 0:
			modTDat = np.array(modTin[:,intersection[0]])
		else:
			modTDat = np.zeros(nt,'f')
	else:
		modTDat = np.zeros(nt,'f')
	Fou.write(modTDat)

#------------------------------------------------
Fin.close()
Fou.close()
