#!/usr/bin/env python
'''
Create a foldmap from multiple orbit swaths
Created by: Dylan Hickson, Colorado School of Mines
Created on: Mar 9, 2022
'''
import rsf.api as rsf
import numpy as np

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
topof   = par.string('topo')
topoWstr = par.string('topoWl')

topoWlist = topoWstr.split(',')

# Local functions
def read_rsf2D(rsffile):
	Fin = rsf.Input(rsffile)

	n1 = Fin.int('n1') # number of coordinates
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

topo = read_rsf2D(topof)

foldmap = np.zeros((4,topo.shape[1]))
foldmap[0:3,:] = topo[0:3,:]

Fou = rsf.Output()
Fou.put('n1',4)
Fou.put('d1',1)
Fou.put('o1',0)
Fou.put('n2',topo.shape[1])
Fou.put('d2',1)
Fou.put('o2',0)

for file in topoWlist:
	topoW = read_rsf2D(file)

	for i in range(topo.shape[1]):
		if foldmap[0,i] in topoW[0,:] and foldmap[1,i] in topoW[1,:] and foldmap[2,i] in topoW[2,:]:
			xind = set(np.nonzero(topoW[0,:]==foldmap[0,i])[0])
			yind = set(np.nonzero(topoW[1,:]==foldmap[1,i])[0])
			zind = set(np.nonzero(topoW[2,:]==foldmap[2,i])[0])

			intersection = list(xind.intersection(yind, zind))

			if len(intersection) > 0:
				foldmap[3,i] += 1

# Write fold map
for i in range(topo.shape[1]):
	writeOut = foldmap[:,i]
	Fou.write(writeOut)

#------------------------------------------------
Fou.close()
