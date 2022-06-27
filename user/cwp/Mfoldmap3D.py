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
def read_rsf(fname):
	fIn = rsf.Input(fname+'.rsf')
	data = fIn.read()
	fIn.close()
	return data.T

topo = read_rsf(topof)

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
	topoW = read_rsf(file)

	for i in range(topo.shape[1]):
		if foldmap[0,i] in topoW[0,:] and foldmap[1,i] in topoW[1,:] and foldmap[2,i] in topoW[2,:]:
			xind = set(np.nonzero(topoW[0,:]==foldmap[0,i])[0])
			yind = set(np.nonzero(topoW[1,:]==foldmap[1,i])[0])
			zind = set(np.nonzero(topoW[2,:]==foldmap[2,i])[0])

			intersection = list(xind.intersection(yind, zind))

			if len(intersection) > 0:
				foldmap[3,i] += 1

# Change all instances where fold number = 0 to 1 to avoid division by 0
foldmap[foldmap==0]=1

# Write fold map
for i in range(topo.shape[1]):
	writeOut = foldmap[:,i]
	Fou.write(writeOut)

#------------------------------------------------
Fou.close()
