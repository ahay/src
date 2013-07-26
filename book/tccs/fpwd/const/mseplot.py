#!/usr/bin/env python
'variance'

import sys
import rsf.api as rsf

try:
	import numpy as np
except Exception, e:
	print 'ERROR : need numpy'
	sys.exit(1)

try:
	from pylab import *
	from matplotlib.colors import LogNorm
except Exception, e:
	print 'ERROR : need numpy'
	sys.exit(1)

pi = rsf.Input()
po = rsf.Output()
par= rsf.Par()

dim=pi.shape()

n1=pi.int("n1")
o1=pi.float("o1")
d1=pi.float("d1")
n2=pi.int("n2")
o2=pi.float("o2")
d2=pi.float("d2")


# x label
label1=par.string("label1")
if label1 == None:
	label1=pi.string("label1")
if label1 != None:
	label=label1
	unit1=par.string("unit1",None)
	if unit1 == None:
		unit1=pi.string("unit1")
	if unit1 != None:
		label =label+' ('+ unit1+')'
	ylabel(label)

# y label
label2=par.string("label2")
if label2 == None:
	label2=pi.string("label2")
if label2 != None:
	label=label2
	unit2=par.string("unit2",None)
	if unit2 == None:
		unit2=pi.string("unit2")
	if unit2 != None:
		label =label+' ('+ unit2+')'
	xlabel(label)

# bar label
barlabel=par.string("label")
if barlabel == None:
	barlabel=pi.string("label")
if barlabel != None:
	barunit=par.string("unit",None)
	if barunit == None:
		barunit=pi.string("unit")
	if barunit != None:
		barlabel =barlabel+' ('+ barunit+')'

font = {'family' : 'serif'}
rc('font', **font)

#x = array(range(2,11)+range(20,101,10)+range(200,501,100))
y = np.linspace(o1, o1+d1*(n1-1), n1)
x = np.linspace(o2, o2+d2*(n2-1), n2)

z = np.zeros((n2,n1),'f')
pi.read(z)

vmin=par.float("vmin", z.min())
vmax=par.float("vmax", z.max())
gray()
pcolormesh(x, y, z.T, vmin=vmin, vmax=vmax)
semilogx()
ax=gca()
ax.xaxis.set_label_position('top') 
ax.xaxis.set_ticks_position('top') 

cb=colorbar()
cb.ax.set_ylabel(barlabel)

xticks( (2,3,10,100,200),('2','3','10','100','200') )

xlim(o2, o2+(n2-1)*d2)
ylim(o1+(n1-1)*d1, o1)

savefig(sys.stdout, format='eps', bbox_inches='tight', transparents=True)

