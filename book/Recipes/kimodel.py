#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# kimodel.py  (Madagascar Recipe)
#
# Purpose: Recipe to Kirchhoff modeling.
# 
# Important!: It should be called from a SConstruct 
#
# Site: http://www.dirackslounge.online
# 
# Version 1.0
#
# Programer: Rodolfo A. C. Neves (Dirack) 04/03/2020
#
# Email: rodolfo_profissional@hotmail.com
#
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

# Selfdoc string
'''
Madagascar recipe to Kirchhoff modeling

Define functions to generate multi layer models and Kirchhoff Newton modeling.
'''

if __name__ == "__main__":
	print(__doc__)
	exit()

# Madagascar package
from rsf.proj import *

__author__="Rodolfo Dirack <rodolfo_profissional@hotmail.com>"
__version__="1.0"

def arr2str(array,sep=' '):
	'''
	
	Convert a tuple into a comma separeted string
	
	:param array: tuple to be converted
	'''
	return string.join(map(str,array),sep)


def multiLayerModel(
	interfaces,
	dipsfile,
	modelfile,
	xmax,
	zmax,
	layers,
	velocities
	):
	'''

	Generate a multi layer model to use in the program sfkirmod_newton

	:out interfaces: RSF filename, interpolated interfaces file
	:out dipsfile: RSF filename, dips of interfaces
	:out modelfile: RSF filename, model for ploting
	:param xmax: interger, max x axis model distance
	:param zmax: interger, max z axis model depth
	:param layers: tuple array, points describing interfaces
	:param velocities: tuple, layers velocities in Km/s
	'''
	vstr=arr2str(velocities,',')
	
	n1 = len(layers[0])
	n2 = len(layers)

	Flow('layers.asc',None,
	     '''
	     echo %s
	     n1=%d n2=%d o1=0 d1=%g
	     data_format=ascii_float in=$TARGET
	     ''' % (string.join(map(arr2str,layers),' '),
		    n1,n2,xmax/(n1-1)))
	Flow('layers','layers.asc','dd form=native')

	d = 0.0101 # non-round for reproducibility

	Flow(interfaces,'layers',
	     'spline o1=0 d1=%g n1=%d' % (d,int(1.5+xmax/d)))
	Flow(dipsfile,interfaces,'deriv scale=y')

	Flow(modelfile,interfaces,
	     '''
	     unif2 d1=%g n1=%d v00=%s
	     ''' % (d,int(1.5+zmax/d),vstr))

def kirchhoffNewtonModeling(
	reflectors, 
	reflectorsDip,
	filename,
	velocities,
	nt,
	dt,
	ns,
	ds,
	nh,
	dh):
	'''
	Kirchhoff modeling for a multi layer model

	:param reflectors: tuple of interfaces, each interface is a tuple of points 
	:param reflectorsDip: RSF filename, interfaces dip
	:param filename: RSF filename, modeled data cube (time x offset x CMP)
	:param velocities: tuple of velocities of each layer
	:param nt: integer, number of time samples
	:param dt: float, time sampling
	:param ns: integer, number of CMP samples
	:param ds: float, CMP sampling
	:param nh: integer, number of offset samples
	:param dh: float, offset sampling
	'''

	# Convert velocity tuple in a comma separated string
	vstr=arr2str(velocities,',')

	# Kirchhoff modeling for multi layer model
	Flow(filename,[reflectors,reflectorsDip],
		'''
		kirmod_newton nt=%i dt=%g freq=10
		ns=%i ds=%g nh=%i dh=%g h0=0 s0=0 verb=y cmp=y
		vstatus=0 velocity=%s debug=n fwdxini=y
		xref=0 zref=0 dip=${SOURCES[1]} |
		put label3="CMP" unit3="Km" label2="Offset" unit2="Km" label1=Time unit1=s
		''' % (nt,dt,ns,ds,nh,dh,vstr))
 
