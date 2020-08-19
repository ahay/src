#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# velocityAnalysis.py (Madagascar Recipe)
# 
# Purpose: Madagascar recipe to do velocity analysis, NMO
# correction and stack.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 01/04/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

from rsf.proj import *
import math

# Selfdoc string
'''
Madagascar recipe to do velocity analysis with automatic velocity scan and picking

Define functions to do a velocity analysis of a (time,offset,CMP) datacube.
'''

if __name__ == "__main__":
	print(__doc__)
	exit()

# Madagascar package
from rsf.proj import *

__author__="Rodolfo Dirack <rodolfo_profissional@hotmail.com>"
__version__="1.0"

def velocityAnalysis(
	dataCube,
	pick,
	stack,
	vrms,
	v0,
	dv,
	nv,
	vel0,
	rect1=15,
	rect2=40,
	rect3=3,
	dt=0.004):
	'''
	
	Velocity analysis and automatic picking and semblance

	:param datacube: RSF filename, (time, offset, CMP) datacube
	:out pick: RSF filename, picking velocity
	:out stack: RSF filename, stacked section
	:out vrms: RSF filename, vrms velocity
	:param dt: float, time sampling
	:param v0: float, inicial scan velocity
	:param dv: float, scan velocity sampling
	:param nv: int, number of scan velocities
	:param vel0: float, initial velocity picking 
	:param rect1: int, smooth time radius to velocity scan
	:param rect2: int, smooth space radius to velocity scan
	:param rect3: int, smooth time radius post picking
	'''
	
	scan = 'scan'+dataCube
	nmo = 'nmo'+dataCube

	# Velocity scan
	Flow(scan,[dataCube],
	     '''
	     vscan semblance=y
	     v0=%g dv=%g nv=%d
	     ''' % (v0,dv,nv))

	# Velocity picking
	Flow(pick,scan,
	     '''
	     pick rect1=%i rect2=%i vel0=%g smooth=y |
	     smooth rect1=%i
	     ''' %(rect1,rect2,vel0,rect3))

	# NMO and stack
	Flow(nmo,[dataCube,pick],
	     'nmo half=n velocity=${SOURCES[1]}')

	Flow(stack,nmo,'stack')

	# generate RMS velocity :
	# sfmul multiply input data by itself to make v**2
	# sfcausint integrate v**2 (this does not include scaling by dt)
	# sfmath scale by dt and divide by t
	Flow(vrms,pick,
	     '''
	     sfmul $SOURCE | 
	     sfcausint |
	     sfmath output="sqrt(input*%g/(x1+%g))" |
	     sfput n3=1 d3=1 o3=0 
	     '''%(dt,dt))
