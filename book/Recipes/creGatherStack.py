#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# creGatherStack.py  (Madagascar Recipe)
#
# Purpose: Recipe to CRE stacking.
# 
# Important!: It should be called from a SConstruct 
#
# Site: https://dirack.github.io
# 
# Version 1.0
#
# Programmer: Rodolfo A. C. Neves (Dirack) 04/03/2020
#
# Email: rodolfo_profissional@hotmail.com
#
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

# Selfdoc string
'''
Madagascar recipe to CRE stacking

Define functions to generate the Zero-offset section using CRE stacking
'''

if __name__ == "__main__":
	print(__doc__)
	exit()

# Madagascar package
from rsf.proj import *

__author__="Rodolfo Dirack <rodolfo_profissional@hotmail.com>"
__version__="1.0"

def creStack(parametersCube,
            dataCube,
            interpolatedDataCube,
            stackedSection,
            nm0,
            om0,
            dm0,
            nt0,
            ot0,
            dt0,
            v0,
            repeat,
            aperture
            ):
    '''

    CRE stacking

    :out parametersCube: RSF file, RN, RNIP and BETA
    :param dataCube: RSf file, seismic data cube
    :param interpolatedDataCube: RSF file, interpolated seismic data cube
    :out stackedSection: RSF file, stacked section
    :param nm0: int, number of m0s
    :param om0: float, m0s axis origin
    :param dm0: float, m0s sampling
    :param nt0: int, number of t0s
    :param ot0: float, t0s axis origin
    :param dt0: float, t0s sampling
    :param v0: float, near surface velocity
    :param repeat: int, how many times to run VFSA
    :param aperture: int, number of offsets to stack
    '''
    Flow(parametersCube,dataCube,
            '''
            vfsacrsnh nm0=%d om0=%g dm0=%g nt0=%d ot0=%g dt0=%g v0=%g repeat=%d
            '''%(nm0,om0,dm0,nt0,ot0,dt0,v0,repeat))

    Flow('creTrajectories',[interpolatedDataCube,parametersCube],
            '''
            cretrajec param=${SOURCES[1]} nm0=%d om0=%g dm0=%g nt0=%d ot0=%g dt0=%g verb=y 
            '''%(nm0,om0,dm0,nt0,ot0,dt0))

    Flow(['cregathers','mhCoordinates'],[interpolatedDataCube,'creTrajectories'],
            '''
            getcregather cremh=${SOURCES[1]} m=${TARGETS[1]} aperture=%d nm0=%g nt0=%g |
            put label1=Time unit1=s label2=Offset unit2=Km label3=t0 unit3=s
            label4=m0 unit4=Km n3=%d d3=%g o3=%g n4=%d d4=%g o4=%g
            ''' % (aperture,nm0,nt0,nt0,dt0,ot0,nm0,dm0,om0))

    Flow('cretimecurves',['mhCoordinates',parametersCube],
            '''
            getcretimecurve param=${SOURCES[1]} nm0=%d om0=%g dm0=%g nt0=%d ot0=%g dt0=%g verb=y v0=%g |
            put label1=Offset unit1=Km label2=t0 unit2=s label3=m0 unit3=Km
            n2=%d d2=%g o2=%g n3=%d d3=%g o3=%g
            '''%(nm0,om0,dm0,nt0,ot0,dt0,v0,nt0,dt0,ot0,nm0,dm0,om0))

    Flow(stackedSection,['cregathers','cretimecurves'],
            '''
            crestack aperture=%d verb=y timeCurves=${SOURCES[1]} |
            put label1=t0 unit1=s label2=m0 unit2=Km
            ''' %(aperture))

