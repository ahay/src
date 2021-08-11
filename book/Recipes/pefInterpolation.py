#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# pefInterpolation.py  (Madagascar Recipe)
#
# Purpose: Recipe to Predictive Adaptative Error Filters (PEF) interpolation
# to double CMP number of samples in a section or data cube.
#
# Site: https://dirack.github.io
# 
# Version 1.0
#
# Programer: 
#		-Saulo S Martins (ssmsaulo) 08/12/2016 - Article reference (2015)
#		 	http://dx.doi.org/10.1190/segam2015-5924518.1
#		 Email: ssmsaulo@gmail.com
#
#		-Hugo Souza (hugo-ss) 22/06/2019 - Original SConstruct code
#		 Email: hugogeof@gmail.com
#
#		-Rodolfo A. C. Neves (Dirack) 27/05/2020 - This recipe
#		 Email: rodolfo_profissional@hotmail.com
#
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

# Selfdoc string
'''
Madagascar recipe to PEF interpolation

Interpolate and increase CMP sampling using Predictive Adaptative Error Filters (PEF) interpolation.
'''

if __name__ == "__main__":
	print(__doc__)
	exit()

__author__="Rodolfo Dirack <rodolfo_profissional@hotmail.com>"
__version__="1.0"

# Madagascar package
from rsf.proj import *

def pefInterpolation(
    dataCube,
    interpolated,
    nm,
    dm,
    nt,
    dt,
    nhi,
    a1,
    a2,
    rect1,
    rect2
    ):
    '''

    PEF interpolation to double CMP sampling of the data cube

    :param dataCube: RSF filename, Seismic data cube to interpolate
    :param interpolated: RSF filename, Interpolated seismic data cube
    :param nm: integer, number of CMPs in the seismic data cube
    :param dm: float, CMP sampling
    :param nt: integer, number of time samples
    :param dt: float, time sampling
    :param nhi: integer, number of constant offsets gathers to interpolate
    :param a1: integer, Number of PEF coeficients in time axis
    :param a2: integer, Number of PEF coeficients in space axis
    :param rect1: integer, Smooth radius in time
    :param rect2: integer, Smooth radius in space
    '''

    # Number of offset gathers should be positive and non zero
    assert nhi>=1, 'Number of offset gathers (nhi) should be >= 1, you used (nhi=%i)' % nhi

    # Divide CMP sampling
    dm = dm/2

    # Define mask file names using input filename
    mask1 = dataCube+'-mask1'
    mask = dataCube+'-mask'
    aa = dataCube+'-aa'
    bb = dataCube+'-bb'
    a = dataCube+'-a'
    b = dataCube+'-b'
    zeroTraceGather = dataCube+'-zeroedGather'
    mask0 = dataCube+'-mask0'


    # Build a mask to interleave zero traces with original data traces
    Flow(aa,None,'spike n1=%i d1=%g o1=0' %(nm,dm))
    Flow(bb,None,'spike n1=%i d1=%g o1=0 mag=0' % (nm,dm))
    Flow(mask1,[bb, aa],
            '''
            interleave axis=1 ${SOURCES[1]} |
            dd type=int
            ''')

    Flow(a,None,'spike n1=%i d1=%g o1=0' % (nm,dm))
    Flow(b,None,'spike n1=%i d1=%g o1=0 mag=0' % (nm,dm))
    Flow(mask,[a, b],
            '''
            interleave axis=1 ${SOURCES[1]} |
            dd type=int
            ''')
    Flow(zeroTraceGather,b,
            '''
            spray axis=2 n=%i d=%g |
            transp |
            put label2=Offset unit2=Km label1=Time unit1=s
            ''' %(nt,dm))

    # Data Mask with double of traces in CMP (half of CMP sampling)
    # Keep the same Time and Offset original data sampling
    Flow(mask0,mask,
         '''
         spray axis=1 n=%i d=%g
         ''' %(nt,dt))

    totalPefIterations = 100
    totalInterpolationIterations = 50

    offsetGathers = []
    for offsetGatherIndex in range(nhi):

            offsetGather = dataCube+"-offsetGather-%i" % offsetGatherIndex
            resampledOffsetGather = dataCube+"-resampledGather-%i" % offsetGatherIndex
            interpolatedOffsetGather = dataCube+"-interpolatedGather-%i" % offsetGatherIndex
            pefCoeficients = dataCube+"-pefCoeficients-%i" % offsetGatherIndex

            Flow(offsetGather,dataCube,
            '''
            window n2=1 f2=%i
            ''' % (offsetGatherIndex))
            
            Flow(resampledOffsetGather,[offsetGather,zeroTraceGather],
            '''
            interleave axis=2 ${SOURCES[1]}
            ''')

            # Calculate adaptive PEF coeficients
            Flow(pefCoeficients,[resampledOffsetGather,mask0],
                    '''
                    apef jump=2 a=%i,%i rect1=%i rect2=%i niter=%g verb=y
                    maskin=${SOURCES[1]}
                    ''' % (a1,a2,rect1,rect2,totalPefIterations))

            # Interpolation
            Flow(interpolatedOffsetGather, 	[resampledOffsetGather,pefCoeficients,mask0,mask1],
                    '''
                    miss4 exact=y filt=${SOURCES[1]} mask=${SOURCES[2]} niter=%g verb=y |
                    put d2=%g
                    ''' % (totalInterpolationIterations,dm))

            offsetGathers.append(interpolatedOffsetGather)

    # Concatenate interpolated sections
    Flow(interpolated,offsetGathers,
            '''
            rcat axis=3 ${SOURCES[1:%d]} |
            transp plane=23
            ''' % nhi)


