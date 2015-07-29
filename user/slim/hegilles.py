#!/usr/bin/env python

'''
Utilities for Python-based Madagascar programs in RSFSRC/user/slim
'''

# Author: G. Hennenfent
#         Seismic Laboratory for Imaging and Modeling
#         Department of Earch & Ocean Sciences
#         The University of British Columbia
#         
# Date  : February, 07

#  Copyright (C) 2006 The University of British Columbia at Vancouver
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy as __np

def killtraces(N,perc=50,maxfactor=1,seed=None):
    '''
    Return a mask to remove random samples from an N-vector using a
    maximum gap size constraint.

    Parameters:
    N          length of the vector
    perc       percentage of samples to retain (0 < perc < 100)
    maxfactor  maximum gap factor (>=1)
    seed       seed for random number generator
    '''
    if seed is not None:
        __np.random.seed(int(seed))

    if maxfactor<1:
        print "ERROR: maxfactor must be >=1"
        return

    if perc<=0 or perc>=100:
        print "ERROR: 0 < perc < 100"
        return        
        
    unifdist = 100./perc
    numtraces = __np.int(round(N*perc/100.))
    maxdist = maxfactor * unifdist
    meandist = unifdist + (maxfactor-1.)*unifdist/2.
    randstart = (__np.random.rand()*(meandist-unifdist))
    initpick = randstart + __np.arange(0,N,meandist)

    randfactor = (__np.random.rand(initpick.size)-.5)*(meandist-unifdist)  
    initpick = __np.round(initpick + randfactor)
    
    initpick = __np.sort(__np.unique(initpick.clip(1,N-1)))
    RemainingTrLocs = __np.setxor1d(initpick,__np.array(range(N)))
    __np.random.shuffle(RemainingTrLocs)

    if (numtraces-len(initpick))>0:
        pickinginds = __np.union1d(initpick,RemainingTrLocs[:(numtraces-len(initpick))])
    else:
        pickinginds = initpick
    picking = __np.zeros(N)
    picking.put(list(pickinginds),__np.ones(len(pickinginds)))

    return picking

def jitter(N,perc=50.,jit=None,seed=None):
    '''
    Return a mask to remove random samples from an N-vector using
    jittered sampling

    Parameters:
    N          length of the vector
    perc       percentage of samples to retain (0 < perc < 100)
    seed       seed for random number generator
    '''
    if seed is not None:
        __np.random.seed(int(seed))
    
    if perc<=0 or perc>=100:
        print "ERROR: 0 < perc < 100"
        return        

    sub = 100/perc

    if jit is None:
        jit = sub

    y = __np.zeros(N)
    inipos = __np.arange(0,N,sub)
    pos = __np.array(__np.round(inipos +__np.float(jit)*__np.random.rand(inipos.size)-jit/2.)%N,dtype='i')
    y.put(pos,__np.ones(pos.size))

    return y
    
# $Id: hegilles.py 5990 2010-05-17 22:11:45Z sfomel $
