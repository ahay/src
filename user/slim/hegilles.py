#!/usr/bin/env python

'''
Utilities for Python-based Madagascar programs in RSFROOT/user/slim
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

import numpy as np

def killtraces(N,perc,maxfactor,seed=None):
    '''
    Return a mask to remove random traces from a 2-D gather using a
    maximum gap size constraint.

    Parameters:
    N          number of traces in the 2-D gather
    perc       percentage of traces to remove (0 < perc < 1)
    maxfactor  maximum gap factor (>=1)
    seed       seed for random number generator
    '''
    if seed is not None:
        np.random.seed(seed)

    if maxfactor<1:
        print "ERROR: maxfactor must be >=1"
        return

    if perc<=0 or perc>=100:
        print "ERROR: 0 < perc < 100"
        return        
        
    unifdist = 100./perc
    numtraces = np.int(round(N*perc/100.))
    maxdist = maxfactor * unifdist
    meandist = unifdist + (maxfactor-1.)*unifdist/2.
    randstart = (np.random.rand()*(meandist-unifdist))
    initpick = randstart + np.arange(0,N,meandist)

    randfactor = (np.random.rand(initpick.size)-.5)*(meandist-unifdist)  
    initpick = np.round(initpick + randfactor)
    
    initpick = np.sort(np.unique(initpick.clip(1,N-1)))
    RemainingTrLocs = np.setxor1d(initpick,np.array(range(N)))
    np.random.shuffle(RemainingTrLocs)

    if (numtraces-len(initpick))>0:
        pickinginds = np.union1d(initpick,RemainingTrLocs[:(numtraces-len(initpick))])
    else:
        pickinginds = initpick
    picking = np.zeros(N)
    picking.put(np.ones(len(pickinginds)),list(pickinginds))

    return picking

# $Id$
