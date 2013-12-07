#!/usr/bin/env python
'testmef'

##   Copyright (C) 2013 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


import numpy
import rsf.api as rsf
from pylab import*

# initialize Madagascar parameters
par = rsf.Par()
 
max_dip_angle=par.float("maxangle",20) # maximum dip angle

time_sample_number=100
mask_width=3./4

mlength=array(range(1,time_sample_number+1))

a=time_sample_number
b=floor(mask_width*max_dip_angle)

width=2*max_dip_angle+1

mleft=floor(b*(1-(mlength)**2./a**2)**.5)
mleft=mleft.astype(int)
mright=-mleft+width

mask=numpy.zeros((time_sample_number,width),'f')
    
for i in range(time_sample_number):
    mask[i,mleft[i]:mright[i]]=1.

mask=mask[::-1]


# output to RSF
fmask=rsf.Output()
fmask.put("n2",time_sample_number)
fmask.put("n1",width)
fmask.write(mask)
