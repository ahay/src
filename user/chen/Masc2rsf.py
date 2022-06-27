#!/usr/bin/env python

##   Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
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

from __future__ import print_function
import sys
try:
	import numpy as np
	import rsf.api as rsf
except:
	print('ERROR : need numpy')
	sys.exit(1)

po = rsf.Output()
par= rsf.Par()
data = np.genfromtxt(sys.stdin)

nk=len(data.shape)

if nk==2:
	(n2, n1)=data.shape
	po.put('o2',0)
	po.put('d2',1)
	po.put('n2',n2)
elif nk==1:
	n1=data.shape[0]

po.put('o1',0)
po.put('d1',1)
po.put('n1',n1)

po.write(data.astype('f'))



