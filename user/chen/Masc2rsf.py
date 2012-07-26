#!/usr/bin/env python

import sys
try:
	import numpy as np
	import rsf.api as rsf
except Exception, e:
	print 'ERROR : need numpy'
	sys.exit(1)

po = rsf.Output()
par= rsf.Par()
data = np.genfromtxt(sys.stdin)
(n2, n1)=data.shape

po.put('o1',0)
po.put('o2',0)
po.put('d1',1)
po.put('d2',1)
po.put('n1',n1)
po.put('n2',n2)

po.write(data.astype('f'))



