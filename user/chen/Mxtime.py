#!/usr/bin/env python
'runtime'

import sys
import os
try:
	import numpy as np
	import rsf.api as rsf
except Exception, e:
	print 'ERROR : need numpy'
	sys.exit(1)

plat=sys.platform
if plat=='darwin': # MAC OS using GNU gtime
	time=' gtime -f "%U" 2>&1 '
else: # RHEL 
	time=' time -f "%U" 2>&1 '

par=rsf.Par()
po=rsf.Output()

strcmd=par.string('cmd')
if strcmd == None:
	sys.stderr.write('cmd must be specified')
	sys.exit(1)
lstcmd=strcmd.split(':')

opt=par.string('opt')

# sys.stderr.write(strcmd+' \n'+opt+'\n')

num=len(lstcmd)
po.put("n1",num)
po.put("o1",0)
po.put("d1",1)

b1=np.zeros(num,'f')

for ii in range(num):
	cmd=time+lstcmd[ii]+opt
	aa=os.popen(cmd).readlines()[-1]
	#sys.stderr.write(aa)
	b1[ii]=float(aa)

po.write(b1)

