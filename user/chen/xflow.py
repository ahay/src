'extra flows'

import string, re
from rsf.proj import *

timer=WhereIs('time')

def Tflow(target, source, command,
	prefix='sf',):
	if type(target) is types.ListType:
		tfiles = target
	else:
		tfiles = string.split(target)
	tfiles.insert(0, tfiles[0]+'_runtime')
	pars=command.split()
	p0=pars.pop(0)

	p1=WhereIs(p0) # sfdip
	if p1==None :
		p1=WhereIs(prefix+p0) # dip
	if re.match(r'[^/]+\.exe$',p0) :
		p1=os.path.join('.',p0) # Mdip.exe

	pars.insert(0,p1)
	cmd=string.join(pars,' ')
	Flow(tfiles, source,
		'''
		( %s -f "%%U" %s <${SOURCES[0]} >${TARGETS[1]} ) 
			>& tmp_timer.txt &&
		(tail -1 tmp_timer.txt; 
		echo in=tmp_timer.rsf n1=1 data_format=ascii_float) 
			> tmp_timer.rsf &&
		dd form=native <tmp_timer.rsf > ${TARGETS[0]} &&
		%s -f tmp_timer.rsf tmp_timer.txt
		'''%(timer, cmd, WhereIs('rm')),
		 stdin=None, stdout=None)


