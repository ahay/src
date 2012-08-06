'extra flows'

import string, re, sys
from rsf.proj import *

if sys.platform == 'darwin':
	timer=WhereIs('gtime')
else:
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
			>& ${TARGETS[0]}.tmp &&
		(tail -1 ${TARGETS[0]}.tmp; 
		echo in=${TARGETS[0]}.tmp.rsf n1=1 data_format=ascii_float) 
			> ${TARGETS[0]}.tmp.rsf &&
		dd form=native <${TARGETS[0]}.tmp.rsf > ${TARGETS[0]} &&
		%s -f ${TARGETS[0]}.tmp.rsf ${TARGETS[0]}.tmp
		'''%(timer, cmd, WhereIs('rm')),
		 stdin=None, stdout=None)


def Mplot(target, cmd, mtx):
	Result(target,
		cmd+''' 
		| %s serifs=y fat=3 color=y label="" tex=y scale=0.8
		| ps2pdf - -
		| pdfjam -q -o /dev/stdout --papersize "{5in,5in}"
		  --trim '6.2cm 2.2cm 7.2cm 7.5cm'
		| pdfjam -q -o /dev/stdout --papersize "{5in,5in}"
		  --nup %s  
		'''%(WhereIs('pspen'),mtx), suffix='.pdf')


