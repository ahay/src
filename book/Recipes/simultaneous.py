import os
from rsf.proj import *
from math import pi

def SSreceiver(rr, zg=2, ng=100, og=0, dg=1):
	Flow(rr+'_z', None, 
		'math n1=%d d1=1 o1=0 output="%d"' % (ng, zg))
	Flow(rr+'_x', rr+'_z',
		'math output="x1*%d+%d"'%(dg, og))
	Flow(rr, [rr+'_z', rr+'_x'],
		'''
		cat ${SOURCES[1]} axis=2
		| transp plane=12 
		| put label2="Receivers"
		| dd type=int
		''')

def SSsource(ss, zs=1, ns=1, os=0, ds=10):
	Flow(ss+'_z', None, 
		'math n1=%d d1=1 o1=0 output="%d"' %(ns, zs))
	Flow(ss+'_x', ss+'_z',
		'math output="x1*%d+%d"'%(ds, os))
	Flow(ss, [ss+'_z', ss+'_x'],
		'''
		cat ${SOURCES[1]} axis=2
		| transp plane=12
		| put label2="Shots"
		| dd type=int
		''')

def SSdelay(dd, st=100, et=4000, n1=10, n2=5, var=0.5, seed=2012 ):
	Flow(dd+'_l', None,
		'''
		math n1=%d o1=0 d1=1 n2=%d output="x1*%g+%d"
		| put n1=%d n2=1
		'''%(n1, n2, (et-et/n1*var-st)/n1, st, n1*n2))
	Flow(dd+'_v', None,
		'''
		signal n=%d o=0 d=1 waveform=rand para=%d
		| sfmath output="input*%g"
		'''%(n1*n2, seed, et/n1*var))
	Flow(dd, [dd+'_l', dd+'_v'],
		'''
		math var=${SOURCES[1]} output="input+var"
		| dd type=int
		''')

def SSfd(vel, rr, ss, dd, 
	freq=10, st=200, et=4000, nt=4000, dt=0.001, jt=4, jtm=100,
	prefix="", wfl=""):

	Flow(prefix+'src', dd, 'ss o1=%g d1=%g n1=%g '%(-st*dt, dt, nt+et))
	Flow(prefix+'wvlt', prefix+'src', ' ricker1 frequency=%g'%(freq))
	#Result(prefix+'wvlt', 'grey pclip=100 title="Simultaneous Sources"')

	Flow([prefix+'dat'], [prefix+'wvlt', vel, ss, rr],
		'''
		wavmod jt=%d jtm=%d verb=y
		vel=${SOURCES[1]} sgrid=${SOURCES[2]} ggrid=${SOURCES[3]}
		'''%(jt, jtm))
	Flow(prefix+'crg', [prefix+'dat', dd],
		'sscrg delay=${SOURCES[1]} jt=%d nt=%d '%(jt, nt/jt))
	Flow(prefix+'csg', prefix+'crg',
		'transp plane=23 ')

	#if wfl != "":
	#	Result(wfl, 'grey title="Wave field" gainpanel=all ');
	Result(prefix+'crg', 'grey title="CRG"');
	Result(prefix+'csg', 'grey title="CSG"');



def SSfdbd(vel, rr, seed=2012,
	n1s=5, n2s=2, n3s=10, ds=10, os=500, 
	freq=10, st=200, et=4000, nt=4000, dt=0.001, jt=4, jtm=100,
	prefix="", wfl=""):
	l1=[]
	l2=[]
	for i3s in range(n3s):
		ss=prefix+'ss%d'%i3s
		dd=prefix+'dd%d'%i3s
		wvlt=prefix+'wvlt%d'%i3s
		dat=prefix+'dat%d'%i3s
		csg=prefix+'csg%d'%i3s

		SSsource(ss, ns=n1s*n2s, os=os+i3s*ds, ds=ds*n3s)
		SSdelay(dd, st=st, n1=n1s, n2=n2s, var=0.8, et=et, seed=seed+i3s)

		Flow(wvlt, dd, 
			'''
			ss o1=%g d1=%g n1=%g 
			| ricker1 frequency=%g
			'''%(-st*dt, dt, nt+et, freq))

		Flow(dat, [wvlt, vel, ss, rr],
			'''
			wavmod jt=%d jtm=%d verb=y
			vel=${SOURCES[1]} sgrid=${SOURCES[2]} ggrid=${SOURCES[3]}
			'''%(jt, jtm))
		Flow(csg, [dat, dd],
			'''
			sscrg delay=${SOURCES[1]} jt=%d nt=%d
			| transp plane=23
			'''%(jt, nt/jt))
		l1.append(csg)
		l2.append(dd)

	Flow(prefix+'csg', l1,
		'''
		cat axis=4 ${SOURCES[1:%d]} 
		| transp plane=34 
		| put n4=1 n3=%d 
		'''%(n3s, n1s*n2s*n3s))

	Flow(prefix+'crg', prefix+'csg',
		'transp plane=23 ')

	Result(prefix+'crg', 'grey title="CRG"');
	Result(prefix+'csg', 'grey title="CSG"');


def SSblend(vel, rr, ss, dd, 
	freq=10, st=200, et=4000, nt=4000, dt=0.001, jt=4, jtm=100,
	prefix="", wfl=""):

	Flow(prefix+'src', dd, 'ss o1=%g d1=%g n1=%g '%(-st*dt, dt, nt+et))
	Flow(prefix+'wvlt', None, 
		'signal o=%g d=%g n=%d para=%g '%(-st*dt, dt, nt, freq))
	Result(prefix+'wvlt', 'graph title="Source"')

	Flow([prefix+'shots', wfl], [prefix+'wvlt', vel, ss, rr],
		'''
		wavmod wfl=${TARGETS[1]} jt=%d jtm=%d verb=y
		vel=${SOURCES[1]} sgrid=${SOURCES[2]} ggrid=${SOURCES[3]}
		'''%(jt, jtm))
	Flow(prefix+'dat', [prefix+'shots', dd],
		'''
		transp plane=23
		| ssblend delay=${SOURCES[1]} jt=%d 
		'''%jt)
	Flow(prefix+'crg', [prefix+'dat', dd],
		'sscrg delay=${SOURCES[1]} jt=%d nt=%d '%(jt, nt/jt))
	Flow(prefix+'csg', prefix+'crg',
		'transp plane=23 ')
	Flow(prefix+'reciev', prefix+'shots',
		'transp plane=23 ')

	#if wfl != "":
	#	Result(wfl, 'grey title="Wave field" gainpanel=y ');
	Result(prefix+'shots', 'grey title="Clean CSG"');
	Result(prefix+'reciev', 'grey title="Clean CRG"');
	Result(prefix+'crg', 'grey title="CRG"');
	Result(prefix+'csg', 'grey title="CSG"');


