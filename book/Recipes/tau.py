# tau domain zero-offset migration

from rsf.proj import *

def init(par):
	global nx,x0,dx,x, \
		   nz,z0,dz,z, \
		   nT,T0,dT,T, \
		   nt,t0,dt,t

	if par.has_key('nx'):
		nx = par['nx']; x0 = par['x0']; dx = par['dx']
		x = map(lambda t : x0 + t * dx,range(nx)); xlen = x[-1] - x[0]
	if par.has_key('nz'):
		nz = par['nz']; z0 = par['z0']; dz = par['dz']
		z = map(lambda t : z0 + t * dz,range(nz)); zlen = z[-1] - z[0]
	if par.has_key('nT'):
		nT = par['nT']; T0 = par['T0']; dT = par['dT']
		T = map(lambda t : T0 + t * dT,range(nT)); Tlen = T[-1] - T[0]
	if par.has_key('nt'):
		nt = par['nt']; t0 = par['t0']; dt = par['dt']
		t = map(lambda t : t0 + t * dt,range(nt)); tlen = t[-1] - t[0]

def compute_tau(tau1,vmap):
	Flow(tau1,vmap,'math output="1./input" | integral1 rule=s')

def compute_sigma(sigma,tau1):
	Flow(sigma,tau1,'transp | deriv scale=y | transp')

def coord_xmid(out):
	'''returns surface midpoint coordinate'''
	Flow([out,out+'.asc'],None,'''
	echo %g %g > ${TARGETS[1]} &&
	echo n1=2 in=${TARGETS[1]} data_format=ascii_float |
	dd form=native
	''' % (0.,x[len(x) / 2]))

def coord_center(out,tau=False):
	'''returns center point coordinate'''
	if tau:
		depth = T[len(T) / 2]
	else:
		depth = z[len(z) / 2]

	Flow([out,out+'.asc'],None,'''
	echo %g %g > ${TARGETS[1]} &&
	echo n1=2 in=${TARGETS[1]} data_format=ascii_float |
	dd form=native
	''' % (depth,x[len(x) / 2]))

def mapping(out,inp,tau1,inv=False):
	'''z to tau : inv=False
	tau to z : inv=True'''
	if inv:
		Flow(out,[inp,tau1],'''
		pseudodepth tau=${SOURCES[1]}
		inv=y linear=y | put label1=Depth unit1=m
		''')
	else:
		Flow(out,[inp,tau1],'''
		pseudodepth tau=${SOURCES[1]}
		inv=n n=%d o=0 d=%g linear=y | put label1=Tau unit1=s
		''' % (nT,dT))

def rtm_vti(imag,wave,data,vver,vnmo,heta,vmap,sigm,cr,n3=1,b=[],c=[],tau=False):
	'''VTI zero-offset migration'''
	if not len(b):
		b = [0,0,0,0]
	if not len(c):
		c = [10,10,10,10]

	if tau: # nonorthogonal tau
		Flow([imag,wave],[data,vver,vnmo,heta,vmap,sigm,cr],'''
		zomvti inv=n opt=y verb=y tau=y
		vz=${SOURCES[1]} vnmo=${SOURCES[2]}
		eta=${SOURCES[3]} vmap=${SOURCES[4]}
		sigm=${SOURCES[5]} cr=${SOURCES[6]}
		wave=${TARGETS[1]}
		bzl=%d bzh=%d bxl=%d bxh=%d
		czl=%g czh=%g cxl=%g cxh=%g n3=%d |
        put label1=Tau unit1=s label2=Distance unit2=m
		''' % (b[0],b[1],b[2],b[3],
			   c[0],c[1],c[2],c[3],n3))
	else:   # cartesian
		Flow([imag,wave],[data,vver,vnmo,heta,cr],'''
		zomvti inv=n opt=y verb=y tau=n
		vz=${SOURCES[1]} vnmo=${SOURCES[2]}
		eta=${SOURCES[3]} cr=${SOURCES[4]}
		wave=${TARGETS[1]}
		bzl=%d bzh=%d bxl=%d bxh=%d
		czl=%g czh=%g cxl=%g cxh=%g n3=%d |
        put label1=Depth unit1=m label2=Distance unit2=m
		''' % (b[0],b[1],b[2],b[3],
			   c[0],c[1],c[2],c[3],n3))

def rtm_iso(imag,wave,data,velo,vmap,sigm,cr,n3=1,b=[],c=[],tau=False):
	'''isotropic zero-offset migration'''
	if not len(b):
		b = [0,0,0,0]
	if not len(c):
		c = [10,10,10,10]

	if tau: # nonorthogonal tau
		Flow([imag,wave],[data,velo,vmap,sigm,cr],'''
		zomiso inv=n opt=y verb=y tau=y
		velo=${SOURCES[1]}
		vmap=${SOURCES[2]}
		sigm=${SOURCES[3]}
		cr=${SOURCES[4]}
		wave=${TARGETS[1]}
		bzl=%d bzh=%d bxl=%d bxh=%d
		czl=%g czh=%g cxl=%g cxh=%g n3=%d |
        put label1=Tau unit1=s label2=Distance unit2=m
		''' % (b[0],b[1],b[2],b[3],
			   c[0],c[1],c[2],c[3],n3))
	else:   # cartesian
		Flow([imag,wave],[data,velo,cr],'''
		zomiso inv=n opt=y verb=y tau=n
		velo=${SOURCES[1]}
		cr=${SOURCES[2]}
		wave=${TARGETS[1]}
		bzl=%d bzh=%d bxl=%d bxh=%d
		czl=%g czh=%g cxl=%g cxh=%g n3=%d |
        put label1=Depth unit1=m label2=Distance unit2=m
		''' % (b[0],b[1],b[2],b[3],
			   c[0],c[1],c[2],c[3],n3))

def plot_sour(sour):
	'''source wavelet'''
	Result(sour,'transp | graph plotfat=10 title="total %g s"' % t[-1])

def compute_vnmo(vnmo,vver,delta):
	Flow(vnmo,[vver,delta],'math d=${SOURCES[1]} output="input*sqrt(1.+2.*d)"')

def compute_eta(eta,epsilon,delta):
	Flow(eta,[epsilon,delta],'math d=${SOURCES[1]} output="(input-d)/(1.+2.*d)"') # clip2 lower=1.e-4 # stablize stiffness tensor













# def init(mod=0):
# 	'''model=0 : isotropic linear velocity gradient + lens
#        model=1 : isotropic marmousi
#        model=2 : VTI hess'''
# 	global nx,x0,dx,nz,z0,dz,nT,T0,dT,nt,t0,dt,j3,n3,d3,x,z,T,t,xlen,zlen,tlen,Tlen,b,c,freq,fsnap,m2km,M2km,model

# 	if 0 == mod:
# 		nx = 501;  x0 = 0.; dx = 5.
# 		nz = 401;  z0 = 0.; dz = 5.
# 		nT = 251;  T0 = 0.; dT = .0042
# 		nt = 2001; t0 = 0.; dt = .0008
# 		bzl = 40;    bzh = 20;    bxl = 20;    bxh = 20
# 		czl = 4.e-5; czh = 1.e-4; cxl = 1.e-4; cxh = 1.e-4
# 		j3 = 10
# 		fsnap = 180
# 		freq = 20
# 	elif 1 == mod:
# 		nx = 751;  x0 = 0.; dx = 4.
# 		nz = 751;  z0 = 0.; dz = 4.
# 		nT = 501;  T0 = 0.; dT = .0024
# 		nt = 12001;t0 = 0.; dt = .00016
# 		bzl = 40;    bzh = 20;    bxl = 20;    bxh = 20
# 		czl = 4.e-5; czh = 1.e-4; cxl = 1.e-4; cxh = 1.e-4
# 		j3 = 50
# 		fsnap = 240
# 		freq = 25
# 	elif 2 == mod:
# 		ft = .3048
# 		nx = 751;  x0 = 0.; dx = 20. * ft
# 		nz = 601;  z0 = 0.; dz = 20. * ft
# 		nT = 501;  T0 = 0.; dT = .00324
# 		nt = 5001; t0 = 0.; dt = .0005
# 		bzl = 40;    bzh = 20;    bxl = 20;    bxh = 20
# 		czl = 4.e-5; czh = 1.e-4; cxl = 1.e-4; cxh = 1.e-4
# 		j3 = 50
# 		fsnap = 100
# 		freq = 25
# 	elif 3 == mod:
# 		nx = 501; x0 = 0.; dx = 5.
# 		nz = 401; z0 = 0.; dz = 5
# 		nT = 301; T0 = 0.; dT = .004
# 		nt = 2001;t0 = 0.; dt = .0005
# 		bzl = 40;    bzh = 20;    bxl = 20;    bxh = 20
# 		czl = 4.e-5; czh = 1.e-4; cxl = 1.e-4; cxh = 1.e-4
# 		j3 = 25
# 		fsnap = 80
# 		freq = 25
# 	else:
# 		pass

# 	b = [bzl,bzh,bxl,bxh]
# 	c = [czl,czh,cxl,cxh]

# 	n3 = 1 + int((nt-1)/j3)
# 	d3 = dt*j3
	
# 	model = mod

# 	m2km = 'put d2=%g unit2=km d1=%g unit1=km' % (.001*dx,.001*dz)
# 	M2km = 'put d2=%g unit2=km' % (.001*dx)

# def get_velo(velo,vmap):
# 	if 0 == model:
# 		Flow(velo,None,'''
# 		math n1=%d o1=0 d1=%g 
# 			 n2=%d o2=0 d2=%g
# 		output="2100+.025*x1+.14*x2-750*exp(-4.e-6*((x1-%g)^2+(x2-%g)^2))" |
# 		put label1=z unit1=m label2=x unit2=m
# 		''' % (nz,dz,nx,dx,z0+.55*zlen,x0+.45*xlen))
# 		Flow(vmap,velo,'math output=input')
# 	elif 1 == model:
# 		Flow(velo,'marmvel.hh','''
# 		dd form=native | window n2=%d f2=1100 |
# 		put o1=0 d1=%g o2=0 d2=%g label1=z unit1=m label2=x unit2=m
# 		''' % (nx,dz,dx))
# 		# smoothed marmousi as mapping velocity
# 		Flow(vmap,velo,'smooth rect1=8 rect2=8 repeat=2')
# 	elif 2 == model:
# 		window  = '''window n1=%d f1=0 n2=%d f2=0 |
# 		put o1=%g d1=%g label1=z unit1=m
# 			o2=%g d2=%g label2=x unit2=m
# 		''' % (nz,nx,z0,dz,x0,dx)

# 		# return vertical velocity only
# 		Flow(velo,'timodel_vp.segy','''
# 		segyread tfile=/dev/null verb=y | scale rscale=.3048 |
# 		''' + window)
# 		# smoothed vertical velocity as mapping velocity
# 		Flow(vmap,velo,'smooth rect1=8 rect2=8 repeat=2')
# 	elif 3 == model:
# 		slow0 = 1./1500
# 		dslowx = -3e-8
# 		dslowz = -8e-8
# 		Flow(velo,None,'''
# 		math n1=%d o1=0 d1=%g
# 			 n2=%d o2=0 d2=%g output="%g+%g*x1+%g*x2" |
# 		math output="1./input" |
# 		put label1=z unit1=m label2=x unit2=m
# 		''' % (nz,dz,nx,dx,slow0,dslowz,dslowx))
# 		# exact velocity as vmap
# 		Flow(vmap,velo,'math output=input')
# 		# tau (analytical)
# 		Flow('tauC',velo,'math output="(%g+%g*x2)*x1+%g*x1^2"' % (slow0,dslowx,dslowz*.5))
# 	else:
# 		pass

# def get_velo2(vnmo,heta,vver):
# 	'''anisotropic velocities'''
# 	window  = '''window n1=%d f1=0 n2=%d f2=0 |
# 	put o1=%g d1=%g label1=z unit1=m
# 	    o2=%g d2=%g label2=x unit2=m
# 	''' % (nz,nx,z0,dz,x0,dx)

# 	Flow('epsi','timodel_epsilon.segy','''
# 	segyread tfile=/dev/null verb=y |
# 	''' + window)
# 	Flow('delt','timodel_delta.segy','''
# 	segyread tfile=/dev/null verb=y |
# 	''' + window)
# 	Flow(vnmo,[vver,'delt'],'math d=${SOURCES[1]} output="input*sqrt(1.+2.*d)"')
# 	Flow(heta,'epsi delt','math d=${SOURCES[1]} output="(input-d)/(1.+2.*d)"')

# def get_spik(spik):
# 	'''synthetic zero-offset section with three spikes'''
# 	spk0 = 1 + int(1./(freq*dt))
# 	dspk = int(.5/ dt)
# 	spk = range(spk0,spk0+3*dspk,dspk)

# 	Flow(spik,None,'''
# 	spike nsp=3 mag=1 n1=%d o1=0 d1=%g k1=%d,%d,%d |
# 	ricker1 frequency=%g | scale axis=1 | reverse which=1 opt=i |
# 	transp | put label1=x unit1=m
# 	''' % (nt,dt,spk[0],spk[1],spk[2],freq)) 


# def plot_grid(grid,tau,veloC,veloT):
# 	'''gridc : tau grid in cartesian domain
# 	gridt : tau grid in tau domain'''

# 	if 0 == model:
# 		ct = 'contour c0=0 dc=.08 allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n screenwd=10 screenht=8'
# 		cx = 'contour c0=0 dc=%g allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n screenwd=10 screenht=8' % (25*dx)
# 	elif 1 == model:
# 		ct = 'contour c0=0 dc=.12 allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n screenwd=10 screenht=8'
# 		cx = 'contour c0=0 dc=%g allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n screenwd=10 screenht=8' % (45*dx)
# 	elif 2 == model:
# 		ct = 'contour c0=0 dc=.12 allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n wanttitle=n screenwd=10 screenht=8'
# 		cx = 'contour c0=0 dc=%g allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n wanttitle=n screenwd=10 screenht=8' % (45*dx)
# 	elif 3 == model:
# 		ct = 'contour c0=0 dc=.05 allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n wanttitle=n screenwd=10 screenht=8'
# 		cx = 'contour c0=0 dc=%g allpos=y plotcol=0 plotfat=5 wantaxis=n wanttitle=n wanttitle=n screenwd=10 screenht=8' % (25*dx)
# 	else:
# 		pass

# 	Plot(grid + 'C_1',tau,ct)
# 	Plot(grid + 'C_2',tau,'math output=x2 | ' + cx)
# 	Plot(grid + 'C_3',veloC,m2km + ' | grey color=j allpos=y bias=1500 title= labelsz=8 labelfat=3 screenwd=10 screenht=8')
# 	Result(grid + 'C',[grid+'C_3',grid+'C_2',grid+'C_1'],'Overlay')

# 	Plot(grid + 'T_1',veloT,'math output=x1 | ' + ct)
# 	Plot(grid + 'T_2',veloT,'math output=x2 | ' + cx)
# 	Plot(grid + 'T_3',veloT,M2km + ' | grey color=j allpos=y bias=1500 title= labelsz=8 labelfat=3 screenwd=10 screenht=8')
# 	Result(grid + 'T',[grid+'T_3',grid+'T_2',grid+'T_1'],'Overlay')

# def plot_veloB(velo,veloB):
# 	'''placeholder only'''
# 	Result(velo + 'B',veloB,m2km + ' | grey color=j allpos=y bias=1500 title= labelsz=8 labelfat=3 screenwd=10 screenht=8')

# def plot_spik(spik,sour):
# 	'''sour : source wavelet'''
# 	pad1 = int((nx-1) / 2)
# 	pad2 = nx - 1 - pad1
# 	Result(sour,'transp | graph plotfat=10 title="total %g sec"' % t[-1])
# 	Result(spik,sour,'pad beg1=%d end1=%d | put o1=%g d1=%g | bandpass fhi=.1 | put o1=%g d1=%g unit1=km | grey title=' % (pad1,pad2,x0,dx,1e-3*x0,1e-3*dx) )

# def plot_snap(snap,tau,waveC,waveT):
# 	'''take snapshots of wavefield'''
# 	if 0 == model:
# 		f1z=49; f1t=42; f2=29
# 	elif 1 == model:
# 		f1z=54; f1t=43; f2=24
# 	elif 2 == model:
# 		f1z=45; f1t=43; f2=24
# 	elif 3 == model:
# 		f1z=25; f1t=25; f2=25 # untested
# 	else:
# 		pass
# 	Flow(snap + 'C',waveC,'window n3=1 f3=%d n1=%d f1=%d n2=%d f2=%d' % (fsnap,nz,f1z,nx,f2))
# 	Flow(snap + 'T',waveT,'window n3=1 f3=%d n1=%d f1=%d n2=%d f2=%d' % (fsnap,nT,f1t,nx,f2))
# 	mapping(snap + 'B',snap + 'T',tau,inv=True)

# 	wgrey = ' | grey grid=y title= labelsz=8 labelfat=3 screenwd=10 screenht=8'
# 	Result(snap + 'C',m2km + wgrey)
# 	Result(snap + 'T',M2km + wgrey)
# 	Result(snap + 'B',m2km + wgrey)
