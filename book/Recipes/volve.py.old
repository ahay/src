try:    from rsf.cluster import *
except: from rsf.proj    import *
import fdmod

# ------------------------------------------------------------
def param():
	par = dict(
		nx=247,  dx=0.050, ox=0, lx='x', ux='km',
		ny=137,  dy=0.050, oy=0, ly='y', uy='km',
		nz=451,  dz=0.010, oz=0, lz='z', uz='km'
	)
	fdmod.param(par)
	dx=par['xmax']-par['xmin']
	dy=par['ymax']-par['ymin']
	par['xyratio']=1.0*dy/dx

	return par

# ------------------------------------------------------------
def pabin():
	par = dict(
		nx=493,  dx=0.025, ox=0, lx='x', ux='km',
		ny=273,  dy=0.025, oy=0, ly='y', uy='km',
		nz=451,  dz=0.010, oz=0, lz='z', uz='km',
		nt=1876, dt=0.004, ot=0, lt='t', ut='s'
	)
	fdmod.param(par)
	dx=par['xmax']-par['xmin']
	dy=par['ymax']-par['ymin']
	par['xyratio']=1.0*dy/dx

	return par
	
# ------------------------------------------------------------
def xygraph(custom,par):
    return '''
           graph title=""
           yreverse=n transp=y 
	   screenratio=%g screenht=8
           min2=%g max2=%g label2=%s unit2=%s
           min1=%g max1=%g label1=%s unit1=%s
           %s %s
            '''%(par['xyratio'],
                 par['xmin'],par['xmax'],par['lx'],par['ux'],
                 par['ymin'],par['ymax'],par['ly'],par['uy'],
                 par['labelattr'],custom)

# ------------------------------------------------------------
def xygrey(custom,par):
    return '''
	   transp | reverse which=1 |   
	   grey pclip=100 title="" wantaxis=n 
	   screenratio=%g screenht=8
           min2=%g max2=%g label2=%s unit2=%s
           min1=%g max1=%g label1=%s unit1=%s
           %s %s 
	   '''%(par['xyratio'],
                 par['xmin'],par['xmax'],par['lx'],par['ux'],
                 par['ymin'],par['ymax'],par['ly'],par['uy'],
                 par['labelattr'],custom)

# ------------------------------------------------------------
def model(par):

    datadir = 'data/volve/'
    MODPAR = ['vp','vs','eps','dlt']

    for i in MODPAR:
        Flow(i,None,
             'segyread tape=%s |'%(datadir+i+'.segy') +
             '''
             put label1=%(lz)s unit1=%(uz)s
             n2=%(nx)g o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s
             n3=%(ny)g o3=%(oy)g d3=%(dy)g label3=%(ly)s unit3=%(uy)s
             '''%par
             )
	Result(i,
	       'byte gainpanel=a pclip=100 mean=y |' +
	       fdmod.cgrey3d('color=j frame1=0 movie=1 dframe=10',par))
        Plot(i,'window n1=1 f1=225 |' + xygrey('mean=y',par))

# ------------------------------------------------------------
def winCRGS():

    nCRGS=2880
    allCRGS=range(1,1+nCRGS)
    EXCLUDE=[21,22,42,146,218,265,293,305,422,437,508,509,698,751,773,785,902,993,994,1059,1106,1178,1253,1265,1370,1382,1586,1658,1733,1745,1850,1862,2138,2213,2225,2342,2546,2618,2693,2705,2822]
    
    winCRGS=[x for x in allCRGS if x not in EXCLUDE]

    return winCRGS

# ------------------------------------------------------------
def tstCRGS():
    return [1,120,240,241,360,480,481,600,720,721,840,960,961,1080,1200,1201,1320,1440,1441,1560,1680,1681,1800,1920,1921,2040,2160,2161,2280,2400,2401,2520,2640,2641,2760,2880]
 
# ------------------------------------------------------------
def rrCRGS(iCRG):

	datadir='data/volve/crgs/'
	
	CRGS=winCRGS()
	index=CRGS.index(iCRG)

	rr3xCRGS=[]
	for line in open(datadir+'x.asc', "r" ).readlines():
		for value in line.split( ' ' ):
			rr3xCRGS.append( value )
	rr3yCRGS=[]
	for line in open(datadir+'y.asc', "r" ).readlines():
		for value in line.split( ' ' ):
			rr3yCRGS.append( value )
	rr3zCRGS=[]
	for line in open(datadir+'z.asc', "r" ).readlines():
		for value in line.split( ' ' ):
			rr3zCRGS.append( value )

	rr3x=float(rr3xCRGS[index])
	rr3y=float(rr3yCRGS[index])
	rr3z=float(rr3zCRGS[index])

	return [rr3x,rr3y,rr3z]

# ------------------------------------------------------------
def binning(bin,trc,ss3,rr3,pab):

	datadir='data/volve/crgs/'

	Flow(rr3+'x',None,'window n1=1 f1=0 <%s'%(datadir+rr3+'.hh'))
	Flow(rr3+'y',None,'window n1=1 f1=1 <%s'%(datadir+rr3+'.hh'))	
	Flow(ss3+'x',None,'window n1=1 f1=0 <%s'%(datadir+ss3+'.hh'))
	Flow(ss3+'y',None,'window n1=1 f1=1 <%s'%(datadir+ss3+'.hh'))
	
	# Receiver indices in the binning grid
	for c in ['x','y']:   
	     Flow(rr3+c+'I',
       	          rr3+c,
	          'math output="(input-%(xmin)g)/%(dx)g" | dd type=int'%pab)
	     Flow(ss3+c+'I',
 	          ss3+c,
		  'math output="(input-%(xmin)g)/%(dx)g" | dd type=int'%pab)         
	
        #Flow(trc,None,'window <%s'%(datadir+trc+'.hh'))

	# key indices
	Flow(trc+'key',[rr3+'xI',rr3+'yI',ss3+'xI',ss3+'yI'],
	     '''
	     cat axis=2 space=n ${SOURCES[1:4]} |
	     transp
	     ''')

	# binning
	Flow(bin+'tmp',[trc,trc+'key'],
	     '''
	     srbin3d verb=y key=${SOURCES[1]}
	     n1=%(nx)d o1=%(ox)g d1=%(dx)g
	     n2=%(ny)d o2=%(oy)g d2=%(dy)g
	     on1=%(nx)d oo1=%(ox)g od1=%(dx)g
	     on2=%(ny)d oo2=%(oy)g od2=%(dy)g
	     '''%pab)

	Flow(bin,bin+'tmp',
	     '''
	     put o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s
                 o3=%(oy)g d3=%(dy)g label3=%(ly)s unit3=%(uy)s
	     '''%pab)
 
