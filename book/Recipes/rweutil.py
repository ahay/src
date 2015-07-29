from rsf.proj import *

# ------------------------------------------------------------

def C2Rdat(mapRC, mapCC,cos,par):
    Flow(  mapRC,[mapCC,cos],
         '''
         window |
         spray axis=2 n=1 o=%(oz)g d=%(dz)g |
         pad beg2=10 n2out=20 |
         c2r adj=n linear=n nsz=3 nsx=3
         rays=${SOURCES[1]} |
         put label1=g label2=t
         ''' % par)

def R2Cdat(mapCC, mapRC,cos,par):    
    par['nzcut']=20
    par['ozcut']=-10*par['dz']
    Flow(  mapCC,[mapRC,cos],
         '''
         c2r adj=y linear=n nsz=3 nsx=3
         rays=${SOURCES[1]}
         a1n=%(nx)d    a1o=%(ox)g    a1d=%(dx)g
         a2n=%(nzcut)d a2o=%(ozcut)g a2d=%(dz)g |
         window n2=1 min2=%(oz)g
         ''' % par)

# ------------------------------------------------------------

def C2Rimg(mapRC, mapCC,cos,par):
    Flow(  mapRC,[mapCC,cos],
         '''
         window |
         c2r adj=n linear=y
         rays=${SOURCES[1]} |
         put label1=g label2=t
         ''' % par)

def R2Cimg(mapCC, mapRC,cos,par):
    Flow(  mapCC,[mapRC,cos],
         '''
         window |
         c2r adj=y linear=n nsz=3 nsx=3
         rays=${SOURCES[1]}
         a1n=%(nx)d a1o=%(ox)g a1d=%(dx)g
         a2n=%(nz)d a2o=%(oz)g a2d=%(dz)g |
         put label1=z label2=x
         ''' % par)

# ------------------------------------------------------------

def SPmod(datRC, wavRC,imgRC,abmRC,abrRC,par):
    Flow( datRC,[wavRC,imgRC,abmRC,abrRC],
         '''
         rwesrmig ntap=25 adj=y verb=y method=1
         nw=%(nw)d ow=%(ow)g dw=%(dw)g 
         img=${SOURCES[1]}
         abm=${SOURCES[2]}
         abr=${SOURCES[3]} |
         put label1=g label2=t label3=w
         ''' % par)

def SPmig(imgRC, wavRC,datRC,abmRC,abrRC,par):
    Flow( imgRC,[wavRC,datRC,abmRC,abrRC],
         '''
         rwesrmig ntap=25 adj=n verb=y method=1
         rwf=${SOURCES[1]}
         abm=${SOURCES[2]}
         abr=${SOURCES[3]} |
         put label1=g label2=t
         ''')

# ------------------------------------------------------------

def ZOmod(datRC, imgRC,abmRC,abrRC,par):
    Flow( datRC,[imgRC,abmRC,abrRC],
         '''
         rwezomig ntap=25 adj=y verb=y method=1
         nw=%(nw)d ow=%(ow)g dw=%(dw)g 
         abm=${SOURCES[1]}
         abr=${SOURCES[2]} |
         put label1=g label2=t label3=w
         ''' % par)

def ZOmig(migRC, datRC,abmRC,abrRC,par):
    Flow( migRC,[datRC,abmRC,abrRC],
         '''
         rwezomig ntap=25 adj=n verb=y method=1
         abm=${SOURCES[1]}
         abr=${SOURCES[2]} |
         put label1=g label2=t
         ''')

# ------------------------------------------------------------

def t2w (wdat,tdat,par):
    Flow(wdat,tdat,
         '''
         fft1 opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g |
         put label1=w label2=x label3=y |
         transp plane=12 |
         transp plane=23 
         ''' % par)


def w2t (tdat,wdat,par):
    Flow(tdat,wdat,
         '''
         transp |
         pad beg1=%d n1out=%d |
         fft1 opt=n inv=y |
         window f1=%d |
         pad n1out=%d |
         put label1=t label2=x o1=0
         ''' % (int(par['ow']/par['dw']),
                par['nT']/2+1,
                par['kT'],
                par['nT']))
    
