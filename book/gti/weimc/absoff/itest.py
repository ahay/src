from rsfproj import *
import sys
from math import *
import spmig

def model(par):
    Flow('pvel',None,
         '''
         spike nsp=1 mag=2000
         n1=%(nz)d o1=%(oz)g d1=%(dz)g
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=z label2=x
         ''' % par)
    Flow('cvel','pvel','math output=1000')
    Flow('pslo','pvel','math output=1/input | transp | spray axis=2 n=1')
    Flow('cslo','cvel','math output=1/input | transp | spray axis=2 n=1')

    Flow('wave',None,
         '''
         spike nsp=1 mag=1 k1=1
         n1=%(nt)d d1=%(dt)g o1=0|
         put label1=t label2=x label3=y 
         ''' % par )    
    
def data(case,dat,DIP,ANG,par):

    dip = 'dip' + DIP + ANG
    Flow(dip,None,'spike n1=%d o1=%g d1=%g mag=%g'
         % (par['nx'],par['ox'],par['dx'],
            tan(pi*int(DIP)/180.)))

    ref = 'ref' + DIP + ANG
    Flow(ref,dip,'math output="%g+x1*input"'
         % (par['zcig']-par['xcig']*tan(pi*int(DIP)/180.) ))
    
    if case == 'p':
        par['vel2']=2000
    else:
        par['vel2']=1000

    par['os'] = par['xcig'] - \
                par['zcig'] * tan(pi*(int(ANG)-int(DIP))/180.)
    print DIP, ANG, par['os']

    Flow(dat,[ref,dip],
         '''
         kirmod vel=2000 vel2=%(vel2)g dip=${SOURCES[1]}
         nt=%(nt)d  dt=%(dt)g freq=15         
         nh=%(nh)d  h0=%(oh)g dh=%(dh)g
         ns=%(ns)d  s0=%(os)g ds=%(ds)g |
         put label1=t label2=h
         ''' % par)

    # wavefields
    sou = 'sou' + dat
    rec = 'rec' + dat
    spmig.wflds(sou,rec,'wave',dat,par)

def migrate(case,imco,dat,img,cig,par):

    locpar = par
    if(imco=='o'): locpar['misc']='itype=o'
    if(imco=='t'): locpar['misc']='itype=t nht=160 oht=-0.200 dht=0.0025             jcx=500'
    if(imco=='x'): locpar['misc']='itype=x hsym=y nhx=50                             jcx=500'
    if(imco=='z'): locpar['misc']='itype=x hsym=y        nhz=50                      jcx=500'
    if(imco=='m'): locpar['misc']='itype=x hsym=y nhx=50 nhz=50                      jcx=500'
    if(imco=='h'): locpar['misc']='itype=h        nhh=50 dhh=10 nha=180 dha=2 oha=0  jcx=500'

    sou = 'sou' + dat
    rec = 'rec' + dat
    # migration
    if case == 'p':
        spmig.imagePW(img,cig,'pslo',       sou,rec,locpar)
    else:
        spmig.imageCW(img,cig,'pslo','cslo',sou,rec,locpar)
        
    

    
