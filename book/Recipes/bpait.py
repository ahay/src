## 
 # BPAIT parameters
 ##

from rsf.proj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def param():
    par = dict(
    nt=8001,  ot=0,      dt=0.00050, lt='t', ut='s',
    nx=5395,  ox=2.4384, dx=0.01250, lx='x', ux='km',
    nz=1911,  oz=0,      dz=0.00625, lz='z', uz='km'
    )

    return par

def fwipar(par):
    par['frq']=10
    par['kt']=100
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jsnap']=500
    par['jdata']=1
    par['wweight']=50
    par['wclip']=0.5

    par['nqz']=par['nz']/2
    par['oqz']=par['oz']
    par['dqz']=par['dz']*2
    
    par['nqx']=par['nx']/2
    par['oqx']=par['ox']
    par['dqx']=par['dx']*2

# ------------------------------------------------------------
def getmigvel(velo,par):

    if(local):
        migvelfile = 'DATA/bpait/bpaitvel.hh'
    else:
        migvelfile = 'bpaitvel.hh'
        Fetch(datafile,'bpait')
        
    Flow(velo,migvelfile,
         '''
         dd form=native |
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s label2=%(lx)s
         o2=%(ox)g d2=%(dx)g  unit1=%(uz)s  unit2=%(ux)s |
         scale rscale=0.001
         ''' % par)
    
