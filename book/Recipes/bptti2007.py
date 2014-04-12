from rsf.proj import *

# =============================================================
# download BP 2007 TTI  models from 'ftp://software.seg.org'
dx = 6.25/1000
dz = 6.25/1000
ifinite = 0

for m in ['Vp','Delta','Epsilon','Theta']:
    sgy = '/data/cjb/bptti2007/%s_Model.sgy' % m
    if m=='Vp':
        Flow(m+'_Model',sgy,
             '''
             segyread read=data | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance label=Velocity unit=km/s
             ''' % (dx,dz))
    else:
        Flow(m+'_Model',sgy,
             '''
             segyread read=data | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance
             ''' % (dx,dz))

Flow('Vx_Model','Vp_Model Epsilon_Model',
     '''
     math e=${SOURCES[1]} output="input/sqrt(1+2*e)"
     ''')

Flow('Eta_Model','Epsilon_Model Delta_Model',
     '''
     math e=${SOURCES[0]} d=${SOURCES[1]} output="((1+2*e)/(1+2*d)-1)/2"
     ''')

# set finite Vs0 (samll nonzero) for stability and suppression S-wave
# ga = vp0^2/vs0^2*(epsilon-delta)
# The smaller ga, the less triplication in the SV wavefront
# The vairation of Vs0 is slower, the less reflected SV-wave
# Fletcher et al.,(2009)
if ifinite==1:
     Flow('Vs_Model','Vp_Model Epsilon_Model Delta_Model',
     '''
     math p=${SOURCES[0]} e=${SOURCES[1]} d=${SOURCES[2]} output="sqrt(p*p*abs(e-d)/0.5)"
     ''')
else:
     Flow('Vs_Model','Vp_Model',
     '''
     math output="input*0.6"
     ''')

xv0=7250*dx
zv0=0

par = dict(
    nxv=1500, ox=xv0, lx='x',  ux='km',
    nzv=1500, oz=zv0, lz='z',  uz='km',
#    nxv=12596, ox=xv0, lx='x',  ux='km',
#    nzv=1801, oz=zv0, lz='z',  uz='km',
    )
def get_model(name):    
    for bb in Split(name):
        Flow(bb,bb+'_Model',
        '''
        put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance unit=m/s |
        window n1=%d min1=%g n2=%d min2=%g
        ''' % (dx, dz, 
               par['nzv'],par['oz'],
               par['nxv'],par['ox']))
