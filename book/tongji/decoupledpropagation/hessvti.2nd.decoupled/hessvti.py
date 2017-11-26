from rsf.proj import *

# =============================================================
# Fetch will not work unless you do the following:
#1. Download Hess VTI models from 'ftp://software.seg.org'
#2. Uncompress the package in ./hessvti

# Depth grid: 1501 
# Horizontal grid: 3617
# Vertical and horizontal spacing are both 20ft.
ft2km = 0.0003048
dx = 20*ft2km
dz = 20*ft2km
ifinite = 0

for m in ['vp','delta','epsilon']:
    sgy = 'hessvti/timodel_%s.segy' % m
    if m=='vp':
        Flow(m+'_hess',sgy,
             '''
             segyread read=data | 
             scale dscale=%g | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance label=Velocity unit=km/s
             ''' % (ft2km,dx,dz))
    else:
        Flow(m+'_hess',sgy,
             '''
             segyread read=data | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance
             ''' % (dx,dz))

# set finite Vs0 (samll nonzero) for stability and suppression S-wave
# ga = vp0^2/vs0^2*(epsilon-delta)
# The smaller ga, the less triplication in the SV wavefront
# The vairation of Vs0 is slower, the less reflected SV-wave
# Fletcher et al.,(2009)
if ifinite==1:
     Flow('vs_hess','vp_hess epsilon_hess delta_hess',
     '''
     math p=${SOURCES[0]} e=${SOURCES[1]} d=${SOURCES[2]} output="sqrt(p*p*abs(e-d)/0.5)"
     ''')
else:
     Flow('vs_hess','vp_hess',
     '''
     math output="input*0.5"
     ''')

name0='''
vp_hess
'''

name00='''
epsilon_hess delta_hess
'''

for ff in Split(name0):
        Result(ff,
        '''
        grey color=j scalebar=y bias=1.5 allpos=n barreverse=y wanttitle=n screenht=5 screenwd=14
        ''')

for gg in Split(name00):
        Result(gg,
        '''
        grey color=j scalebar=y allpos=n barreverse=y wanttitle=n screenht=5 screenwd=14
        ''')

#xv0=1200*dx
xv0=9.0
#xv0=0
zv0=0

par = dict(
    nxv=1400, ox=xv0, lx='x',  ux='km',
    nzv=1320, oz=zv0, lz='z',  uz='km',
    )
#    nxv=1801, ox=xv0, lx='x',  ux='km',
#    nxv=3600, ox=xv0, lx='x',  ux='km',

def get_model(name):    
    for bb in Split(name):
        Flow(bb,bb+'_hess',
        '''
        put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance unit=km/s |
        window n1=%d min1=%g n2=%d min2=%g
        ''' % (dx, dz, 
               par['nzv'],par['oz'],
               par['nxv'],par['ox']))
