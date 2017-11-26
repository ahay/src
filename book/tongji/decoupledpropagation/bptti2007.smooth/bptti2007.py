from rsf.proj import *

# Fetch will not work unless you do the following:
# 1. Download ModelParams.tar.gz from http://www.freeusp.org/2007_BP_Ani_Vel_Benchmark/
# 2. Put in under $DATAPATH/BP
tgz = 'ModelParams.tar.gz'

Fetch(tgz,'BP',top=os.environ.get('DATAPATH'),server='local')

pars = Split('Epsilon Delta Vp Theta')

sgy = {}
for par in pars:
    sgy[par] = os.path.join('ModelParams',par + '_Model.sgy')

zcat = WhereIs('gzcat') or WhereIs('zcat')

Flow(sgy.values(),tgz,
     zcat + ' $SOURCE | tar -xvf -',stdin=0,stdout=-1)

# Vertical and horizontal spacing are both 20ft.
dx = 6.25/1000
dz = 6.25/1000
ifinite = 0

for m in pars:
    if m=='Vp':
        Flow(m+'_Model',sgy[m],
             '''
             segyread read=data | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance label=Velocity unit=km/s
             ''' % (dx,dz))
    else:
        Flow(m+'_Model',sgy[m],
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
     Flow('Vs_Model','Vp_Model Epsilon_Model Delta_Model',
     '''
     math p=${SOURCES[0]} e=${SOURCES[1]} d=${SOURCES[2]} output="sqrt(p*p*abs(e-d)/0.5)"
     ''')
else:
     Flow('Vs_Model','Vp_Model',
     '''
     math output="input*0.5"
     ''')

name0='''
Vp_Model
'''

name00='''
Epsilon_Model Delta_Model Theta_Model
'''

for ff in Split(name0):
        Result(ff,
        '''
        grey color=j scalebar=y bias=1.5 allpos=n barreverse=y wanttitle=n screenht=5 screenwd=30
        ''')

for gg in Split(name00):
        Result(gg,
        '''
        grey color=j scalebar=y allpos=n barreverse=y wanttitle=n screenht=5 screenwd=30
        ''')

xv0=4800.0*dx
zv0=200*dx

par = dict(
    nxv=1600, ox=xv0, lx='x',  ux='km',
    nzv=1600, oz=zv0, lz='z',  uz='km',
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
