from rsf.proj import *

#Depth: 1501, horizontal grid: 3617
#The spacing is: Depth 20ft, Horizontal: 20 ft

ft2km = 0.0003048
d = 20*ft2km

for m in ['vp','delta','epsilon','crho']:
    sgy = 'timodel_%s.segy' % m
    sgyz = sgy + '.gz'
    Fetch(sgyz,dir='Hess_VTI',server='ftp://software.seg.org',top='pub/datasets/2D')
    Flow(sgy,sgyz,'zcat $SOURCE',stdin=0)
    if m=='vp':
        Flow(m+'_hess',sgy,
             '''
             segyread read=data | 
             scale dscale=%g | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance label=Velocity unit=km/s
             ''' % (ft2km,d,d))
    else:
        Flow(m+'_hess',sgy,
             '''
             segyread read=data | 
             put d1=%g d2=%g unit1=km label1=Depth unit2=km label2=Distance
             ''' % (d,d))

Flow('vx_hess','vp_hess epsilon_hess',
     '''
     math e=${SOURCES[1]} output="input*sqrt(1+2*e)"
     ''')
Flow('eta_hess','epsilon_hess delta_hess',
     '''
     math e=${SOURCES[0]} d=${SOURCES[1]} output="((1+2*e)/(1+2*d)-1)/2"
     ''')

def get_model(name):    
    for par in Split(name):
        Flow(par,par+'_hess','cp')

sgy = ('timodel_shot_data_II_shot001-320.segy',
       'timodel_shot_data_II_shot321-720.segy')

for s in range(2):
    sgyz = sgy[s]+'.gz'
    Fetch(sgyz,dir='Hess_VTI',server='ftp://software.seg.org',top='pub/datasets/2D')
    Flow(sgy[s],sgyz,'zcat $SOURCE',stdin=0)

    shot = 'shot%d' % s
    Flow([shot,'t'+shot],sgy[s],
         '''
         segyread tfile=${TARGETS[1]}
         ''')

Flow('shot','shot0 shot1','cat axis=2 ${SOURCES[1]}')
Flow('tshot','tshot0 tshot1','cat axis=2 ${SOURCES[1]}')

# shot number
Flow('sn','tshot','dd type=float | headermath output=sy/100 | dd type=int')

# offset number
Flow('on','tshot','dd type=float | headermath output=offset/40 | dd type=int')
Flow('zomask','on','mask min=0 max=0')

Flow('head','on sn','cat axis=1 ${SOURCES[1]}')

def get_shots(shots):
    Flow(shots,'shot head',
         '''
         intbin head=${SOURCES[1]} xkey=0 ykey=1 |
         put
         o2=0 d2=%g label2=Offset unit2=km 
         o3=0 d3=%g label3=Shot   unit3=km
         ''' % (40*ft2km,100*ft2km))

def get_zodata(zo):
    Flow(zo,'shot zomask','headerwindow mask=${SOURCES[1]}')

