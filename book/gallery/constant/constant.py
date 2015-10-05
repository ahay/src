from rsf.proj import *
import rsf.gallery 
import math

method = rsf.gallery.method()

z = (1,1.5,2)
a = (20,0,0)
x = (0.5,1)
v = (1,1.1,0.9)

n1 = 7501

for ref in range(3):
    a0 = a[ref]*math.pi/180

    rdip = 'rdip%d' % ref
    Flow(rdip,None,
         '''
         spike n1=%d d1=0.001 o1=-3 label1=Distance mag=%g
         ''' % (n1,math.tan(a0)))
    
    refl = 'refl%d' % ref
    Flow(refl,rdip,
         'math output="%g+x1*input"' % z[ref])

    ampl = 'ampl%d' % ref
    if ref==1:
        Flow(ampl,rdip,'spike k1=3501 mag=200')
    else:
        Flow(ampl,rdip,'math output=1')

for case in ('rdip','refl','ampl'):
    Flow(case,[case+'0',case+'1',case+'2'],'cat axis=2 ${SOURCES[1:3]}')

def get_zodata(data):
    Flow(data,'refl rdip ampl',
         '''
         kirmod cmp=y dip=${SOURCES[1]} refl=${SOURCES[2]}
         nh=1  dh=0.1  h0=0 twod=y
         ns=351 ds=0.01 s0=-1
         freq=10 dt=0.004 nt=1500
         vel=1 verb=y | window | put label2=Distance unit2=km
         ''',split=[1,n1], reduce='add')

def get_cmps(data):
    Flow(data,'refl rdip ampl',
         '''
         kirmod cmp=y dip=${SOURCES[1]} refl=${SOURCES[2]}
         nh=201 dh=0.01  h0=0 twod=y
         ns=351 ds=0.01 s0=-1
         freq=10 dt=0.004 nt=1500
         vel=1 verb=y | put label2=Offset unit2=km label3=Midpoint unit3=km
         ''',split=[1,n1], reduce='add')

def get_impulse(imp,data):
    Flow(imp,data,'spike k1=501 k2=175 | ricker1 frequency=10')

def zo_image(image):
    Result(image,'window min2=0 max2=1.5 | grey title="Zero-Offset %s" ' % method)

def impulse_response(image):
    Result(image,'window n1=751 | grey title="%s Impulse Response" ' % method)

