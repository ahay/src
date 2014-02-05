from rsf.proj import *
import rsf.gallery

method = rsf.gallery.method()

ft2km = 0.3048

def getvel(vel,veltype):
    segy = dict(migvel='sigsbee2a_migvel.sgy',
                strvel='sigsbee2a_stratigraphy.sgy')
    
    o2 = dict(migvel=10.025,
              strvel=10.000)
    d2 = dict(migvel=0.0375,
              strvel=0.0250)
    d1 = 0.0250

    Fetch(segy[veltype],'sigsbee')
    Flow(vel,segy[veltype],
         '''
         segyread read=data |
         scale rscale=0.001 |
         scale rscale=%g |
         put label=Velocity unit=km/s
         o1=0  d1=%g label1=Depth    unit1=km
         o2=%g d2=%g label2=Distance unit2=km
         ''' % (ft2km,d1*ft2km,o2[veltype]*ft2km,d2[veltype]*ft2km))

segy = 'sigsbee2a_nfs.sgy'
Fetch(segy,'sigsbee')

Fetch('sigexp.rsf','sigsbee')
    
Flow('data tdata',segy,'segyread tfile=${TARGETS[1]}')

Flow('zomask','tdata','window n1=1 f1=11 | mask max=0')

s0 = 10.925*ft2km
ds = 0.15*ft2km
dh = 0.075*ft2km

def getzo(zodata):
    Flow(zodata,'data zomask',
         '''
         headerwindow mask=${SOURCES[1]} | cut max1=1 |
         reverse which=2 | costaper nw2=10 |
         put o2=%g d2=%g label2=Distance unit2=km
         ''' % (s0,ds))

Flow('offset','tdata','window n1=1 f1=11 | dd type=float | math output=input/75 | dd type=int')
Flow('head','tdata offset','window n1=1 f1=2 | cat axis=2 ${SOURCES[1]} | transp')

def gethrzo(zodata):
    Flow(zodata,'sigexp.rsf','dd form=native')

def zoimage(image):
    Result(image,'grey title="Zero-Offset %s" label1=Depth label2=Distance' % method)

def getshots(shots):
    Flow(shots,'data head',
         '''
         intbin head=${SOURCES[1]} xkey=1 ykey=0 |
         put label1=Time unit1=s
         d2=%g o2=0  label2=Offset unit2=km
         d3=%g o3=%g label3=Shot   unit3=km |
         mutter half=n tp=0.3 v0=1.5
         ''' % (dh,ds,s0))

def psimage(image):
    Result(image,'grey title="Prestack %s" label1=Depth label2=Distance' % method)
