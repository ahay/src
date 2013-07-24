from rsf.proj import *

velsegy = 'vel_z6.25m_x12.5m_exact.segy'

Fetch(velsegy+'.gz',dir='2004_BP_Vel_Benchmark',
      server='ftp://software.seg.org',top='pub/datasets/2D')
Flow(velsegy,velsegy+'.gz','zcat $SOURCE',stdin=0)

def getvel(vel):
    Flow(vel,velsegy,
         '''
         segyread read=data |
         put label=Velocity unit=m/s
         o1=0 d1=6.25 label1=Depth    unit1=m
         o2=0 d2=12.5 label2=Distance unit2=m
         ''')

shots = [
    'shots0001_0200',
    'shots0201_0400',
    'shots0401_0600',
    'shots0601_0800',
    'shots0801_1000',
    'shots1001_1200',
    'shots1201_1348'
    ]

for shot in shots:
    Fetch(shot+'.segy.gz',dir='2004_BP_Vel_Benchmark',
          server='ftp://software.seg.org',top='pub/datasets/2D')
    Flow(shot+'.segy',shot+'.segy.gz','zcat $SOURCE',stdin=0)
    Flow([shot,'t-'+shot,shot+'.asc',shot+'.bin'],shot+'.segy',
         '''
         segyread tfile=${TARGETS[1]} hfile=${TARGETS[2]} bfile=${TARGETS[3]}
         ''')
Flow('bpshots',shots,'cat axis=2 ${SOURCES[1:%d]}' % len(shots))
Flow('tbpshots',map(lambda x: 't-'+x,shots),
     'cat axis=2 ${SOURCES[1:%d]}' % len(shots))

def getshots(data):
    Flow(data,'bpshots tbpshots',
         '''
         intbin head=${SOURCES[1]} xk=tracf yk=fldr |
         put label3=Shot   o3=0.05 d3=0.05   unit3=km
             label2=Offset o2=-15  d2=0.0125 unit2=km
         ''')

def zodata(zo):
    Flow(zo,'bpshots tbpshots',
         '''
         intbin head=${SOURCES[1]} xk=tracf yk=fldr |
         window n2=1 f2=-1 | cut max1=0.1 | 
         put label2=Distance o2=0.05 d2=0.05   unit2=km
         ''')
