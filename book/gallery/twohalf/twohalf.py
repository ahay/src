from rsf.proj import *

Fetch('model.segy',dir='1997_2.5D',server='http://software.seg.org',top='datasets/2D')

def getvel(model):
    Flow(model,'model.segy',
         '''
         segyread read=data |
         put
         d1=12.5 o1=0 label1=Depth unit1=m
         d2=12.5 o2=0 label2=Distance unit2=m
         label=Velocity unit=m/s
         ''')

shots = '1997_2.5D_shots.segy'
shotsgz = shots+'.gz'

Fetch(shotsgz,dir='1997_2.5D',server='ftp://software.seg.org',top='/pub/datasets/2D')
Flow(shots,shotsgz,'zcat $SOURCE',stdin=0)

def getshots(data):
    Flow(data,shots,
         '''
         segyread read=data |
         put
         n2=256 d2=25 o2=0 label2=Offset unit2=m
         n3=385 d3=50 o3=0 label3=Shot   unit3=m
         ''')


