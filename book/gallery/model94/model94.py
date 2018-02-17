from rsf.proj import *

Fetch('velocity.segy.gz',dir='bpmodel94',
      server='https://s3.amazonaws.com',top='open.source.geoscience/open_data')
Flow('velocity.segy','velocity.segy.gz','gunzip -c $SOURCE',stdin=0)

def get_vel(velocity):
    Flow(velocity,'velocity.segy',
         '''
         segyread read=data | scale dscale=0.001 | 
         put o1=-2 unit1=km label1=Depth 
         d2=0.015 unit2=km label2=Distance label=Velocity unit=km/s
         ''')

shots = 'Model94_shots.segy'

Fetch(shots+'.gz',dir='bpmodel94',
      server='https://s3.amazonaws.com',top='open.source.geoscience/open_data')

Flow(shots,shots+'.gz','gunzip -c $SOURCE',stdin=0)
Flow('shots94 tshots94',shots,
     'segyread tfile=${TARGETS[1]}')

def get_shots(shots):
    Flow(shots,'shots94',
         '''
         intbin xk=tracf yk=fldr | 
         put label2=Offset unit2=km d2=0.015 o2=-3.6 
         label3=Shot unit3=km d3=0.09 o3=0
         ''')

def get_shot_headers(heads):
    Flow(heads,'tshots94',
         '''
         intbin xk=tracf yk=fldr head=$SOURCE | 
         put label2=Offset unit2=km d2=0.015 o2=-3.6 
         label3=Shot unit3=km d3=0.09 o3=0
         ''')


