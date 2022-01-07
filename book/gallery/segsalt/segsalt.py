from rsf.proj import *
from rsf.recipes.beg import server

tgz = 'Salt_Model_3D.tar.gz'

Fetch(tgz,'seg_eage_models_cd',
      server='https://s3.amazonaws.com',
      top='open.source.geoscience/open_data')

zcat = WhereIs('gzcat') or WhereIs('zcat')

saltaa = ['./Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/saltaa.h',
          './Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/saltaa@']

Flow(saltaa,tgz,
     zcat + ' $SOURCE | tar -xvf - $TARGETS',stdin=0,stdout=-1,suffix='')

# Get 2-D velocity model
# Fetch(['saltaa.h','saltaa@'],'segsalt',server)
Flow('saltaa',saltaa,
     '''
     (cat ${SOURCES[0]} ; echo in=${SOURCES[1]} data_format=xdr_float) |
     put n1=780 |
     dd form=native | transp | 
     put label1=Z unit1=m label2=X unit2=m label=Velocity unit=m/s
     ''',stdin=0)

def getvel2D(vel):
    Flow(vel,'saltaa',
         '''
         math output=1500 |
         cat axis=1 $SOURCE |
         put o1=-4014.281 |
         remap1 o1=0 d1=20 n1=201 
         ''')
    
# Get 3-D velocity model
salt = ['./Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/salt.h',
        './Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/SALTF.ZIP',
        './Saltf@@']
Flow(salt[:2],tgz,
     zcat + ' $SOURCE | tar -xvf - $TARGETS',stdin=0,stdout=-1,suffix='')
Flow(salt[2],salt[1],'unzip $SOURCE',stdin=0,stdout=-1,suffix='')

#Fetch(['salt.h','Saltf@@'],'segsalt',server)  
Flow('salt',[salt[0],salt[2]],
     '''
     (cat ${SOURCES[0]} ; echo in=${SOURCES[1]} data_format=xdr_float) | 
     dd form=native | 
     put label1=X unit1=m label2=Y unit2=m label3=Z unit3=m
     label=Velocity unit=m/s
     ''',stdin=0)

def getvel3D(vel):
    Flow(vel,'salt',
         '''
         transp plane=23 memsize=1000 | 
         transp plane=12 memsize=1000 |
         put d1=0.02 unit1=km d2=0.02 unit2=km d3=0.02 unit3=km |
         scale dscale=0.001 
         ''')

