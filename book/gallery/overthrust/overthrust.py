from rsf.proj import *
# from rsf.recipes.beg import server

# fetch data and resample from m to km
# Fetch(['overthrust.vites.h','overthrust.vites'],'overthrust',server)

tgz = 'Overthrust_3D_CD1.tar.gz'

Fetch(tgz,'seg_eage_models_cd',
      server='https://s3.amazonaws.com',
      top='open.source.geoscience/open_data')

overthrust = ['./Overthrust_3D_CD1/3-D_Overthrust_Model_Disk1/3D-Velocity-Grid/overthrust.vites.h',
              './Overthrust_3D_CD1/3-D_Overthrust_Model_Disk1/3D-Velocity-Grid/overthrust.vites'] 

zcat = WhereIs('gzcat') or WhereIs('zcat')

Flow(overthrust,tgz,
     zcat + ' $SOURCE | tar -xvf - $TARGETS',stdin=0,stdout=-1,suffix='')

Flow('overthrust',overthrust,
     '''
     (cat ${SOURCES[0]} ; echo in=${SOURCES[1]} data_format=xdr_float) |
     dd form=native |
     scale dscale=0.001 |
     put label1=X label2=Y label3=Z unit1=km unit2=km unit3=km
     label=Velocity unit=km/s d1=0.025 d2=0.025 d3=0.025
     ''',stdin=0)

def getvel(vel):
     Flow(vel,'overthrust','transp plane=23 memsize=1000 | transp plane=12 memsize=1000')

def getvel2D(vel):
    # 2-D slice
    Flow(vel,'overthrust','window n2=1 f2=400 | transp')

