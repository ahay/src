from rsf.proj import *

def getvel(vel):
    Fetch('marmvel.hh','marm')

    Flow(vel,'marmvel.hh',
         '''
         dd form=native | 
         scale rscale=.001 | put d1=0.004 d2=0.004
         label1=Depth label2=Distance unit1=km unit2=km
         label=Velocity unit=km/s
         ''')

def get_zodata(data):
    Fetch('marmexp.hh','marm')

    Flow(data,'marmexp.hh','dd form=native')

def get_shots(data):
    Fetch('marmrefl.hh','marm')

    Flow(data,'marmrefl.hh',
         '''
         dd form=native | 
         put label1=Time unit1=s 
         d2=0.025 o2=-2.575 label2=Offset unit2=km 
         d3=0.025  o3=3     label3=Shot   unit3=km
         ''')

def get_ffd_shots(data):
    Fetch('ffd_data.rsf','marm')

    Flow(data,'ffd_data.rsf','dd form=native | put label2=Receiver label3=Shot')

    
