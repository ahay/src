from rsf.proj import *

def slow(name,par):
    Flow('slow_'+name,'vel_'+name,
         'window n2=%(nx)d min2=%(ox)g | math output=1/input ' % par)
    
