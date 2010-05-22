from rsf.proj import *

def Bay(name,f2=500,n2=1600,f1=400,n1=700):
    '''Extracts the Bay Area dataset (digital topography) and windows it'''
    Fetch('bay.h','bay')

    Flow(name,'bay.h',
         '''
         dd form=native |
         window f2=%d n2=%d f1=%d n1=%d |
         reverse which=1 |
         costaper nw1=50 nw2=50
         ''' % (f2,n2,f1,n1))
