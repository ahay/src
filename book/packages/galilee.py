from rsfproj import *

def Galilee(name,nx=280,ny=440,interp=1):
    '''Extracts the Sea of Galiee dataset and bins it'''

    Fetch('galilee.h','galilee')

    base = -212

    Flow('data','galilee.h','dd form=native')
    Flow('mask','data','window n1=1 f1=2 | mask max=%g' % base)
    Flow('triplets','data mask','headerwindow mask=${SOURCES[1]}')

    Flow(name,'triplets',
         '''
         window n1=1 f1=2 |
         math output=%g-input |
         bin interp=%d xkey=0 ykey=1 head=$SOURCE
         xmin=198 xmax=212 ymin=233 ymax=257
         nx=%d ny=%d
         ''' % (base,interp,nx,ny))
 
    
