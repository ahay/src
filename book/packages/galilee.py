from rsf.proj import *

base = -212

Fetch('galilee.h','galilee')
Flow('data','galilee.h','dd form=native')
Flow('mask','data','window n1=1 f1=2 | mask max=%g' % base)
Flow('triplets','data mask','headerwindow mask=${SOURCES[1]}')

def Galilee(name,nx=280,ny=440,interp=1,
            xmin=198,xmax=212,ymin=233,ymax=257,shape='',extra=''):
    '''Extracts the Sea of Galiee dataset and bins it'''
    global base

    Flow(name,'triplets',
         '''
         window n1=1 f1=2 |
         math output=%g-input |
         %sbin interp=%d xkey=0 ykey=1 head=$SOURCE
         xmin=%d xmax=%d ymin=%d ymax=%d nx=%d ny=%d %s
         ''' % (base,shape,interp,xmin,xmax,ymin,ymax,nx,ny,extra))
 
