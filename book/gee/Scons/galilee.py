def Galilee(name,grad=1):
    '''Extracts the Sea of Galiee dataset and bins it'''
    
    datdir  = 'ftp://begpc132.beg.utexas.edu/data/galilee'
    Fetch('galilee.H',datdir)

    base = -212

    Flow('data','galilee.H','dd data_format=native_float')
    Flow('mask','data','window n1=1 f1=2 | mask max=%g' % base)
    Flow('triplets',['data','mask'],'headerwindow mask=${SOURCES[1]}')

    Flow(name,'triplets',
         '''
         window n1=1 f1=2 |
         math output=%g-input |
         bin interp=1 xkey=0 ykey=1 head=$SOURCE
         xmin=198 xmax=212 ymin=233 ymax=257
         nx=280 ny=440
         ''' % base)
 
    
