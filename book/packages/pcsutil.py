from rsfproj import *
from math import *

# make a model with a dipping linear interface
# defined by coordinates [start](s1,s2) and [end](e1,e2)
def dipline(mod,s1,s2,e1,e2,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    layers = (
        ((s2,s1),(e2,e1)),
        ((min2,max1),(max2,max1))
        )

    ra = (e1-s1)/(e2-s2)
    vels = "%s,%s,%s" %(vi,vt,vt)
    drvs = "%s,%s" %(tan(ra),tan(ra))

    dim1 = 'd1=%g o1=%g n1=%d' % (d2,o2,n2)
    dim2 = 'd2=%g o2=%g n2=%d' % (d1,o1,n1)

    for i in range(len(layers)):
        inp = mod+'inp%d.rsf' % (i+1)
        Flow('./'+inp,None,
             'echo %s in=$TARGET data_format=ascii_float n1=2 n2=%d' % \
             (string.join(map(lambda x: string.join(map(str,x)),layers[i]),
                          ' '),
              len(layers[i])))
        
    Flow(mod+'lay1',mod+'inp1','dd form=native | spline %s fp=%s' % (dim1,drvs))
    Flow(mod+'lay2',mod+'inp2','dd form=native | spline %s fp=%s' % (dim1,drvs))

    Flow(    mod+'layers',[mod+'lay1',mod+'lay2'],'cat axis=2 ${SOURCES[1:2]}')
    Flow(mod,mod+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,n1,d1,o1) )

def polar(cart,polar,par):
    
    Flow(cart+'rad',None,
         '''
         math
         o1=-90 d1=0.6 n1=301
         o2=-90 d2=0.6 n2=301
         output="sqrt(x1*x1+x2*x2)"
         ''')
    Flow(cart+'ang',None,
         '''
         math
         o1=-90 d1=0.6 n1=301
         o2=-90 d2=0.6 n2=301
         output="x2&x1"
         ''')
    Flow(cart+'coord',[cart+'rad',cart+'ang'],
         'cat ${SOURCES[1]} | transp plane=13')

    Flow(cart,[polar,cart+'coord'],
         '''
         put d2=.0174444 o2=-3.14 |
         inttest2 coord=${SOURCES[1]} interp=cub nw=4 |
         put o1=-1 d1=.006666 o2=-1 d2=.006666
         ''')


