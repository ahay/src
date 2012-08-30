try:    from rsf.cluster import *
except: from rsf.proj    import *
from math import *

# ------------------------------------------------------------
def retiredplane(mod,s1,s2,aa,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    e2 = max2
    e1 = s1 + aa * (e2-s2)

    layers = (
        ((s2,s1),(e2,e1)),
        ((min2,min1),(max2,max1))
        )

    vels = "%s,%s,%s" %(vi,vt,vt)
    drvs = "%s,%s" %(aa,aa)

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

    Flow(mod+'1',mod+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,n1,d1,o1) )
    Flow(mod+'2',mod+'1',
         '''
         pad beg1=1 | window n1=%d
         ''' % (n1) )
    Flow(mod,[mod+'1',mod+'2'],'add ${SOURCES[1]} scale=1,-1 | scale axis=123')

# ------------------------------------------------------------
def plane(mod,s1,s2,aa,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    e2 = max2
    e1 = s1 + aa * (e2-s2)

#    print min1,min2,max1,max2
#    print s1,s2,e1,e2,aa,aa*(e2-s2)

    vels = "%s,%s,%s" %(vi,vt,vt)
    drvs = "%s,%s" %(aa,aa)

    dim1 = 'd1=%g o1=%g n1=%d' % (d2,o2,n2)
    dim2 = 'd2=%g o2=%g n2=%d' % (d1,o1,n1)

    Flow(mod+'lay2',None,
         '''
         spike nsp=4 mag=%g,%g,%g,%g
               n1=4 n2=1 k1=1,2,3,4 |
         put n1=2 n2=2 |
         dd form=native |
         spline %s fp=%s
         '''%(min2,min1,max2,max1,dim1,drvs))

    Flow(mod+'lay1',None,
         '''
         spike nsp=4 mag=%g,%g,%g,%g
               n1=4 n2=1 k1=1,2,3,4 |
         put n1=2 n2=2 |
         dd form=native |
         spline %s fp=%s
         '''%(s2,s1,e2,e1,dim1,drvs))
     
    Flow(    mod+'layers',[mod+'lay1',mod+'lay2'],'cat axis=2 ${SOURCES[1:2]}')

    Flow(mod+'1',mod+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,n1,d1,o1) )
    Flow(mod+'2',mod+'1',
         '''
         pad beg1=1 | window n1=%d
         ''' % (n1) )
    Flow(mod,[mod+'1',mod+'2'],'add ${SOURCES[1]} scale=1,-1 | scale axis=123')

# ------------------------------------------------------------
# make a model with a dipping linear interface
# defined by coordinates [start](s1,s2) and [end](e1,e2)
def dipline(mod,s1,s2,e1,e2,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    layers = (
        ((s2,s1),(e2,e1)),
        ((min2,min1),(max2,max1))
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

# ------------------------------------------------------------
def angline(mod,s1,s2,aa,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    ll = 1000*sqrt(d1*d1 + d2*d2)

    ra = aa/180.*pi
    e1 = s1 + ll*sin(ra)
    e2 = s2 + ll*cos(ra)

    layers = (
        ((s2,s1),(e2,e1)),
        ((min2,min1),(max2,max1))
        )

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

# ------------------------------------------------------------
# make a polar plot of a 2D Cartesian map
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
         o1=+90 d1=-0.6 n1=301
         o2=+90 d2=-0.6 n2=301
         output="-x2&x1"
         ''')
    Flow(cart+'coord',[cart+'rad',cart+'ang'],
         'cat ${SOURCES[1]} | transp plane=13')

    Flow(cart,[polar,cart+'coord'],
         '''
         put d2=.0174444 o2=-3.14 |
         inttest2 coord=${SOURCES[1]} interp=cub nw=4 |
         put o1=-1 d1=.006666 o2=-1 d2=.006666
         ''')

# ------------------------------------------------------------
def plotpolar(cc,cz,cx,custom):

    Flow(cc,[cz,cx],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} |
         transp |
         dd type=complex
         ''', stdin=0)
    Plot(cc,
         '''
         graph labelrot=n wantaxis=n title="" yreverse=y
         min1=-90 max1=+90
         min2=-90 max2=+90
         plotcol=4 symbol="." symbolsz=3
         screenratio=1
         %s
         ''' % custom)

# ------------------------------------------------------------
def polargrid(cc,custom,par):
    
    Flow(cc+'t_',None,'math n1=60  d1=1 o1=30 output="x1"')
    for k in range(0,360,30):
        Flow(cc+'tz'+str(k),cc+'t_','math output="input*sin(%g)" ' % (k*pi/180.))
        Flow(cc+'tx'+str(k),cc+'t_','math output="input*cos(%g)" ' % (k*pi/180.))
        plotpolar(cc+'t'+str(k),cc+'tz'+str(k),cc+'tx'+str(k),custom+'')
        
    allt =  map(lambda x: cc+'t%d.vpl' % x,range(0,360,30))
    Plot(cc+'t',allt,'Overlay')

    Flow(cc+'p_',None,'math n1=360 d1=1 o1=0  output="x1*%g"' % (pi/180.))
    for l in range(30,91,30):
        Flow(cc+'pz'+str(l),cc+'p_','math output="%d*sin(input)" ' % l)
        Flow(cc+'px'+str(l),cc+'p_','math output="%d*cos(input)" ' % l)
        plotpolar(cc+'p'+str(l),cc+'pz'+str(l),cc+'px'+str(l),custom+'')

    allp =  map(lambda x: cc+'p%d.vpl' % x,range(30,91,30))
    Plot(cc+'p',allp,'Overlay')
    
    Plot(cc,[cc+'t',cc+'p'],'Overlay')
