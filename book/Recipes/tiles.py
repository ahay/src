from rsf.proj import *

def plot(plot,data,n1,n2,f1,f2,j1,j2,w1,w2,custom,par):

    Flow(plot+'_wind',
         data,
         '''
         window n1=%d j1=%d f1=%d n2=%d j2=%d f2=%d |
         transp plane=23 |
         transp plane=12 |
         transp plane=34
         ''' % (n1/j1,j1,f1,
                n2/j2,j2,f2) )

    Flow(plot+'_tile',
         plot+'_wind',
         '''
         put
         n1=%d n2=%d n3=1 n4=1
         o1=%d o2=%d o3=0 o4=0
         d1=%g d2=%g d3=1 d4=1
         ''' % (n1/j1*w1,
                n2/j2*w2,
                f1/j1,
                f2/j2,
                1./w1,
                1./w2))

    ratio = 1.0*(n1/j1*w1)/(n2/j2*w2);
    height= ratio*14
    
    Result(plot,plot+'_tile',
           '''
           grey title="" pclip=100 labelsz=4
           label1="" unit1=" " label2="" unit2=" "
           grid=y g1num=1 g2num=1 gridcol=0
           screenratio=%g screenht=%g
           %s
           ''' % (ratio,height,custom))
