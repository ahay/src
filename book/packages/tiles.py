from rsfproj import *

def plot(data,n1,n2,j1,j2,nh1,nh2,custom):

    Flow(data+'_wind',
         data,
         '''
         window j1=%d f1=%d j2=%d f2=%d |
         transp plane=23 |
         transp plane=12 |
         transp plane=34
         ''' % (j1,0,j2,0) )

    Flow(data+'_tile',
         data+'_wind',
         '''
         put n1=%d n2=%d n3=1 n4=1
         o1=0  o2=0  o3=0 o4=0
         d1=%g d2=%g d3=1 d4=1
         ''' % (n1/j1*nh1,
                n2/j2*nh2,
                1./nh1,
                1./nh2))

    zmin=0
    zmax=zmin+n1/j1
    xmin=0
    xmax=xmin+n2/j2
    
    ratio = 1.0*(zmax-zmin)/(xmax-xmin);
    height= ratio*14

    Result(data,data+'_tile',
           '''
           grey title=""
           label1="" unit1="" label2="" unit2=""
           pclip=99.9
           grid=y g1num=1 g2num=1 gridcol=0
           min1=%g max1=%g
           min2=%g max2=%g
           screenratio=%g screenht=%g
           %s
           ''' % (zmin,zmax,xmin,xmax,ratio,height,custom))
