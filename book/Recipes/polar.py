from rsf.proj import *

## 
 # mapping from polar to Cartesian coordinates
 # assumes angles measured in degrees!
 # uses the trigonometric convention (counter-clockwise angle variation)
 ##

def p2c(pol,car,cco):

    # radius
    Flow(car+'_rad',None,
         '''
         math
         n1=%(n)d o1=%(o)g d1=%(d)g
         n2=%(n)d o2=%(o)g d2=%(d)g
         output="sqrt(x1*x1+x2*x2)"
         '''%cco)

    # angle (in degrees)
    Flow(car+'_ang',None,
         '''
         math
         n1=%(n)d o1=%(o)g d1=%(d)g
         n2=%(n)d o2=%(o)g d2=%(d)g
         output="-(x2&x1)/3.1415*180"
         '''%cco)
    
    Flow(car+'_coo',[car+'_rad',car+'_ang'],
         'cat ${SOURCES[1]} | transp plane=13 | put o1=0 d1=1')
    
    Flow(car,[pol,car+'_coo'],
         'inttest2 coord=${SOURCES[1]} interp=cub nw=4')

## 
 # TODO: write c2p(cart,polar,cco)
 ##

def ovl(ovl,jc,jr,custom,cco):

    min=cco['o']
    max=cco['o']+(cco['n']-1)*cco['d']

    for ic in range(jc,(cco['n']-1)/2+jc,jc):
        ctag='%02d'%ic

        Flow(ovl+'_cx'+ctag,None,
             'math n1=%d d1=%g o1=%g output="%g*cos(3.1415*x1/180)"'
             % (361,1,0,ic) )
        Flow(ovl+'_cz'+ctag,None,
             'math n1=%d d1=%g o1=%g output="%g*sin(3.1415*x1/180)"'
             % (361,1,0,ic) )
        Flow(ovl+ctag,[ovl+'_cx'+ctag,ovl+'_cz'+ctag],
             '''
             cat axis=2 space=n
             ${SOURCES[0]} ${SOURCES[1]} | transp
             ''', stdin=0)
        Plot(ovl+ctag,
             '''
             dd type=complex |
             window |
             graph title=""
             wantaxis=n yreverse=y screenratio=1 plotcol=2
             min1=%d max1=%d min2=%d max2=%d
             %s
             ''' %(min,max,min,max,
                   custom))
    
    for ir in range(0,360,jr):
        rtag='%03d'%ir

        Flow(ovl+'_rx'+rtag,None,
             'math n1=%d d1=%g o1=%g output="x1*cos(3.1415*%g/180)"'
             % (cco['n'],cco['d'],cco['o'],ir) )
        Flow(ovl+'_rz'+rtag,None,
             'math n1=%d d1=%g o1=%g output="x1*sin(3.1415*%g/180)"'
             % (cco['n'],cco['d'],cco['o'],ir) )
        Flow(ovl+rtag,[ovl+'_rx'+rtag,ovl+'_rz'+rtag],
             '''
             cat axis=2 space=n
             ${SOURCES[0]} ${SOURCES[1]} | transp
             ''', stdin=0)
        
        Flow(ovl+rtag,[ovl+'_rx'+rtag,ovl+'_rz'+rtag],
             '''
             cat axis=2 space=n
             ${SOURCES[0]} ${SOURCES[1]} | transp
             ''', stdin=0)
        Plot(ovl+rtag,
             '''
             dd type=complex |
             window |
             graph title=""
             wantaxis=n yreverse=y screenratio=1 plotcol=2
             min1=%d max1=%d min2=%d max2=%d
             %s
             ''' %(min,max,min,max,
                   custom))
        
    Plot(ovl,
         map(lambda x: ovl+'%02d' % x,range(jc,(cco['n']-1)/2+jc,jc))+
         map(lambda x: ovl+'%03d' % x,range(0,360,jr)),
         'Overlay')
