try:    from rsf.cluster import *
except: from rsf.proj    import *
import math
## 
 # mapping from polar to Cartesian coordinates
 # assumes angles measured in degrees!
 # uses the trigonometric convention (counter-clockwise angle variation)
 ##

# ------------------------------------------------------------
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

# ------------------------------------------------------------
def ovl(ovl,jc,jr,custom,cco):

    min=cco['o']
    max=cco['o']+(cco['n']-1)*cco['d']

    # circles
    for ic in range(jc,int((cco['n']-1)/2*cco['d']+jc),jc):
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
             wantaxis=n yreverse=y screenratio=1 plotcol=3
             min1=%d max1=%d min2=%d max2=%d
             %s
             ''' %(min,max,min,max,
                   custom))

        Plot(ovl+ctag+'l',ovl+ctag,
             'box x0=%g y0=%g label="%s\^o\_" xt=%g yt=%g lab_fat=1 lab_color=3 boxit=n'%
             ((5.6+2.75*ic/90.),(5.1-2.75*ic/90.),"%s"%ic,0,0),stdin=0)

    # radii
    for ir in range(0,360,jr):
        rtag='%03d'%ir
        
        Flow(ovl+'_rx'+rtag,None,
             'math n1=%d d1=%g o1=%g output="x1*cos(3.1415*%g/180)"'
             % ( (cco['n']-1)/2+1-jc/cco['d'],cco['d'],jc,ir) )
        Flow(ovl+'_rz'+rtag,None,
             'math n1=%d d1=%g o1=%g output="x1*sin(3.1415*%g/180)"'
             % ( (cco['n']-1)/2+1-jc/cco['d'],cco['d'],jc,ir) )
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
             wantaxis=n yreverse=y screenratio=1 plotcol=3
             min1=%d max1=%d min2=%d max2=%d
             %s
             ''' %(min,max,min,max,
                   custom))

    Plot(ovl+'-ann000',ovl+'000',
         'box x0=%g y0=%g label="%s" xt=%g yt=%g lab_fat=1 boxit=n'%
         (9.65,5,"0\^o\_",0,0),stdin=0)
    Plot(ovl+'-ann090',ovl+'000',
         'box x0=%g y0=%g label="%s" xt=%g yt=%g lab_fat=1 boxit=n'%
         (5.55,9.25,"90\^o\_",0,0),stdin=0)
    Plot(ovl+'-ann180',ovl+'000',
         'box x0=%g y0=%g label="%s" xt=%g yt=%g lab_fat=1 boxit=n'%
         (1.2,5,"180\^o\_",0,0),stdin=0)
    Plot(ovl+'-ann270',ovl+'000',
         'box x0=%g y0=%g label="%s" xt=%g yt=%g lab_fat=1 boxit=n'%
         (5.5,1,"270\^o\_",0,0),stdin=0)

    Plot(ovl+'-center',ovl+'000',
         'box x0=%g y0=%g label="%s" xt=%g yt=%g lab_fat=1 boxit=n'%
         (5.6,5.10,"0",0,0),stdin=0)
    
    Plot(ovl,
         map(lambda x: ovl+'%02d' % x,range(jc,int((cco['n']-1)/2*cco['d']+jc),jc))+
         map(lambda x: ovl+'%03d' % x,range(0,360,jr))+
         map(lambda x: ovl+'%02dl'% x,range(jc,int((cco['n']-1)/2*cco['d']+jc),jc))+
         [ovl+'-ann000',ovl+'-ann090',ovl+'-ann180',ovl+'-ann270'],
         'Overlay')

#    Plot(ovl,
#         map(lambda x: ovl+'%02d' % x,range(jc,int((cco['n']-1)/2*cco['d']+jc),jc))+
#         map(lambda x: ovl+'%03d' % x,range(0,360,jr)),
#         'Overlay')

# ------------------------------------------------------------
def polplot(custom,par):
    return '''
    grey title="" pclip=99.9 allpos=n
    screenratio=0.25 screenht=3.5
    label1="\F10 q\F3" unit1="\^o\_"
    label2="\F10 f\F3" unit2="\^o\_"
    %s
    ''' % (par['labelrot']+custom)

def carplot(custom,par):
    return '''
    grey title="" pclip=99.9 allpos=n
    screenratio=1 wantaxis=n
    %s
    ''' %(par['labelrot']+custom)

# ------------------------------------------------------------
def ann(ann,tht,phi,custom,cco):

    min=cco['o']
    max=cco['o']+(cco['n']-1)*cco['d']

    Flow(ann,None,
         'spike nsp=2 n1=2 k1=1,2 mag=%g,%g'
         %(+tht*math.cos(3.1415*phi/180),
           -tht*math.sin(3.1415*phi/180)))   
    Plot(ann+'bg',ann,
         '''
         dd type=complex |
         graph title=""
         wantaxis=n yreverse=y screenratio=1 plotcol=0
         symbol=. symbolsz=25 plotfat=20
         min1=%d max1=%d min2=%d max2=%d
         %s plotcol=7
         ''' %(min,max,min,max,
               custom))
    Plot(ann+'fg',ann,
         '''
         dd type=complex |
         graph title=""
         wantaxis=n yreverse=y screenratio=1 plotcol=0
         symbol=. symbolsz=20 plotfat=10
         min1=%d max1=%d min2=%d max2=%d
         %s
         ''' %(min,max,min,max,
               custom))
    Plot(ann,[ann+'bg',ann+'fg'],'Overlay')
    
