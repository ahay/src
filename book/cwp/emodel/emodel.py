from rsfproj import *
import fdmod

def execute(vp,par):

    # ------------------------------------------------------------
    Result(vp,fdmod.cgrey('wantaxis=n color=j allpos=y',par))
    Plot  (vp,fdmod.cgrey('wantaxis=n color=j allpos=y',par))
    # ------------------------------------------------------------

    # HWT
    Flow('hwt',vp,
         '''
         smooth rect1=10 rect2=10 repeat=5 |
         hwt2d verb=n xsou=%(xsou)g zsou=%(zsou)g
         nt=%(nt)d ot=%(ot)g dt=%(dt)g
         ng=%(ng)d og=%(og)g dg=%(dg)g
         ''' % par)

    fdmod.rayplot('hwt',10,50,2,50,par)
    Result('hwtvel',[vp,'hwt'],'Overlay')
    
#    Plot('ray','hwt','window j1=40 j2=50 | transp |'+ fdmod.cgraph('wantaxis=n plotcol=0 plotfat=3',par))
#    Plot('wft','hwt','window j1=1  j2=500 |'        + fdmod.cgraph('wantaxis=n plotcol=0 plotfat=3 symbol=.',par))
#    Result('hwt',[vp,'ray','wft'],'Overlay')
    
    # FME
    Flow('fme',vp,'eikonal zshot=0 yshot=%(xsou)g' % par)    
    Plot('fme',fdmod.ccont('plotcol=0 dc=0.25',par))
    Result('fme',[vp,'fme'],'Overlay')

    # ------------------------------------------------------------
    # salt mask
    Flow(  'smask',vp,'mask min=4.2 | dd type=float')
    Result('smask',fdmod.cgrey('allpos=y',par))

    # water mask
    Flow(  'wmask',vp,'mask max=1.5 | dd type=float')
    Result('wmask',fdmod.cgrey('allpos=y',par))
    
    Flow(  'amask',vp,'mask min=2.0 max=2.2 | dd type=float')
    Result('amask',fdmod.cgrey('allpos=y',par))
    
    Flow(  'bmask',vp,'mask min=3.0 max=3.5 | dd type=float')
    Result('bmask',fdmod.cgrey('allpos=y',par))    

    # ------------------------------------------------------------
    # vp/vs ratio
    Flow('vpvs_',None,
         '''
         math n1=%(nz)d o1=%(oz)g d1=%(dz)g output="3-0.15*x1" |
         spray axis=2 n=%(nx)d o=%(ox)g d=%(dx)g
         ''' % par)
    Flow('vpvs','vpvs_ smask wmask',
         '''
         math output="input*(1-s)*(1-w)+1.73*s+w"
         s=${SOURCES[1]}
         w=${SOURCES[2]}
         ''' % par)
    Result('vpvs',fdmod.cgrey('allpos=y',par))  

    # ------------------------------------------------------------
    # vs
    Flow('vs',[vp,'vpvs','wmask'],
         '''
         math output="input/r*(1-w)"
         r=${SOURCES[1]}
         w=${SOURCES[2]}
         ''')
    Result('vs',fdmod.cgrey('color=j allpos=y',par))

    # ------------------------------------------------------------
    # ro
    Flow('ro',[vp,'smask','wmask'],
         '''
         math output="(2200000+(input-1.5)/3.6*100000)*(1-w)*(1-s)+s*1730000+w*1000000"
         s=${SOURCES[1]}
         w=${SOURCES[2]}
         ''')
    Result('ro',fdmod.cgrey('color=j allpos=y',par))

    # ------------------------------------------------------------
    # epsilon
    Flow('epsilon',['amask','bmask'],
         '''
         math output="input*0.1+b*0.25"
         b=${SOURCES[1]}
         ''')
    Result('epsilon',fdmod.cgrey('color=j allpos=y',par))

    # ------------------------------------------------------------
    # delta
    Flow('delta',['amask','bmask'],
         '''
         math output="input*0.05+b*0.1"
         b=${SOURCES[1]}
         ''')
    Result('delta',fdmod.cgrey('color=j allpos=y',par))

    # ------------------------------------------------------------
    # stiffness tensor
    fdmod.isotropic(  'ci',vp,'vs','ro',par)
    fdmod.anisotropic('ca',vp,'vs','ro','epsilon','delta',par)

    # ------------------------------------------------------------
    fdmod.wavelet('wav_',par['frq'],par)
    # ------------------------------------------------------------
    # acoustic source
    Flow(  'wava','wav_','transp')
    Result('wava','transp | window n1=1000 |' + fdmod.waveplot('',par))
    # ------------------------------------------------------------
    Flow('ver','wav_','math output=input*1')
    Flow('hor','wav_','math output=input*1')
    Flow('wave',['ver','hor'],
         '''
         cat axis=2 space=n ${SOURCES[1:2]} |
         transp plane=12 |
         transp plane=23 |
         transp plane=12
         ''')
    Plot('ver','wave','window n2=1 f2=0 n3=1000 |' + fdmod.waveplot('',par))
    Plot('hor','wave','window n2=1 f2=1 n3=1000 |' + fdmod.waveplot('',par))
    Result('wave','ver hor','Movie')
    
    # ------------------------------------------------------------
    # source/receiver coordinates
    fdmod.point('ss',par['xsou'],par['zsou'],1,par)
    fdmod.horizontal('rr',par['oz'],par)

    Plot('rr','window n1=2 | dd type=complex | window j2=10 | '
    + fdmod.cgraph('symbol=o plotcol=1',par))
    Plot('ss','window n1=2 | dd type=complex | window | '
    + fdmod.cgraph('symbol=x plotcol=2',par))

    # ------------------------------------------------------------
    # acoustic modeling
    fdmod.awefd('da','wa','wava',vp,'ro','ss','rr','',par)
    Result('wa',fdmod.wgrey('wantaxis=n pclip=99',par))
    Result('da','window j2=8 | transp |'
           + fdmod.dgrey('wantaxis=n pclip=97',par))

    fdmod.wom('wam','wa','vp',1.5,par)
    Result(   'wam',fdmod.wgrey('wantaxis=n pclip=99',par))

    # ------------------------------------------------------------
    # elastic modeling
    fdmod.ewefd('di','wi','wave','ci','ro','ss','rr','ssou=y opot=n',par)
    fdmod.ewefd('de','we','wave','ca','ro','ss','rr','ssou=y opot=n',par)

    for k in ('i'):
        
        Flow('w1'+k,'w'+k,
             'window n3=1 f3=0 | window min1=%g max1=%g min2=%g max2=%g'
             % (par['zmin'],par['zmax'],par['xmin'],par['xmax']))
        Flow('w2'+k,'w'+k,
             'window n3=1 f3=1 | window min1=%g max1=%g min2=%g max2=%g'
             % (par['zmin'],par['zmax'],par['xmin'],par['xmax']))
        
        fdmod.wom('w1'+k+'m','w1'+k,'vp',1.0,par)
        fdmod.wom('w2'+k+'m','w2'+k,'vp',1.0,par)
        
        Result('w1'+k+'m',fdmod.wgrey('wantaxis=n pclip=99',par))
        Result('w2'+k+'m',fdmod.wgrey('wantaxis=n pclip=99',par))
        
        Result('d1'+k,'d'+k,'window n2=1 f2=0 j3=8 | transp | bandpass flo=1 |'
               + fdmod.dgrey('wantaxis=n pclip=98',par))
        Result('d2'+k,'d'+k,'window n2=1 f2=1 j3=8 | transp | bandpass flo=1 |'
               + fdmod.dgrey('wantaxis=n pclip=98',par))
