from rsf.proj import *
import fdmod

def param():
    par = {
        'nx':501,  'ox':0, 'dx':0.002,  'lx':'x', 'ux':'km',
        'nz':501,   'oz':0, 'dz':0.002,  'lz':'z', 'uz':'km',
        'nt':5000,  'ot':0, 'dt':0.0002, 'lt':'t', 'ut':'s',
        'kt':150,
        'jsnap':200,
        'nb':100,
        'frq':45,
        'ssou':'y',
        'vscale':1,
        'hscale':0
        }
    fdmod.param(par)

    return par

# ------------------------------------------------------------
# acoustic modeling
def amodel(dat,wfl,  wav,vel,den,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [dat,wfl],[wav,vel,den,sou,rec],
         '''
         awefd2d
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# elastic modeling
def emodel(dat,wfl,  wav,ccc,den,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [dat,wfl],[wav,ccc,den,sou,rec],
         '''
         ewefd2d
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# elastic modeling with free surface
def emodelfs(dat,wfl,wav,ccc,den,sou,rec,custom,par):
    par['fdcustom'] = custom

    Flow( [dat,wfl],[wav,ccc,den,sou,rec],
         '''
         ewefd2d
         verb=y free=y snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def test(vp,vs,ro,epsilon,delta,ss,rr,par):
    # ------------------------------------------------------------
    # source/receiver coordinates
    Plot(rr,'window n1=2 | dd type=complex | window j2=10 | '
         + fdmod.cgraph('wantscalebar=y symbol=o plotcol=1',par))
    Plot(ss,'window n1=2 | dd type=complex | window | '
         + fdmod.cgraph('wantscalebar=y symbol=x plotcol=2',par))

    # ------------------------------------------------------------
    # acoustic source
    fdmod.wavelet('wava0',par['frq'],par)
    Flow(  'wava','wava0','transp')
    Result('wava','transp | window n1=500 |' + fdmod.waveplot('title="Acoustic source"',par))
    
    # ------------------------------------------------------------
    # elastic source
    fdmod.wavelet('hor0',par['frq'],par)
    fdmod.wavelet('ver0',par['frq'],par)
    Flow('hor','hor0','math output=input*%(hscale)g' % par)
    Flow('ver','ver0','math output=input*%(vscale)g' % par)

    Flow('wave0','ver hor','cat axis=2 space=n ${SOURCES[1:2]}')
    Flow('wave','wave0',
         '''
         transp plane=12 |
         transp plane=23 |
         transp plane=12
         ''')
    
    Plot('ver','wave','window n2=1 f2=0 | window n1=500 |' + fdmod.waveplot('title="Elastic vertical source"',par))
    Plot('hor','wave','window n2=1 f2=1 | window n1=500 |' + fdmod.waveplot('title="Elastic horizontal source"',par))
    Result('wave','hor ver','Movie')

    # ------------------------------------------------------------
    Plot(vp,     fdmod.cgrey('wantscalebar=y allpos=y bias=1.0    pclip=100',par))
    Plot(vs,     fdmod.cgrey('wantscalebar=y allpos=y bias=1.0    pclip=100',par))
    Plot(ro,     fdmod.cgrey('wantscalebar=y allpos=y bias=100000 pclip=100',par))
    Plot(epsilon,fdmod.cgrey('wantscalebar=y allpos=y             pclip=100',par))
    Plot(delta,  fdmod.cgrey('wantscalebar=y allpos=y             pclip=100',par))

    Result(vp,     [vp,     ss,rr],'Overlay')
    Result(vs,     [vs,     ss,rr],'Overlay')
    Result(ro,     [ro,     ss,rr],'Overlay')
    Result(epsilon,[epsilon,ss,rr],'Overlay')
    Result(delta,  [delta,  ss,rr],'Overlay')

    fdmod.anisotropic('cc','vp','vs','ro','epsilon','delta',par)

    # ------------------------------------------------------------
    # acoustic modeling
    amodel('da','wa','wava',vp,ro,ss,rr,'',par)
    
    Flow('waw','wa',
         '''
         window min1=%g max1=%g min2=%g max2=%g |
         scale axis=123
         ''' % (par['zmin'],par['zmax'],par['xmin'],par['xmax']))
    
    Result('wa',          fdmod.wgrey('pclip=99 title="Acoustic wavefield"',par))
    Result('da','transp | window f1=%(kt)d | put o1=%(ot)g | pad end1=%(kt)d |' % par
           + fdmod.dgrey('pclip=99 title="Acoustic data" grid=y',par))

    # elastic modeling
    emodel('de','we','wave','cc',ro,ss,rr,'ssou=%(ssou)s opot=n' % par,par)
    
    for i in range(2):
        Flow('we'+str(i+1),'we',
             '''
             window n3=1 f3=%d |
             window min1=%g max1=%g min2=%g max2=%g |
             scale axis=123
             ''' % (i,par['zmin'],par['zmax'],par['xmin'],par['xmax']))
        
        Result('we'+str(i+1),
               fdmod.wgrey('title=u%s pclip=99' % str(i+1),par))
        Result('de'+str(i+1),'de',
               '''
               window n2=1 f2=%d |
               transp |
               window f1=%d | put o1=%g | pad end1=%d |
               ''' % (i,par['kt'],par['ot'],par['kt'])
               + fdmod.dgrey('title=u%s pclip=99 grid=y' %str(i+1),par))
    
    Flow(  'weall','we1 we2','cat axis=1 space=n ${SOURCES[1]}')
    Result('weall',
           '''
           grey title="Elastic Wavefields" wantaxis=y screenratio=%f screenht=8
           gainpanel=a pclip=99
           grid1=y grid2=y g1num=0.25 g2num=0.25
           ''' % (2*par['ratio']) )

    Flow('wall','waw we1 we2','cat axis=1 space=n ${SOURCES[1:3]}')
    Result('wall',
           '''
           grey title="" wantaxis=y screenratio=%f screenht=10
           gainpanel=a pclip=99
           grid1=y grid2=y g1num=0.1 g2num=0.1
           ''' % (3*par['ratio']) )


    # wavefield movie frames
    for j in range(0,par['nt']/par['jsnap'],1):
        fdmod.wframe('wa-' +str(j),'wa', j,'pclip=99.9',par)
        fdmod.wframe('we1-'+str(j),'we1',j,'pclip=99.9',par)
        fdmod.wframe('we2-'+str(j),'we2',j,'pclip=99.9',par)
        

