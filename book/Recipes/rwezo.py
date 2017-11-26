try:    from rsf.cluster import *
except: from rsf.proj    import *
import pplot,fdmod

def param(par):
    if(not par.has_key('nsz')): par['nsz']=3
    if(not par.has_key('nsx')): par['nsx']=3

    par['dw']=1/(par['nT']*par['dT'])
#    par['nw']=par['nT']/500*100

    # ------------------------------------------------------------
    par['xmin']=par['ox']                            
    par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    par['zmin']=par['oz']                            
    par['zmax']=par['oz'] + (par['nz']-1) * par['dz']
    par['tmin']=par['ot']
    par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    par['gmin']=par['og']
    par['gmax']=par['og'] + (par['ng']-1) * par['dg']
    # ------------------------------------------------------------

    if(not par.has_key('ratio')):    par['ratio']=1.0*(par['zmax']-par['zmin'])/(par['xmax']-par['xmin'])
    if(not par.has_key('height')):   par['height']=par['ratio']*14
    if(par['height']>10): par['height']=10
    
    if(not par.has_key('ntap')):      par['ntap']=10
    if(not par.has_key('prefix')):    par['prefix']=''
    if(not par.has_key('scalebar')):  par['scalebar']='n'
    if(not par.has_key('labelattr')): par['labelattr']=' labelsz=6 labelfat=3 titlesz=12 titlefat=3 '
    # parallel2=n

def cgrey(custom,par):
    return '''
    window j1=2 j2=2 |
    grey parallel2=n labelrot=n wantaxis=y wanttitle=y
    title="" pclip=100
    min1=%g max1=%g label1=%s unit1=%s 
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],
           custom)

def rgrey(custom,par):
    return '''
    window j1=2 j2=2 |
    grey parallel2=n labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['gmin'],par['gmax'],par['lg'],par['ug'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def dgrey(custom,par):
    return '''
    grey parallel2=n labelrot=n wantaxis=y wanttitle=n
    title="" pclip=99
    min2=%g max2=%g %s
    %s
    ''' % (par['xmin'],par['xmax'],par['lx'],
           custom)

def cgraph(custom,par):
    return '''
    graph labelrot=n
    yreverse=y wantaxis=n title=" " 
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'], 
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],
           custom)

# plot coordinate system
def cos(cos,jray,jwft,custom,par):
    wft = cos + '-wft'
    ray = cos + '-ray'

    par['jray']=jray
    par['jwft']=jwft

    par['jray2']=jray//5
    par['jwft2']=jwft//5

    Plot(ray,cos,'window j1=%(jray)d j2=%(jwft2)d | transp |' % par
         + fdmod.cgraph('plotcol=1 '+custom,par))
    Plot(wft,cos,'window j2=%(jwft)d j1=%(jray2)d |' % par
         + fdmod.cgraph('plotcol=2 '+custom,par))
    
    Plot(cos,[ray,wft],'Overlay')

# slowness in CC and RC
def slow(sloCC,sloRC,vel,cos,par):

    # slowness in CC
    Flow(sloCC,vel,
         '''
         math output="1/input" |
         transp |
         spray axis=2 n=1 o=0 d=1 |
         put label1=x label2=y label3=z
         ''')
    
    # slowness in RC
    Flow(sloRC,[sloCC,cos],
         '''
         window |
         c2r rays=${SOURCES[1]} adj=n |
         put label1=g label2=t label3=
         ''')

def abm(abmRC,abrRC,sloRC,cos,par):
    Flow([abmRC+'tmp',abrRC],[cos,sloRC],
         '''
         rweab naref=1 nbref=1
         slo=${SOURCES[1]} abr=${TARGETS[1]}
         ''')

    Flow(abmRC+'a',abmRC+'tmp','window n3=1 f3=0')
    Flow(abmRC+'b',abmRC+'tmp','window n3=1 f3=1 | smooth rect1=11 rect2=1')
    Flow(abmRC+'m',abmRC+'tmp','window n3=1 f3=2')

    Flow(abmRC,[abmRC+'a',abmRC+'b',abmRC+'m'],
         'cat axis=3 space=n ${SOURCES[1:3]}')

def abmplot(abmRC,par):
    Result(par['prefix']+abmRC+'a',abmRC,'window n3=1 f3=0 | transp |' % par
         + rgrey('pclip=100 bias=1 scalebar=y',par))
    Result(par['prefix']+abmRC+'b',abmRC,'window n3=1 f3=1 | transp |' % par
         + rgrey('pclip=100 allpos=y scalebar=y',par))
    Result(par['prefix']+abmRC+'m',abmRC,'window n3=1 f3=2 | transp |' % par
         + rgrey('pclip=100 allpos=y',par))

def frq(frqRC,frqCC,datCC,cos,par):
    Flow(frqCC,datCC,
         '''
         fft1 opt=n inv=n |
         window squeeze=n min1=%(ow)g n1=%(nw)d |
         transp |
         spray axis=2 n=1 o=0 d=1 | 
         put label1=x label2=y label3=w
         ''' % par)

    Flow(frqRC,[frqCC,cos],
         '''
         window |
         spray axis=2 n=1 o=%(oz)g d=%(dz)g |
         pad beg2=10 n2out=20 |
         c2r rays=${SOURCES[1]} adj=n linear=n nsz=%(nsz)d nsx=%(nsx)d |
         put label1=g label2=t label3=w
         ''' % par)
    Plot(par['prefix']+frqRC,frqRC,'window j3=10 | real | transp |' % par
         + rgrey('title= gainpanel=a',par))

# run migration
def mig(migCC,migRC,frqRC,abmRC,abrRC,cos,par):
    par['j2'] = par.get('j2',1)
    for i in (['SSF','FFD','PSC',
               'F15','F45','F60']):
        sfx = '-' + i
        
        if(i=='F15'): method='method=0 c1=0.50   c2=0.00'
        if(i=='F45'): method='method=0 c1=0.50   c2=0.25'
        if(i=='F60'): method='method=0 c1=0.4761 c2=0.3767'
        if(i=='SSF'): method='method=1'
        if(i=='FFD'): method='method=2 c1=0.50   c2=0.25'
        if(i=='PSC'): method='method=3 c1=0.50   c2=0.25'

        Flow(migRC+sfx,[frqRC,abmRC,abrRC],
             '''
             rwezomig ntap=%d adj=n verb=y %s
             abm=${SOURCES[1]}
             abr=${SOURCES[2]} |
             put label1=g label2=t
             ''' % (par['ntap'],method))

        Flow(migCC+sfx,[migRC+sfx,cos],
             '''
             c2r rays=${SOURCES[1]} adj=y linear=n nsz=%(nsz)d nsx=%(nsx)d
             a2n=%(nz)d a2o=%(oz)g a2d=%(dz)g
             a1n=%(nx)d a1o=%(ox)g a1d=%(dx)g |
             put label1=z label2=x
             ''' % par)

        Plot  (migRC+sfx,'window | transp |'
               + rgrey('pclip=99.9',par))
        Result(par['prefix']+migRC+sfx,migRC+sfx,'window j2=%(j2)d | transp |' % par
               + rgrey('pclip=99.9',par))
        
        Plot(migCC+sfx,'window j1=2 j2=2 | transp |'
             + fdmod.cgrey('pclip=100',par))
        Result(par['prefix']+migCC+sfx,[migCC+sfx,cos],'Overlay' % par)

# run modeling
def mod(modCC,modRC,migRC,abmRC,abrRC,cos,par):
    for i in (['SSF','FFD','PSC',
               'F15','F45','F60']):
        sfx = '-' + i
        
        if(i=='F15'): method='method=0 c1=0.50   c2=0.00'
        if(i=='F45'): method='method=0 c1=0.50   c2=0.25'
        if(i=='F60'): method='method=0 c1=0.4761 c2=0.3767'
        if(i=='SSF'): method='method=1'
        if(i=='FFD'): method='method=2 c1=0.50   c2=0.25'
        if(i=='PSC'): method='method=3 c1=0.50   c2=0.25'
        
        Flow(modRC+sfx,[migRC+sfx,abmRC,abrRC],
             '''
             rwezomig ntap=%d adj=y verb=y %s
             nw=%d ow=%g dw=%g 
             abm=${SOURCES[1]}
             abr=${SOURCES[2]} |
             put label1=g label2=t label3=w
             ''' % (par['ntap'],method,par['nw'],par['ow'],par['dw']) )
        Result(par['prefix']+modRC+sfx,modRC+sfx,'window j3=10 | real | transp |' % par
               + rgrey('title= gainpanel=a',par))

        par['nzcut']=20
        par['ozcut']=-10*par['dz']

        Flow('cutCC'+sfx,[modRC+sfx,cos],
             '''
             c2r rays=${SOURCES[1]} adj=y linear=n nsz=%(nsz)d nsx=%(nsx)d
             a1n=%(nx)d    a1o=%(ox)g    a1d=%(dx)g
             a2n=%(nzcut)d a2o=%(ozcut)g a2d=%(dz)g |
             window n2=1 min2=%(oz)g
             ''' % par)
        
        Flow(modCC+sfx,'cutCC'+sfx,
             '''
             transp |
             pad beg1=%d n1out=%d |
             fft1 opt=n inv=y |
             put label1=t label2=x
             ''' % (int(par['ow']/par['dw']),par['nT']/2+1))
        
        Result(par['prefix']+modCC+sfx,modCC+sfx,
               'grey pclip=100 wanttitle=n grid=y label1=t unit1=s label2=x unit2=km' % par)

# produce plots
def plots(par):
    pplot.p3x2(par['prefix']+'iCCall',
               'migCC-SSF','migCC-PSC','migCC-FFD',
               'migCC-F15','migCC-F45','migCC-F60',
               0.3,0.3,-8,-12)
    pplot.p3x2(par['prefix']+'iRCall',
               'migRC-SSF','migRC-PSC','migRC-FFD',
               'migRC-F15','migRC-F45','migRC-F60',
               0.3,0.3,-10,-12)

    Result(par['prefix']+'icomp',['imgCC','migCC-FFD','migCC-F60'],'Movie')
    Result(par['prefix']+'imgCC',['imgCC','cos'],'Overlay')

    Plot  ('imgCC',
           'window | transp |' + fdmod.cgrey('pclip=100',par))
    Plot  ('imgRC','migCC-FFD',
           '         transp |' + fdmod.cgrey('pclip=100',par))
    Plot('imgRC-ovl',['imgRC','cos'],'Overlay')
    pplot.p2x1(par['prefix']+'CCvsRC','imgRC-ovl','imgCC',0.5,0.5,-9)
    
#    Result('CCvsRC',['imgRC-ovl','imgCC'],'OverUnderIso')

        
