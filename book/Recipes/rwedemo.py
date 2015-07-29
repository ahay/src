try:    from rsf.cluster import *
except: from rsf.proj    import *
import pplot

def cgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y wanttitle=y
    title="" pclip=100 label1="z" unit1=km label2="x" unit2=km %s
    min1=%g max1=%g min2=%g max2=%g screenratio=%g screenht=%g
    ''' % (custom,par['zmin'],par['zmax'],par['xmin'],par['xmax'],par['ratio'],par['height'])
def rgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y wanttitle=y
    title="" pclip=100 label1="t" unit1=s label2="g"
    %s
    min1=%g max1=%g min2=%g max2=%g
    ''' % (custom,par['tmin'],par['tmax'],par['gmin'],par['gmax'])

def cgraph(custom,par):
    return '''
    graph labelrot=n  %s
    yreverse=y wantaxis=n title=" " 
    min1=%g max1=%g min2=%g max2=%g  screenratio=%g screenht=%g
    ''' % (custom,par['xmin'],par['xmax'],par['zmin'],par['zmax'],par['ratio'],par['height'])

def param(par):
    par['dw']=1/(par['nT']*par['dT'])
    par['nw']=par['nT']/500*100

    # ------------------------------------------------------------
    par['xmin']=par['ox']                             -0.1
    par['xmax']=par['ox'] + (par['nx']-1) * par['dx'] +0.1
    par['zmin']=par['oz']                             -0.1
    par['zmax']=par['oz'] + (par['nz']-1) * par['dz'] +0.1
    par['tmin']=par['ot']
    par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    par['gmin']=par['og']
    par['gmax']=par['og'] + (par['ng']-1) * par['dg']
    # ------------------------------------------------------------

    par['ratio']=(par['zmax']- par['zmin'])/(par['xmax']- par['xmin'])
    par['height']=par['ratio']*14
    
# plot coordinate system
def cos(cos,jray,jwft,par):
    wft = cos + '-wft'
    ray = cos + '-ray'

    par['jray']=jray
    par['jwft']=jwft

    Plot(ray,cos,'window j1=%(jray)d | transp |' % par
         + cgraph('plotcol=1',par))
    Plot(wft,cos,'window j2=%(jwft)d |' % par
         + cgraph('plotcol=2',par))
    
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
    Flow([abmRC,abrRC],[cos,sloRC],
         '''
         rweab naref=1 nbref=1
         slo=${SOURCES[1]} abr=${TARGETS[1]}
         ''')

def abmplot(abmRC,par):
    Plot('aRC',abmRC,'window n3=1 f3=0 | transp |'
         + rgrey('pclip=100 bias=1 scalebar=y',par))
    Plot('bRC',abmRC,'window n3=1 f3=1 | transp |'
         + rgrey('pclip=100 allpos=y scalebar=y',par))
    Plot('mRC',abmRC,'window n3=1 f3=2 | transp |'
         + rgrey('pclip=100 allpos=y',par))

    Result(abmRC,['aRC','bRC'],'SideBySideIso')
#    Result(abmRC,['aRC','bRC','mRC'],'OverUnderAniso')


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
         c2r rays=${SOURCES[1]} adj=n linear=n nsz=3 nsx=3 |
         put label1=g label2=t label3=w
         ''' % par)
    Result(frqRC,'window j3=10 | real | transp |'
           + rgrey('title= gainpanel=a',par))

# run migration
def mig(migCC,migRC,frqRC,abmRC,abrRC,cos,par):
    for i in (['SSF','FFD','PSC',
               'F15','F45','F60']):
        sfx = '-' + i
        
        if(i=='F15'): method='method=0 c1=0.50   c2=0.00'
        if(i=='F45'): method='method=0 c1=0.50   c2=0.25'
        if(i=='F60'): method='method=0 c1=0.4761 c2=0.3767'
        if(i=='SSF'): method='method=1'
        if(i=='FFD'): method='method=2 c1=0.50   c2=0.00'
        if(i=='PSC'): method='method=3 c1=0.50   c2=0.00'

        Flow(migRC+sfx,[frqRC,abmRC,abrRC],
             '''
             rwezomig ntap=10 adj=n verb=y %s
             abm=${SOURCES[1]}
             abr=${SOURCES[2]} |
             put label1=g label2=t
             ''' % method)

        Flow(migCC+sfx,[migRC+sfx,cos],
             '''
             c2r rays=${SOURCES[1]} adj=y linear=n nsz=3 nsx=3
             a2n=%(nz)d a2o=%(oz)g a2d=%(dz)g
             a1n=%(nx)d a1o=%(ox)g a1d=%(dx)g |
             put label1=z label2=x
             ''' % par)

        Plot  (migRC+sfx,'window | transp |'
               + rgrey('title=%s pclip=99',par) % i)
        Result(migRC+sfx,'window | transp |'
               + rgrey('title=%s',par) % i)
        
        Plot(migCC+sfx,'window | transp |'
             + cgrey('title=%s pclip=99',par) % i)
        Result(migCC+sfx,[migCC+sfx,cos],'Overlay')

# run modeling
def mod(modCC,modRC,migRC,abmRC,abrRC,cos,par):
    for i in (['SSF','FFD','PSC',
               'F15','F45','F60']):
        sfx = '-' + i
        
        if(i=='F15'): method='method=0 c1=0.50   c2=0.00'
        if(i=='F45'): method='method=0 c1=0.50   c2=0.25'
        if(i=='F60'): method='method=0 c1=0.4761 c2=0.3767'
        if(i=='SSF'): method='method=1'
        if(i=='FFD'): method='method=2 c1=0.50   c2=0.00'
        if(i=='PSC'): method='method=3 c1=0.50   c2=0.00'
        
        Flow(modRC+sfx,[migRC+sfx,abmRC,abrRC],
             '''
             rwezomig ntap=10 adj=y verb=y %s
             nw=%d ow=%g dw=%g 
             abm=${SOURCES[1]}
             abr=${SOURCES[2]} |
             put label1=g label2=t label3=w
             ''' % (method,par['nw'],par['ow'],par['dw']) )
        Result(modRC+sfx,'window j3=10 | real | transp |'
               + rgrey('title= gainpanel=a',par))

        par['nzcut']=20
        par['ozcut']=-10*par['dz']

        Flow('cutCC'+sfx,[modRC+sfx,cos],
             '''
             c2r rays=${SOURCES[1]} adj=y linear=n nsz=3 nsx=3
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
        
        Result(modCC+sfx,
               'grey pclip=100 wanttitle=n grid=y label1=t unit1=s label2=x unit2=km')

# produce plots
def plots(par):
    Result('icomp','imgCC migCC-FFD migCC-F60','Movie')
    Result('imgCC',['imgCC','cos'],'Overlay')

    pplot.p3x2('iCCall',
               'migCC-SSF','migCC-PSC','migCC-FFD',
               'migCC-F15','migCC-F45','migCC-F60',
               0.3,0.3,-8,-12)
    pplot.p3x2('iRCall',
               'migRC-SSF','migRC-PSC','migRC-FFD',
               'migRC-F15','migRC-F45','migRC-F60',
               0.3,0.3,-10,-12)

    Plot  ('imgCC',
           'window | transp |' + cgrey('pclip=99',par))
    Plot  ('imgRC','migCC-FFD',
           '         transp |' + cgrey('pclip=99',par))
    Plot('imgRC-ovl',['imgRC','cos'],'Overlay')

    pplot.p2x1('CCvsRC','imgRC-ovl','imgCC',0.5,0.5,-9)
    
#    Result('CCvsRC',['imgRC-ovl','imgCC'],'OverUnderIso')

        
