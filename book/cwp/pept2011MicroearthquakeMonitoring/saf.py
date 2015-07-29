from rsf.proj import *
import fdmod, gfield, pplot

def wingrey(custom,par):
    return '''
    grey title="" pclip=99.99 wheretitle=t titlesz=20 titlefat=2
    wantaxis=n screenratio=1 screenht=8 min1=%g max1=%g min2=%g max2=%g
    %s
    ''' %(par['zminq'],par['zmaxq'],par['xminq'],par['xmaxq'],custom)

def winss(custom,par):
    return fdmod.ssplot(
        '''
        wantaxis=n screenratio=1 screenht=8
        min2=%g max2=%g min1=%g max1=%g
        ''' % (par['zminq'],par['zmaxq'],par['xminq'],par['xmaxq']),
        par)

def radgrey(custom,par):
    return '''
    grey title="" pclip=100
    labelsz=5 labelfat=3
    screenratio=1 screenht=10
    label2=%s  unit2=%s
    label1="p" unit1=""
    %s
    ''' %(par['lx'],par['ux'],custom)

def souplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=o plotcol=2 plotfat=5 wantaxis=n %s' % custom,par)

# ------------------------------------------------------------
def param(xtarget,ztarget):
    par = {
        'nt':5000, 'ot':0,  'dt':0.0004, 'lt':'Time',    'ut':'s',
        'nx':1251, 'ox':-2, 'dx':0.004,  'lx':'Position','ux':'km',
        'nz':801,  'oz':-1, 'dz':0.004,  'lz':'Depth',   'uz':'km',
        'kt':200,
        'nb':100,
        'iweight':2,
        'wweight':20,
        'jsnap':100
        }
    fdmod.param(par)

    par['xsou']=xtarget
    par['zsou']=ztarget

    par['nzq']=240
    par['ozq']=par['zsou']-par['nzq']/2*par['dz']
    par['dzq']=par['dz']

    par['nxq']=240
    par['oxq']=par['xsou']-par['nxq']/2*par['dx']
    par['dxq']=par['dx']
    
    par['xminq']=par['oxq']
    par['xmaxq']=par['oxq'] + (par['nxq']-1) * par['dxq']
    par['zminq']=par['ozq']
    par['zmaxq']=par['ozq'] + (par['nzq']-1) * par['dzq']

    par['ff']=0
    par['aa']=2.0
    par['ru']=0.05
    par['rv']=0.05

    return par

# ------------------------------------------------------------
def iparam(par):

    pari=par.copy()
    pari['tpad']=110
    pari['twin']=110
    
    pari['nht']=2
    pari['nhx']=5
    pari['nhz']=5
    
    pari['nt']=par['nt']+pari['tpad']

    pari['nqz']=par['nzq']
    pari['oqz']=par['ozq']
    pari['dqz']=par['dzq']

    pari['nqx']=par['nxq']
    pari['oqx']=par['oxq']
    pari['dqx']=par['dxq']

    return pari

# ------------------------------------------------------------
def onesou(ss,qq,wav,par):
    
    fdmod.point(ss,par['xsou'],par['zsou'],par)

    fdmod.wavelet(wav+'_',50,par) 
    Flow(wav,wav+'_','transp')
    
# ------------------------------------------------------------
def thrsou(ss,wav,par):

    par['xsou0']=par['xsou']-0.25
    par['zsou0']=par['zsou']-0.25

    par['xsou1']=par['xsou']
    par['zsou1']=par['zsou']
    
    par['xsou2']=par['xsou']+0.25
    par['zsou2']=par['zsou']+0.25
    
    for j in range(3):
        fdmod.point(ss+str(j),par['xsou'+str(j)],par['zsou'+str(j)],par)
        Flow(ss,[ss+'0',ss+'1',ss+'2'],'cat axis=2 space=n ${SOURCES[1:3]}')

    fdmod.wavelet(wav+'0_',50,par)
    Flow(wav+'0',
         wav+'0_',
         'scale rscale=1.0 | window f1=100 | pad end1=100 | put o1=%(ot)g n2=1' %par)
        
    fdmod.wavelet(wav+'1_',50,par) 
    Flow(wav+'1',
         wav+'1_',
         'window' %par)
    
    fdmod.wavelet(wav+'2_',50,par)
    Flow(wav+'2',
         wav+'2_',
         'scale rscale=1.0 | pad beg1=100 | window n1=%(nt)d | put o1=%(ot)g n2=1' %par)
    
    Flow(wav,[wav+'0',wav+'1',wav+'2'],
         '''
         cat axis=2 space=n ${SOURCES[1:3]} |
         transp
         ''')

# ------------------------------------------------------------
def verticalbh(rr,xrec,par):
    fdmod.vertical(rr,xrec,par)

def safodbh(rr,par):

    rgx_su = 'safod_mh_geom_corr_x_dswp.su'
    Fetch(rgx_su,'safod')
    Flow(rr+'x',rgx_su,'suread')

    rgz_su = 'safod_mh_geom_corr_z_dswp.su'
    Fetch(rgz_su,'safod')
    Flow(rr+'z',rgz_su,'suread')
    Flow(rr,[rr+'x',rr+'z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} |
         transp |
         put d1=1 d2=1 |
         scale rscale=0.001
         ''',stdin=0)

# ------------------------------------------------------------
def safodmodel(vel,den,par):

    # velocity
    vfile_su = 'safod_velmod_transp_dswp.su'
    Fetch(vfile_su,'safod')
    Flow(vel+'_',vfile_su,'suread')

    # density
    dfile_su = 'safod_denmod_transp_sm_dswp.su'
    Fetch(dfile_su,'safod')
    Flow(den+'_',dfile_su,'suread')
    
    for i in ([vel,den]):
        Flow(i,i+'_',
             '''
             put
             o1=-1.000 d1=0.025 label1=z unit1=km
             o2=-1.992 d2=0.025 label2=x unit2=km |
             remap1 o1=%(oz)g n1=%(nz)d d1=%(dz)g |
             transp |
             remap1 o1=%(ox)g n1=%(nx)d d1=%(dx)g |
             transp |
             scale rscale=0.001
             ''' % par)

def geometry(ss,rr,qq,par):
    # source coordinates
    Plot(ss,fdmod.ssplot('',par))
    Plot('win'+ss,ss,winss('',par))
    
    # receiver coordinates
    Plot(rr,fdmod.rrplot('plotfat=12',par))
    
    # image coordinates
    fdmod.makebox(qq,par['zminq'],par['zmaxq'],par['xminq'],par['xmaxq'],par)
    Plot(qq,fdmod.bbplot('plotcol=1',par))
    
# ------------------------------------------------------------
def modeling(wav,vel,den,ss,rr,qq,dat,wfl,par):
    fdmod.awefd(dat+'-tmp',
                wfl,
                wav,
                vel,
                den,
                ss,rr,'free=y',par)

    Flow(dat,
         dat+'-tmp',
         'window f2=%(kt)d | pad end2=%(kt)d | put o2=%(ot)g' %par)

    # wavefield movie
    fdmod.wom(wfl+'om',wfl,vel,4.0,par)
    Plot(     wfl+'om',
              fdmod.wgrey('pclip=99',par),view=1)
    
    for f in range(0,par['nt']/par['jsnap'],2):
        tag = "-%02d"%f
        fdmod.wframe(   wfl+'om-'+tag,wfl+'om',f,'pclip=99',par)
        Result(wfl+tag,[wfl+'om-'+tag,ss,rr,qq],'Overlay')

# ------------------------------------------------------------
def migration(dat,vel,den,ss,rr,qq,uuu,wig,cic,iic,par):

    Result(dat,
           '''
           reverse which=2 opt=i | window f2=%(tpad)d | 
           put d1=1 o2=%(ot)g |
           grey title=""
           label1="Receiver #" unit1="" label2="Time" unit2="s"
           screenratio=0.3 screenht=4
           labelsz=5 labelfat=3 %(labelattr)s
           '''% par)

    fdmod.awefd(uuu+'-out',
                uuu+'-bck',
                dat,
                vel,
                den,
                ss,
                rr,
                'free=y jsnap=1' + fdmod.qqbox2d(par),par)

    Flow(uuu,
         uuu+'-bck',
         '''
         window n3=%d f3=%d |
         put o3=%g label3=t unit3=s
         ''' % (2*par['twin'],
                par['nt']-par['tpad']-par['twin'],
                -par['twin']*par['dt']))

    # CIC
    Flow(cic,uuu,
         'window n3=1 f3=%(twin)d' % par)
    
    # IIC
    Flow(wig,uuu,
         'wigner verb=y nh1=%(nhz)d nh2=%(nhx)d nh3=%(nht)d' %par)
    Flow(iic,wig,
         'window n3=1 f3=%(twin)d' %par)

    for k in ([uuu,wig]):
        wflplot(k,'winss',par)

    for k in ([cic,iic]):

        # image
        Plot(  k,fdmod.cgrey('pclip=99.99',par))
        Result(k,[k,ss],'Overlay')
        
        # image (window)
        Plot(  'win'+k,k,wingrey('pclip=99.99',par))
        Result('win'+k,['win'+k,'winss'],'Overlay')
        
        # slope decomposition
        slope(k+'-ssk',k,par)

# ------------------------------------------------------------
def common(wav,vel,den,ss,rr,qq,par):
    
    # ------------------------------------------------------------
    # wavelet
    Result( 'wav','window n2=1000 | transp |' + fdmod.waveplot('',par))

    # ------------------------------------------------------------
    # plot model geometry
    geometry(ss,rr,qq,par)
    
    # ------------------------------------------------------------
    # velocity and density
    safodmodel(vel,den,par)
    Plot(vel,fdmod.cgrey('bias=1.9 allpos=y ',par))
    Plot(den,fdmod.cgrey('bias=1.72',par))
    Result(vel,[vel,ss,rr,qq],'Overlay')

# ------------------------------------------------------------
# random fluctuations
def random(gg,par):
    gfield.execute(gg,3.1415,par)
    Result(gg,fdmod.cgrey('color=j',par))

# random velocity model
def randomvel(vsm,vrn,gg,scale,ss,rr,qq,par):
    Flow(vrn+'-gg',gg,
         '''
         window n1=%d n2=%d |
         math "output=input*%g"
         ''' % (par['nz'],par['nx'],scale))
    
    Flow(vrn,[vsm,vrn+'-gg'],
         '''
         math n=${SOURCES[1]}
         output="input/sqrt(1+n)"
         ''' % par)

    Plot(vrn,fdmod.cgrey('bias=1.9 allpos=y',par))
    Result(vrn,[vrn,ss,rr,qq],'Overlay')
        
# ------------------------------------------------------------
# wavefield plots
def wflplot(wfl,ss,par):

    Flow([wfl+'-plt',wfl+'-bar'],wfl,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=99.99')
        
    for j in range(9):
        tag = "-%02d" %j

        jdt = 25
        fdt = par['twin']-4*jdt
        time = "%2g ms" % ((j-4)*jdt*par['dt']*1000)

        Plot(wfl+tag+'_',
             [wfl+'-plt',wfl+'-bar'],
             'window f3=%d j3=%d | window n3=1 f3=%d |' % (fdt,jdt,j)
             + wingrey('title="%s" wheretitle=t titlesz=20 titlefat=2 bar=${SOURCES[1]}' % time,par)
             )
        Plot(wfl+tag,[wfl+tag+'_',ss],'Overlay')
        
    pplot.p1x5(wfl,
               wfl+'-00',wfl+'-02',
               wfl+'-04',
               wfl+'-06',wfl+'-08',0.40,0.40,-6.25)
    
# ------------------------------------------------------------
# slope decomposition
def slope(ssk,img,par):
    par['np']=400
    par['p0']=-2
    par['dp']=0.01
    
    Flow(img+'_',
         img,
         '''
         transp |
         put o2=%g
         ''' % (-par['nzq']/2*par['dzq']))
    Flow(ssk,
         img+'_',
         '''
         slant adj=y np=%(np)d p0=%(p0)g dp=%(dp)g verb=y |
         transp
         ''' %par)

    Result(ssk,radgrey('',par))

# ------------------------------------------------------------
# focusing measure
def focus(foc,uuu,ss,par):

    par['fdims']=3
    par['frec1']=1
    par['frec2']=1
    par['frec3']=10

    Flow(foc+'-raw',
         uuu,
         '''
         scale axis=123 |
         focus dim=%(fdims)d rect1=%(frec1)d rect2=%(frec2)d rect3=%(frec3)d niter=20
         ''' % par)

    Flow(foc,
         foc+'-raw',
         'math output="abs(input)"')

    wflplot(foc,ss,par)
    
# ------------------------------------------------------------
def isparse(ss,wav,rr,n,par):
    pari = iparam(par)
    
    # model and geometry
    common('wav','vsm','den',ss,rr,'qq',par)

    # ------------------------------------------------------------
    # FDmodeling
    modeling('wav',
             'vsm','den',
             ss,rr,'qq',
             'dat','wfl',par)

    Flow('rev',
         'dat',
         'reverse which=2 opt=i verb=y | pad end2=%(tpad)d' % pari)
    
    # ------------------------------------------------------------
    for i in range(n):
        j = int(pow(2,i))
        l = "%01d" %i

        Flow(rr+l,rr,'window j2=%d' %j)
        Plot(rr+l,fdmod.rrplot('plotfat=12',par))
        Result('vsm'+l,['vsm',ss,rr+l,'qq'],'Overlay')

        for f in range(0,par['nt']/par['jsnap'],2):
            tag = "-%02d"%f
            Result('wfl'+l+tag,['wflom-'+tag,ss,rr+l,'qq'],'Overlay')
       
        Flow('rev'+l,'rev','window j1=%d' %j)

        # FDmigration      
        migration('rev'+l,
                  'vsm','den',
                  rr+l,rr,'qq',
                  'uuu'+l,'wig'+l,
                  'cic'+l,'iic'+l,pari)

# ------------------------------------------------------------
def irandom(ss,wav,rr,n,par):
    pari = iparam(par)
    
    # model and geometry
    common('wav','vsm','den',ss,rr,'qq',par)
    random('gg',par)

    for i in range(n):
        l = "%01d" %i

        vscale = i * 0.10
        randomvel('vsm','vrn'+l,'gg',vscale,ss,rr,'qq',par)

        # FDmodeling
        modeling('wav',
                 'vrn'+l,'den',
                 ss,rr,'qq',
                 'dat'+l,'wfl'+l,par)

        Flow('rev'+l,
             'dat'+l,
             'reverse which=2 opt=i verb=y | pad end2=%(tpad)d' % pari)

        migration('rev'+l,
                  'vsm','den',
                  rr,rr,'qq',
                  'uuu'+l,'wig'+l,
                  'cic'+l,'iic'+l,pari)

#        focus('fcic'+l,'uuu'+l,'win'+ss,pari)
#        focus('fiic'+l,'wig'+l,'win'+ss,pari)
        
# ------------------------------------------------------------
def imixt(ss,wav,rr,n,par):
    pari = iparam(par)
    
    # model and geometry
    common('wav','vsm','den',ss,rr,'qq',par)
    random('gg',par)

    randomvel('vsm','vrn','gg',0.30,ss,rr,'qq',par)

    # FDmodeling
    modeling('wav',
             'vrn','den',
             ss,rr,'qq',
             'dat','wfl',par)
    
    Flow('rev',
         'dat',
         'reverse which=2 opt=i verb=y | pad end2=%(tpad)d' % pari)
    
    i = n-1
    j = int(pow(2,i))
    l = "%01d" %i
    
    Flow(rr+l,rr,'window j2=%d' %j)
    Plot(rr+l,fdmod.rrplot('plotfat=12',par))
    Result('vsm'+l,['vsm',ss,rr+l,'qq'],'Overlay')
    
    Flow('rev'+l,'rev','window j1=%d' %j)
    
    # FDmigration      
    migration('rev'+l,
              'vsm','den',
              rr+l,rr,'qq',
              'uuu'+l,'wig'+l,
              'cic'+l,'iic'+l,pari)
    
