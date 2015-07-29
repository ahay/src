from rsf.proj import *
import fdmod,gfield,math

# ------------------------------------------------------------
def param():
    par = {
        'co':3000, # velocity (m/s)
        'fo':200,  # central frequency (Hz)
        'nt':3000,  'ot':0, 'dt':0.0001, 'lt':'Time',     'ut':'s', 'nht':10,
        'nx':1201,  'ox':0, 'dx':1.0,    'lx':'Position', 'ux':'m', 'nhz':3,
        'nz':801,   'oz':0, 'dz':1.0,    'lz':'Depth',    'uz':'m', 'nhx':3,
        'kt':100,  # wavelet delay
        'nb':200, 
        'jsnap':250,
        'nexp':2,
        'tpad':200,
        'tcut':20,
        'labelattr':'labelsz=7 labelfat=2',
        'xlo':16,
        'zlo':14,
        'mintplot':0.1,
        'vmin':1275,
        'vbias':1500
        }

    par['lo']=par['co']/par['fo']

    return par
# ------------------------------------------------------------

# ------------------------------------------------------------
# target coordinates
def target(nxwin,nzwin,qq,par):
    if(not par.has_key('zt')): par['zt'] = par['oz']+0.85*par['nz']*par['dz']
    if(not par.has_key('xt')): par['xt'] = par['ox']+0.50*par['nx']*par['dx']
    
    par['nqx'] = nxwin *     par['lo']/par['dx']
    par['oqx'] = par['xt'] - par['nqx']/2 * par['dx']
    par['dqx'] =                            par['dx']
    par['nqz'] = nzwin *     par['lo']/par['dz']
    par['oqz'] = par['zt'] - par['nqz']/2 * par['dz']
    par['dqz'] =                            par['dz']

    par['nqz']=int(par['nqz'])
    par['nqx']=int(par['nqx'])

    par['wzmin'] = par['oqz']
    par['wzmax'] = par['oqz'] + (par['nqz']-1)*par['dqz']
    par['wxmin'] = par['oqx']
    par['wxmax'] = par['oqx'] + (par['nqx']-1)*par['dqx']
    par['wratio']=(par['wzmax']-par['wzmin'])/(par['wxmax']-par['wxmin'])
    if(not par.has_key('wheight')): par['wheight']=par['height']

    # ------------------------------------------------------------
    fdmod.boxarray(qq,
                   par['nqz'],par['oqz'],par['dqz'],
                   par['nqx'],par['oqx'],par['dqx'],par)
    Plot(qq,'window j2=53 |' + fdmod.qqplot('',par))

# ------------------------------------------------------------
def receivers(rr,par):    
    fdmod.horizontal(rr,par['oz']+par['lo'],par)
    Plot(rr,'window |' + fdmod.rrplot('',par))

# ------------------------------------------------------------
def points(sp,sa,wp,xm,qq,par):
    
    # ------------------------------------------------------------
    ls = 2*par['lo']    

    par['zs' ] = par['oz'] + par['lo']
    par['zs0'] = par['zt']+ls
    par['zs1'] = par['zt']
    par['zs2'] = par['zt']-ls

    par['xs' ] = par['xt']
    par['xs0'] = par['xt']
    par['xs1'] = par['xt']-ls
    par['xs2'] = par['xt']+ls

    par['jzs0'] = par['zs0'] / par['dz']
    par['jzs1'] = par['zs1'] / par['dz']
    par['jzs2'] = par['zs2'] / par['dz']
    
    par['jxs0'] = par['xs0'] / par['dx']
    par['jxs1'] = par['xs1'] / par['dx']
    par['jxs2'] = par['xs2'] / par['dx']
    
    # ------------------------------------------------------------
    # source positions
    fdmod.point('ss0',par['xs0'],par['zs0'],par)
    fdmod.point('ss1',par['xs1'],par['zs1'],par)
    fdmod.point('ss2',par['xs2'],par['zs2'],par)
    Flow(sp,'ss0 ss1 ss2','cat axis=2 space=n ${SOURCES[1:3]} | window n1=2')
    Plot(sp,'window |' + fdmod.ssplot('plotcol=2',par))

    fdmod.point(sa,par['xs'],par['zs'],par)    
    Plot(sa,'window |' + fdmod.ssplot('',par))

    # ------------------------------------------------------------
    # source wavelet
    fdmod.wavelet(wp+'_',par['fo'],par)
    Flow(wp,wp+'_','spray axis=2 n=3 o=0 d=1 | transp')

    # ------------------------------------------------------------
    # scatterrers
    Flow(xm,None,
         '''
         spike nsp=3 mag=1,1,1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=%(jzs0)d,%(jzs1)d,%(jzs2)d l1=%(jzs0)d,%(jzs1)d,%(jzs2)d
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=%(jxs0)d,%(jxs1)d,%(jxs2)d l2=%(jxs0)d,%(jxs1)d,%(jxs2)d
         ''' % par)

    Plot(  xm,'smooth rect1=3 rect2=3 repeat=3 |'
           +fdmod.cgrey('pclip=100',par))
    Result(xm,[xm,qq],'Overlay')

# ------------------------------------------------------------
def point(sp,sa,wp,xm,qq,par):

    # ------------------------------------------------------------
    ls = 2*par['lo']    

    par['zs' ] = par['oz'] + par['lo']
    par['zs0'] = par['zt']

    par['xs' ] = par['xt']
    par['xs0'] = par['xt']

    par['jzs0'] = (par['zs0']-par['oz']) / par['dz']
    par['jxs0'] = (par['xs0']-par['ox']) / par['dx']
    
    # ------------------------------------------------------------
    # source positions
    fdmod.point(sp,par['xs0'],par['zs0'],par)
    Plot(sp,'window |' + fdmod.ssplot('',par))

    fdmod.point(sa,par['xs'],par['zs'],par)    
    Plot(sa,'window |' + fdmod.ssplot('',par))

    # ------------------------------------------------------------
    # source wavelet
    fdmod.wavelet(wp+'_',par['fo'],par)
    Flow(wp,wp+'_','transp')

    # ------------------------------------------------------------
    # scatterrers
    Flow(xm,None,
         '''
         spike nsp=3 mag=1,1,1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=%(jzs0)d l1=%(jzs0)d
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=%(jxs0)d l2=%(jxs0)d
         ''' % par)

    Plot(  xm,'smooth rect1=3 rect2=3 repeat=3 |'
           +fdmod.cgrey('pclip=100',par))
    Result(xm,[xm,qq],'Overlay')

# ------------------------------------------------------------
def segments(sp,sa,wp,xm,qq,par):

    par['ff']=30

    par['zs' ] = par['oz'] + par['lo']
    par['zs0'] = par['zt']

    par['xs' ] = par['ox']+0.35*par['nx']*par['dx']
    par['xs0'] = par['xt']

    par['jzs0'] = 0.85*par['nz']
    par['jxs0'] = 0.50*par['nx']
    
    n = 200
    m = 26
    par['nspk'] = n*m

    par['allx'] = ''
    par['allz'] = ''
    par['allr'] = ''
    par['alli'] = ''
    
    ind = 0
    for j in range(m):
        for i in range(n):
            xc = par['jxs0'] + (i-n/2) - (j-m/2+0.5) * 7*par['lo']/ par['dx']
            zc = par['jzs0'] + (i-n/2) * math.tan(math.radians(par['ff']))
            par['allx'] += '%d,'%xc
            par['allz'] += '%d,'%zc
            par['allr'] += '1,'

            ind +=1
            par['alli'] += '%d,'%ind

    # ------------------------------------------------------------
    # source positions
    Flow(sp+'_x',None,'spike nsp=%(nspk)s mag=%(allx)s n1=%(nspk)s k1=%(alli)s' %par)
    Flow(sp+'_z',None,'spike nsp=%(nspk)s mag=%(allz)s n1=%(nspk)s k1=%(alli)s' %par)
    Flow(sp+'_r',None,'spike nsp=%(nspk)s mag=%(allr)s n1=%(nspk)s k1=%(alli)s' %par)
    Flow(sp,[sp+'_x',sp+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0:2]} | transp
         ''', stdin=0)
    Plot(sp,'window |' + fdmod.ssplot('',par))
    
    fdmod.point(sa,par['xs'],par['zs'],par)
    Plot(sa,'window |' + fdmod.ssplot('',par))

    # ------------------------------------------------------------
    # source wavelet
    fdmod.wavelet(wp+'_',par['fo'],par)
    Flow(wp,wp+'_','transp' % par)
    #    Flow(wp,wp+'_','spray axis=2 n=%(nspk)d o=0 d=1 | transp' % par)

    # ------------------------------------------------------------
    # scatterers
    Flow(xm,None,
         '''
         spike nsp=%(nspk)s mag=%(allr)s
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=%(allz)s
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=%(allx)s
         ''' % par)
    Plot(  xm,'smooth rect1=5 rect2=5 repeat=1 |'
           +fdmod.cgrey('pclip=98',par))
    Result(xm,[xm,qq],'Overlay')

# ------------------------------------------------------------
def reflectivity(refl,sp,sa,wp,xm,qq,par):

    par['zs' ] = par['lo']
    par['zs0'] = par['zt']
    
    par['xs' ] = par['ox']+0.50*par['nx']*par['dx']
    par['xs0'] = par['xt']
    
    par['jzs0'] = 0.85*par['nz']
    par['jxs0'] = 0.50*par['nx']

    Flow(refl+'_mask',refl,
         '''
         window f1=300 |
         mask min=-0.01 max=0.01 |
         dd type=float |
         scale rscale=-1 |
         add add=1 |
         dd type=integer |
         pad beg1=300
         ''')
    Result(refl+'_mask',
           'sfdd type=float |' + fdmod.cgrey('pclip=100',par))

    par['nspk']=par['nz']*par['nx']
    for x in ('x1','x2'):
        if(x=='x1'): o=sp+'_z'
        if(x=='x2'): o=sp+'_x'
        
        Flow(o,[refl,refl+'_mask'],
             '''
             math output=%s |
             put n1=1 n2=%d |
             headerwindow mask=${SOURCES[1]} |
             window
             ''' % (x,par['nspk']))
    Flow(sp,[sp+'_x',sp+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0:2]} | transp
         ''', stdin=0)
    Plot(sp,'window |' + fdmod.ssplot('',par))

    fdmod.point(sa,par['xs'],par['zs'],par)
    Plot(sa,'window |' + fdmod.ssplot('',par))

    # ------------------------------------------------------------
    # source wavelet
    Flow(sp+'_r',[refl,refl+'_mask'],
         '''
         put n1=1 n2=%d |
         headerwindow mask=${SOURCES[1]} |
         window
         '''  % (par['nspk']) )
    Flow(sp+'_mask',sp+'_r',
         '''
         spray axis=2 n=%(nt)d o=%(ot)g d=%(dt)g
         ''' % par)
    
    fdmod.wavelet(wp+'_',par['fo'],par)
    Flow(wp,[wp+'_',sp+'_mask'],
         '''
         spray axis=2 n=32914 o=0 d=1 |
         transp |
         math m=${SOURCES[1]} output="input*m"
         ''' % par)

    # ------------------------------------------------------------
    # scatterers
    Flow(xm,refl,'window')
    Plot(  xm,'smooth rect1=3 rect2=3 repeat=3 |'
           +fdmod.cgrey('pclip=98',par))
    Result(xm,[xm,qq],'Overlay')

# ------------------------------------------------------------
def wavelet(wav,par):

    fdmod.wavelet(wav+'_',par['fo'],par)
    Flow(wav,wav+'_','transp')

# ------------------------------------------------------------
# generate random field
def random(seed,gg,mask,ff,aa,ru,rv,par):

    part=par.copy()

    # random field
    part['ff']=ff # ???
    part['aa']=aa # angle    
    part['ru']=ru # characteristic length
    part['rv']=rv # characteristic length
    gfield.execute(gg+'_',seed,part)
    Flow(gg,[gg+'_',mask],'add mode=p ${SOURCES[1]}')
    Result(gg,fdmod.cgrey('color=F',par))    
    
# ------------------------------------------------------------
def model(vo,vv,rm,gg,gm,par):
    if(not par.has_key('vbias')): par['vbias']=0
    
    Flow(gg+'-nu', gg,'math output="input*%g"' % gm)

    # random velocity
    Flow(vv,[vo,gg+'-nu'],
         '''
         math v=${SOURCES[0]} n=${SOURCES[1]}
         output="v/sqrt(1+n)"
         ''' % par)
    
    Plot(vo,fdmod.cgrey('allpos=y bias=1.45',par))
    Plot(vv,fdmod.cgrey('allpos=y bias=1.45',par))
    Result(vo,vo,'Overlay')
    Result(vv,vv,'Overlay')
    
    # density
    Flow(rm,vo,'math output=1')
    Plot(rm,fdmod.cgrey('allpos=y pclip=100',par))

# ------------------------------------------------------------
# traveltimes
def ttimes(vel,sou,par):
    Flow(vel+'-hwt',vel,
         '''
         smooth rect1=20 rect2=20 |
         hwt2d verb=n xsou=%g zsou=%g
         nt=%d ot=%g dt=%g
         ng=%d og=%g dg=%g
         ''' % (par['xs0'],par['zs0'],
                par['nt']/10,par['ot'],par['dt']*10,
                1441,90,.125) )

    fdmod.rayplot(vel+'-hwt',20,20,5,20,'',par)
    Result(vel+'-hwt',[vel,vel+'-hwt',sou],'Overlay')

    Flow(vel+'-fme',vel,
         '''
         smooth rect1=20 rect2=20 |
         eikonal zshot=%g yshot=%g
         ''' %(par['zs0'],par['xs0']) )    
    Plot(vel+'-fme',fdmod.ccont('plotcol=2 dc=0.1',par))
    Result(vel+'-fme',[vel,vel+'-fme',sou],'Overlay')

# ------------------------------------------------------------
def beam(vel,hwt,xsou,zsou,sou,par):

    Flow(hwt,vel,
         '''
         smooth rect1=20 rect2=20 |
         hwt2d verb=n xsou=%g zsou=%g
         nt=%d ot=%g dt=%g
         ng=%d og=%g dg=%g
         ''' % (
        xsou,zsou,
        par['nt'],par['ot'],par['dt'],
        21,-25,.125) )

    Plot(hwt+'_ray',hwt,'transp | window f1=00 j1=25 n2=10 |'
         + fdmod.cgraph('plotcol=4 plotfat=5',par))

    Result(hwt,[vel,hwt+'_ray',sou],'Overlay')

# ------------------------------------------------------------
# passive array modeling
def pmodel(dat,wfl,wav,vel,den,sou,rec,ico,par):
    if(not par.has_key('mintplot')): par['mintplot']=0

    Result(wav,'window n1=1 n2=400 |' + fdmod.waveplot('',par))
    Result('p'+vel,[vel,ico,rec,sou],'Overlay')
    
    fdmod.awefd(dat,wfl,
                wav,vel,den,sou,rec,'',par)

#    Result(dat, 'window j2=5 | transp |' +
#           fdmod.dgrey('pclip=99.9 min1=%(mintplot)g labelsz=4'%par, par))
    Result(dat, 'window j2=5 | transp |' +
           fdmod.dgrey('pclip=99.5 screenratio=1.5 min1=%(mintplot)g labelsz=3 labelfat=3 '%par, par))

# ------------------------------------------------------------
def movie(wfl,vel,par):
    if(not par.has_key('vbias')): par['vbias']=1.5

    fdmod.wom(wfl+'m',wfl,vel,par['vbias'],par)
    Plot(   wfl+'m',fdmod.wgrey('pclip=99',par))

    for i in range(0,par['nt']/par['jsnap'],1):
        fdmod.wframe(wfl+'-'+str(i),wfl+'m',i,'pclip=99',par)
        Result(wfl+'-'+str(i),'Overlay')
        
# ------------------------------------------------------------
# active array modeling
def amodel(dat,wfl,wav,vel,den,ref,sou,rec,ico,par):
    if(not par.has_key('wscale')): par['wscale']=5
    if(not par.has_key('vbias')): par['vbias']=1500

    Result(wav,'window n1=1 n2=400 |' + fdmod.waveplot('format1=%3.2f',par))
    Result('a'+vel,[vel,ico,rec,sou],'Overlay')

    fdmod.lwefd1(dat+'o',wfl+'o',
                 dat+'x',wfl+'x',
                 wav,vel,den,ref,sou,rec,'',par)
    Flow(dat,[dat+'o',dat+'x'],'add ${SOURCES[1]} scale=1,%(wscale)g' % par)
    Flow(wfl,[wfl+'o',wfl+'x'],'add ${SOURCES[1]} scale=1,%(wscale)g' % par)

    Result(dat, 'window j2=5 | transp |' +
           fdmod.dgrey('pclip=99.5 screenratio=1.5 min1=%(mintplot)g labelsz=3 labelfat=3 '%par, par))

# ------------------------------------------------------------
#def wigner(din,dou,nht,nhx,par):
#
#    Flow(dou,din,
#         'wdf verb=y nh1=15 nh2=15 nh3=0')
#    Result(dou, 'window j2=5 | transp |' +
#           fdmod.dgrey('pclip=99.9 screenratio=1.5 min1=%(mintplot)g labelsz=4'%par, par))
    
# ------------------------------------------------------------
# passive array imaging
def pimage(cic,iic,
           dat,vel,den,rec,ico,par):

    # reverse data
    Flow(dat+'-rev',dat,'reverse which=2 opt=i verb=y | pad end2=%(tpad)d' % par)

    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    
    # backpropagate
    fdmod.awefd(dat+'-'+vel+'-bck',
                dat+'-'+vel+'-wfl',
                dat+'-rev',
                vel,den,
                rec,rec,iwindow,par)

    # cut wavefield around t=0
    Flow(dat+'-'+vel+'-cut',
         dat+'-'+vel+'-wfl',
         '''
         window n3=%d f3=%g |
         reverse which=4 opt=i verb=y |
         put o3=%g label3=t unit3=s
         ''' % (2*par['tcut']+1,
                par['nt']-par['kt']-par['tcut'],
                -par['tcut']*par['dt']))
         
    # compute WDF
    Flow(dat+'-'+vel+'-wig',dat+'-'+vel+'-cut',
         'wdf verb=y nh1=%(nhz)d nh2=%(nhx)d nh3=%(nht)d' % par)
    #    Result(dat+'-'+vel+'-wig','grey gainpanel=a')

    # imaging condition
    Flow(cic,dat+'-'+vel+'-cut',
         'window n3=1 f3=%d' % par['tcut'])
    Flow(iic,dat+'-'+vel+'-wig',
         'window n3=1 f3=%d' % par['tcut'])

    for img in ([cic,iic]):
        Plot(img,fdmod.cgrey('pclip=100',par))
        Result( img,[img,rec],'Overlay')
        Result('win'+img,img,
               fdmod.cgrey('min1=%g max1=%g min2=%g max2=%g screenratio=%g wantaxis=y' %
                           (par['wzmin'],par['wzmax'],par['wxmin'],par['wxmax'],par['wratio']),par) )

# ------------------------------------------------------------
# active array imaging
def aimage(cic,iic,
           dat,wav,vel,den,sou,rec,ico,par):

    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '

    # ------------------------------------------------------------
    # source wavefield
    fdmod.awefd(dat+'-'+vel+'-for',
                dat+'-'+vel+'-sou',
                wav,
                vel,den,sou,rec,iwindow,par)
    
    # ------------------------------------------------------------
    # receiver wavefield
    Flow(dat+'-rev',dat,'reverse which=2 opt=i verb=y' % par)
    fdmod.awefd(dat+'-'+vel+'-bck',
                dat+'-'+vel+'-rwf',
                dat+'-rev',
                vel,den,
                rec,rec,iwindow,par)

    Flow(dat+'-'+vel+'-rec',
         dat+'-'+vel+'-rwf',
         '''
         reverse which=4 opt=i verb=y
         ''' )

    # compute assymptotic Wigner distribution
    Flow(dat+'-'+vel+'-wig',
         dat+'-'+vel+'-rec',
         'wdf verb=y nh1=%(nhz)d nh2=%(nhx)d nh3=%(nht)d' % par)
    
    # ------------------------------------------------------------
    # imaging condition
    Flow(cic,[dat+'-'+vel+'-sou',dat+'-'+vel+'-rec'],
         'xcor2d uu=${SOURCES[1]} verb=y nbuf=100 axis=3')
    Flow(iic,[dat+'-'+vel+'-sou',dat+'-'+vel+'-wig'],
         'xcor2d uu=${SOURCES[1]} verb=y nbuf=100 axis=3')

    # ------------------------------------------------------------
    # WDF on image
    Flow(cic+'-wdf',cic,
         'scale axis=123 | wdf verb=y nh1=%(nhz)d nh2=%(nhx)d' % par)
    
    for img in ([cic,iic,cic+'-wdf']):
        Plot(img,fdmod.cgrey('pclip=100',par))
        Result( img,[img,rec,sou],'Overlay')
        Result('win'+img,img,
               fdmod.cgrey('pclip=100 min1=%g max1=%g min2=%g max2=%g screenratio=%g screenht=%g wantaxis=y' %
               (par['wzmin'],par['wzmax'],par['wxmin'],par['wxmax'],par['wratio'],par['wheight']),par))

# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------

# different realizations, same parameters
def realization(gg,mask,ff,aa,ru,rv,
                vo,vv,rm,gm,
                dat,wfl,wav,sou,rec,ico,
                cic,iic,
                par):

    for k in range(0,3):
        ktag = "-r%01d" % k

        # generate random velocity model
        random(112009+k,gg+ktag,mask,ff,aa,ru,rv,par)
        model(vo,vv+ktag,rm+ktag,gg+ktag,gm,par)

        # passive-array modeling and migration
        pmodel(dat+ktag,wfl+ktag,wav,vv+ktag,rm+ktag,sou,rec,ico,par)
        pimage(cic+ktag,iic+ktag,dat+ktag,vo,rm+ktag,rec,ico,par)

# ------------------------------------------------------------
# same realization, different magnitude
def magnitude(gg,
              vo,vv,rm,gm,
              dat,wfl,wav,sou,rec,ico,
              cic,iic,
              par):
    
    for k in range(0,3):
        ktag = "-m%01d" % k
        gmk = gm * (k+1)

        # generate random velocity model
        Flow(gg+ktag,gg,'window')
        model(vo,vv+ktag,rm+ktag,gg+ktag,gmk,par)

        # passive-array modeling and migration
        pmodel(dat+ktag,wfl+ktag,wav,vv+ktag,rm+ktag,sou,rec,ico,par)
        pimage(cic+ktag,iic+ktag,dat+ktag,vo,rm+ktag,rec,ico,par)

# ------------------------------------------------------------
# same magnitude, different anomaly size
def dimension(gg,mask,ff,aa,ru,rv,
              vo,vv,rm,gm,
              dat,wfl,wav,sou,rec,ico,
              cic,iic,
              par):

    for k in range(0,3):
        ktag = "-d%01d" % k

        ruk = ru * (k+1)
        rvk = rv * (k+1)
        
        # generate random velocity model
        random(gg+ktag,mask,ff,aa,ruk,rvk,par)
        model(vo,vv+ktag,rm+ktag,gg+ktag,gm,par)

        # passive-array modeling and migration
        pmodel(dat+ktag,wfl+ktag,wav,vv+ktag,rm+ktag,sou,rec,ico,par)
        pimage(cic+ktag,iic+ktag,dat+ktag,vo,rm+ktag,rec,ico,par)

# ------------------------------------------------------------
def sampling(vo,rm,
             dat,rec,ico,
             cic,iic,
             par):

    for k in range(0,3):
        ktag = "-s%01d" % k
        knum = pow(2,k+1)

        Flow(rec+ktag,rec,'window j2=%d' % (knum))
        Plot(rec+ktag,'window |' + fdmod.rrplot('',par))
        
        Flow(dat+ktag,dat,'window j1=%d' % (knum))
        Result(dat+ktag, 'window j2=5 | transp |' +
               fdmod.dgrey('pclip=99 screenratio=1.5 min1=%(mintplot)g'%par, par))
        
        pimage(cic+ktag,iic+ktag,dat+ktag,vo,rm,rec+ktag,ico,par)

# ------------------------------------------------------------
def ugrey(custom,par):
    return '''
    grey3 flat=n frame1=%d frame2=%d frame3=%d title=""
    label1=z label2=x label3=t
    unit1=m  unit2=m  unit3=s
    point1=0.75 point2=0.75 screenratio=1
    labelsz=6 labelfat=2 wantaxis=y framelabel=n
    %s
    ''' % (par['nqz']/2,par['nqx']/2,par['tcut'],par['labelattr']+custom)

# ------------------------------------------------------------
def wdfic(cii,
          uxx,uyy,
          wxx,wyy,
          dat,vel,den,rec,ico,par):
    
    # ------------------------------------------------------------
    par['jrec']=5
    par['nrec']=320
    par['orec']=0
    par['jrec']=10
    par['nrec']=160
    
    receivers = range(par['orec'],par['orec']+par['nrec']*par['jrec'],par['jrec'])
    # ------------------------------------------------------------

    for k in receivers:
        ktag = "-%04d" % k
        
        # source coordinates
        Flow(rec+ktag,rec,'window n2=1 f2=%g' % k )
        Plot(rec+ktag,    'window |' + fdmod.rrplot('',par))
        
        # velocity (overlay)
        Plot(vel+ktag,[vel,rec+ktag,ico],'Overlay')
    allvxx = map(lambda x: vel+'-%04d' % x,receivers)
    Plot('allvxx',allvxx,'Movie')

    # ------------------------------------------------------------ 
    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # image all traces at once

    # wavefield
    # z-x-t
    Flow(uyy,
         dat+'-'+vel+'-wfl',
         '''
         window n3=%d f3=%g |
         reverse which=4 opt=i verb=y |
         put o3=%g label3=t unit3=s
         ''' % (2*par['tcut']+1,
                par['nt']-par['kt']-par['tcut'],
                -par['tcut']*par['dt']))
    Result(uyy,'byte gainpanel=a pclip=100 |' + igrey('',par))

    # WDF over y
    # z-x-t
    Flow(wyy,
         uyy,
         'wdf verb=y nh1=%(nhz)d nh2=%(nhx)d nh3=%(nht)d' % par)
    Result(wyy,'byte gainpanel=a pclip=100 |' + igrey('',par))
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # image each trace separately
    for k in receivers:
        ktag = "-%04d" % k
        
        # data trace
        Flow(dat+ktag,
             dat,'window squeeze=n n1=1 f1=%g' % k )
        Flow(dat+ktag+'-rev',
             dat+ktag,'reverse which=2 opt=i verb=y | pad end2=%(tpad)d' % par)
        
        fdmod.awefd(dat+ktag+'-bck',
                    dat+ktag+'-wfl',
                    dat+ktag+'-rev',
                    vel,den,
                    rec+ktag,rec+ktag,iwindow,par)

        # z-x-t
        Flow(uyy+ktag,
             dat+ktag+'-wfl',
             '''
             window n3=%d f3=%g |
             reverse which=4 opt=i verb=y |
             put o3=%g label3=t unit3=s
             ''' % (2*par['tcut']+1,
                    par['nt']-par['kt']-par['tcut'],
                    -par['tcut']*par['dt']))
        Result(uyy+ktag,'byte gainpanel=a pclip=100 |' + igrey('',par))

        # z*x-.-t
        Flow(uxx+ktag,
             uyy+ktag,
             'put n1=%d n2=1' % (par['nqz']*par['nqx']) )
        
    # collect traces at the image point for all receiver locations
    # wavefield
    # z*x-xs-t
    Flow(uxx+'-all',
         map(lambda x: uxx+'-%04d' % x,receivers),
         'cat axis=2 space=n ${SOURCES[1:%d]}' % (par['nrec']))

    # WDF over x
    # z*x-xs-t
    Flow(wxx+'-all',
         uxx+'-all',
         '''
         wdf verb=y nh1=0 nh2=100 nh3=%(nht)d
         ''' % par)

    for k in (uxx,wxx):

        # all traces around the target
        # z-x-t
        Flow(k+'-cub',
             k+'-all',
             '''
             stack |
             transp plane=23 |
             put
             n1=%d o1=%g d1=%g label1=z unit1=''
             n2=%d o2=%g d2=%g label2=x unit2=''
             ''' % (par['nqz'],par['oqz'],par['dqz'],
                    par['nqx'],par['oqx'],par['dqx']))

        # one trace at the target
        # xs-t
        Flow(k,
             k+'-all',
             'window n1=1 f1=%d | transp' % (par['nqz']*par['nqx']/2) )        
        
        Result(k,
               '''
               put o2=%g d2=%g | 
               grey title="" pclip=98 labelsz=6 labelfat=2
               label1=%s unit1=%s
               label2=%s unit2=%s
               screenratio=0.3 screenht=4
               %s
               ''' % (par['ox']+par['orec']*par['dx'],
                      par['jrec']*par['dx'],
                      par['lt'],par['ut'],
                      par['lx'],par['ux'],
                      par['labelattr']) )

    # z-x
    Flow(cii,
         wxx+'-cub',
         'window n3=1 f3=%d' % par['tcut'])
    Plot(cii,fdmod.cgrey('pclip=100',par))
    Result(cii,[cii,'rr'],'Overlay')        
    Result('win'+cii,
           cii,
           fdmod.cgrey('min1=%g max1=%g min2=%g max2=%g screenratio=%g wantaxis=y' %
                       (par['wzmin'],par['wzmax'],par['wxmin'],par['wxmax'],par['wratio']),par))
    
    
# ------------------------------------------------------------
def igrey(custom,par):
    return '''
    window n1=%d min1=7.25 n2=%d min2=12.9 |
    grey3 flat=n frame1=%d frame2=%d frame3=%d title=""
    label1=z label2=x label3=t
    unit1=km unit2=km unit3=s
    point1=0.75 point2=0.75 screenratio=1
    labelsz=6 labelfat=2 wantaxis=y framelabel=n
    %s
    ''' % (par['nqz'],par['nqx'],par['nqz']/2,par['nqx']/2,par['tcut'],par['labelattr']+custom)

