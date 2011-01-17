try:
    from rsf.cluster import *
except:
    from rsf.proj import *
import sigsbee,fdmod,spmig,encode,adcig,zomig

# ------------------------------------------------------------
def param(par):

    # anomaly magnitude (fraction)
    par['globalpct']=0.10
    par['localpct']=0.20

    # model params
    par['jximg']=3
    par['dx']=par['dx']*par['jximg']
    par['nx']=1024
    par['nz']=par['nz']+par['nzpad']

    # migration params
    par['verb']='y'
    par['eps']=0.1
    par['nrmax']=3
    par['dtmax']=0.00005
    par['tmx']=16

    # shot params
    par['fs']=20 # shots first
    par['ns']=40 # shots number
    par['js']=10 # shots jump
    
    par['jds']=par['ds']*par['js']
    par['jos']=par['os']+par['fs']*par['ds']

    # I.C. params
    par['nht']=120
    par['oht']=-0.60
    par['dht']=0.01
    par['misc']='itype=t nht=%(nht)d oht=%(oht)g dht=%(dht)g' % par

    # frequency params
    par['fw']=24
    par['jw']=1
    par['dw']=1/(par['nt']*par['dt'])
    par['kw']=par['nt']/2+1
    par['ow']=par['fw']*par['dw']
    par['nw']=160

    par['nwmva']=160
    par['fwmva']=0
    par['jwmva']=1

    # migration domain params
    par['oximg']=par['ox']+20*par['dx']
    par['nximg']=1000
    par['dximg']=par['dx']
    
    par['ozimg']=par['oz']+par['nzdtm']*par['dz']
    par['jzimg']=2
    par['nzimg']=1+(par['nz']-par['nzdtm'])/par['jzimg']
    par['dzimg']= par['dz']*par['jzimg']

    par['xmin']=par['oximg']
    par['xmax']=par['oximg']+(par['nximg']-1)*par['dximg']
    par['zmin']=par['ozimg']
    par['zmax']=par['ozimg']+(par['nzimg']-1)*par['dzimg']

    # CIG tile params
    par['jcig']=100
    par['fcig']=50

    par['nxtile']=par['nximg']/par['jcig']*(par['nht'])
    par['dxtile']=(par['nximg']-1)*par['dximg']/par['nxtile']
    par['oxtile']=par['oximg']
    
    # picking parameters
    par['npck1']=50
    par['npck2']=50
    par['scpck']=1.0

# ------------------------------------------------------------
def setup(par):

    # ------------------------------------------------------------
    param(par)
    fdmod.param(par) 
    
    # ------------------------------------------------------------
    # shot positions
    # ------------------------------------------------------------
    shotcoord('ss',par)

    # ------------------------------------------------------------
    # slowness models
    # ------------------------------------------------------------
    importvels(par)
    
# ------------------------------------------------------------
def run(par):

    for vtag in (['C','L','H']):

        # ------------------------------------------------------------
        # migration
        # ------------------------------------------------------------
        migration('img'+vtag,
                  'cig'+vtag,
                  'drv'+vtag,
                  'dxs',
                  'dxr',
                  'wfs'+vtag,
                  'wfr'+vtag,
                  'slo'+vtag,
                  'ss',par)
        
    for vtag in (['L','H']):
        
        # ------------------------------------------------------------
        # image perturbation
        # ------------------------------------------------------------
        deltaimage('dim'+vtag,
                   'cig'+vtag,
                   'drv'+vtag,
                   'slo'+vtag,
                   'ss',par)
        
        # ------------------------------------------------------------
        # slowness backprojection
        # ------------------------------------------------------------
        adjoint('bsl'+vtag,
                'dim'+vtag,
                'wfs'+vtag,
                'wfr'+vtag,
                'slo'+vtag,
                'ss',par)
        
        # ------------------------------------------------------------
        # image perturbation
        # ------------------------------------------------------------
        forward('dsl'+vtag,
                'ipt'+vtag,
                'wfs'+vtag,
                'wfr'+vtag,
                'slo'+vtag,
                'ss',par)
        
        # ------------------------------------------------------------
        # slowness backprojection
        # ------------------------------------------------------------
        adjoint('aip'+vtag,
                'ipt'+vtag,
                'wfs'+vtag,
                'wfr'+vtag,
                'slo'+vtag,
                'ss',par)


# ------------------------------------------------------------
def importvels(par):
    # prepare velocity
    sigsbee.getmigvel('velo',par)
    sigsbee.getstrvel('vstr',par)
    Result('vstr',fdmod.cgrey('color=j allpos=y bias=1.5',par))
    Result('velo',fdmod.cgrey('color=j allpos=y bias=1.5',par))

    # padding in z
    Flow('vpad','vstr',
         '''
         window n1=1 f1=1200 | 
         spray axis=1 n=143 |
         smooth rect2=250 repeat=5
         ''' )

    Flow('vsed','vstr','window n1=1180' )

    Flow('velC','vsed vpad','cat axis=1 ${SOURCES[1]}')
    Flow('velCt','velC','transp')

    # masks
    sigsbee.makemask('velC','smask','wmask','lmask',par)

# ------------------------------------------------------------
def importref(par):
    # prepare reflectivity
    sigsbee.getreflect('ref_',par)
    Result('ref',fdmod.cgrey('pclip=99',par))

    # padding in z
    Flow('rpad','ref_',
         '''
         window n1=1 f1=1200 | 
         spray axis=1 n=143 |
         smooth rect2=250 repeat=5
         ''' )

    Flow('sub',None,
         '''
         spike nsp=1 mag=1
         n1=3201 d1=0.00762 o1=3.048 k1=0 l1=%(nx)d
         n2=1180 d2=0.00762 o2=0     k2=1179 l2=1179 |
         put label1=x label2=z unit1=km unit2=km | transp 
         ''' % par)

    Flow('rsed','ref_','window n1=1180' )
    Flow('rsed_','sub rsed','add ${SOURCES[1]}')
    Flow('ref','rsed_ rpad','cat axis=1 ${SOURCES[1]}')

# ------------------------------------------------------------
def globalmask(amask,par):
    Flow(amask,None,
         '''
         spike nsp=1 mag=1
         o1=%(oz)g d1=%(dz)g n1=%(nz)d
         o2=%(ox)g d2=%(dx)g n2=%(nx)d
         ''' % par, stdin=0)

# ------------------------------------------------------------
def localmask(amask,par):
    Flow(amask,None,
         '''
         spike nsp=1 mag=1
         o1=%(oz)g d1=0.00762 n1=1343 k1=500 l1=900
         o2=3.048  d2=0.02286 n2=1024 k2=000 l2=400
         ''' % par, stdin=0)

# ------------------------------------------------------------
def slownesses(amask,magn,s,ds,ss,par):

    # anomaly shape
    Flow('shapex',amask,'transp')
    
    Flow('shape1',[amask,'velC'],         'remap1 pattern=${SOURCES[1]}')
    Flow('shape2','shape1 velCt','transp | remap1 pattern=${SOURCES[1]} | transp')
    Flow('anom','shape2 lmask',
         '''
         smooth rect1=100 rect2=100 repeat=2 |
         scale axis=123 |
         math s=${SOURCES[1]} output="input*s"
         ''')
    Result('anom',fdmod.cgrey('allpos=y',par))
    
    # incorrect velocity
    Flow('velL',['velC','anom'],
         'math v=${SOURCES[0]} a=${SOURCES[1]} output="v*(1.0-(%g)*a)"' % magn)
    Flow('velH',['velC','anom'],
         'math v=${SOURCES[0]} a=${SOURCES[1]} output="v*(1.0+(%g)*a)"' % magn)
    # velocity plots
    for i in ('velC','velL','velH'):
        Result(i,fdmod.cgrey('color=j allpos=y bias=1.5',par))
        
    # ------------------------------------------------------------
    # prepare slowness
    spmig.slowness('svelC','velC',par)
    spmig.slowness('svelL','velL',par)
    spmig.slowness('svelH','velH',par)
    Flow('sdt','svelC','window squeeze=n n3=2 j3=%(nzdtm)d' %par)

    Flow(s+'C','svelC','window squeeze=n      f3=%(nzdtm)d j3=%(jzimg)d'%par) # correct velocity
    Flow(s+'L','svelL','window squeeze=n      f3=%(nzdtm)d j3=%(jzimg)d'%par) # low     velocity
    Flow(s+'H','svelH','window squeeze=n      f3=%(nzdtm)d j3=%(jzimg)d'%par) # high    velocity
    
    Flow('slows',[s+'C',s+'L',s+'H'],
         '''
         cat axis=2 space=n ${SOURCES[1:3]} |
         byte pclip=100 allpos=y gainpanel=a bias=0.221
         ''')

    Plot(s+'C','slows','window n2=1 f2=0 | transp |'+ fdmod.cgrey('color=j',par))
    Plot(s+'L','slows','window n2=1 f2=1 | transp |'+ fdmod.cgrey('color=j',par))
    Plot(s+'H','slows','window n2=1 f2=2 | transp |'+ fdmod.cgrey('color=j',par))

    for stag in (['L','H','C']):
        Result(s+stag,[s+stag,ss],'Overlay')

    for stag in (['L','H']):
        Flow(ds+stag,[s+stag,s+'C',amask,'shapex'],
             '''
             add scale=-1,1 ${SOURCES[1]} | window | transp |
             remap1 pattern=${SOURCES[2]} |
             transp |
             remap1 pattern=${SOURCES[3]} |
             window n2=%(nz)d min2=%(oz)g j2=%(jzimg)d |
             spray axis=2 n=1 o=0 d=1 |
             window squeeze=n min1=%(ox)g n1=%(nx)d d1=%(dx)g | rtoc
             ''' % par)
        Plot(  ds+stag,'real | window | transp |'+ fdmod.cgrey('color=e',par))
        Result(ds+stag,[ds+stag,ss],'Overlay')

# ------------------------------------------------------------
def shotcoord(ss,par):

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp
    
        xsou = par['os'] + (par['fs'] + iexp * par['js']) * par['ds']
        fdmod.point('ss'+etag,xsou,par['ozimg'],par)
        Plot(ss+etag,fdmod.ssplot('plotcol=5',par))
        
    Plot(ss,map(lambda x: ss+'-e%03d'%x,range(par['ns'])),
         'cat axis=3 space=n ${SOURCES[1:%d]} |' % par['ns']
         + fdmod.ssplot('plotcol=5',par))

# ------------------------------------------------------------
def originaldata(sdat,rdat,sd,par):

    Flow(  'wav',None,'spike nsp=1 mag=1 k1=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g' % par)
    Result('wav','window n1=500 |' + fdmod.waveplot('',par))

    sigsbee.getdata('data',par)
    sigsbee.makeshots('sall','data',par)
        
    Flow('swin',
         'sall',
         '''
         window n3=%(ns)d f3=%(fs)d j3=%(js)d |
         bandpass flo=2 fhi=10
         ''' % par)
    
    for file in (['swin','sall']):
        for iexp in range(par['ns']):
            etag = "-e%03d" % iexp

            Plot(file+etag,file,
                 'window n3=1 f3=%d | put o2=%g |' % (iexp,par['jos']+iexp*par['jds'])
                 + fdmod.dgrey('pclip=99',par))
            Result(file,map(lambda x: file+'-e%03d' % x,range(par['ns'])),'Movie')


    encode.shot2grid('sdata-t','rdata-t','wav','swin',par)
    Flow('S-sou-t','sdata-t','window | transp' %par)
    Flow('S-rec-t','rdata-t','window | transp' %par)

    encode.time2freq('S-sou-t','S-sou-w',par)
    encode.time2freq('S-rec-t','S-rec-w',par)
    
    # datum through water
    spmig.datum3('S-dfs',
                 'S-dfr',
                 sd,
                 'S-sou-w',
                 'S-rec-w',par) 

    # window datumed data (all shots)
    Flow(sdat,'S-dfs','window squeeze=n min1=%(oximg)g n1=%(nximg)g' % par)
    Flow(rdat,'S-dfr','window squeeze=n min1=%(oximg)g n1=%(nximg)g' % par)

    Result(sdat,fdmod.fgrey('pclip=99',par))
    Result(rdat,fdmod.fgrey('pclip=99',par))

# ------------------------------------------------------------
def simulateddata(sdat,rdat,slo,sd,ss,par):

    # single reflector
    par['kxref']=50
    par['lxref']=par['nximg']-par['kxref']
    par['kzref']=par['nzimg']-par['nzpad']/par['jzimg']-10
    Flow('ref',None,
         '''
         spike nsp=1 mag=1 
         n1=%(nximg)d d1=%(dximg)g o1=%(oximg)g k1=%(kxref)d l1=%(lxref)d
         n2=%(nzimg)d d2=%(dzimg)g o2=%(ozimg)g k2=%(kzref)d |
         smooth rect1=25 repeat=3 |
         spray axis=2 n=1 o=0 d=1 |
         put label1=x label2=y label3=z 
         ''' % par )    
    Plot('ref','window | transp | smooth rect1=3 |'
         + fdmod.cgrey('pclip=100',par))
    Result('ref',['ref',ss],'Overlay')

    # source wavelet
    par['frq']=5
    par['kt']=50
    Flow('wvl',None,
         '''
         spike nsp=1 mag=1
         n1=%(nt)d d1=%(dt)g o1=0   k1=%(kt)d |
         ricker1 frequency=%(frq)g |
         scale axis=123 |
         fft1 |
         window squeeze=n n1=%(nw)d min1=%(ow)g
         ''' % par)

    # source data on the surface
    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp

        xsou = par['os'] + (par['fs'] + iexp * par['js']) * par['ds']
        isou = (xsou-par['oximg'])/par['dximg']
        Flow('spk'+etag,'wvl',
             '''
             pad beg2=%d n2out=%d |
             put label1=w label2=x label3=y o2=%g d2=%g |
             transp memsize=250 plane=12 |
             transp memsize=250 plane=23 
             ''' % (isou,
                    par['nximg'],
                    par['oximg'],
                    par['dximg'],) )
    Flow('spk',
         map(lambda x: 'spk-e%03d'  % x,range(par['ns'])),
         'cat space=n axis=4 ${SOURCES[1:%d]}'%par['ns'])

    # datumed source data
    zomig.Cdtone3(sdat,'spk',sd,par)
    Result(sdat,fdmod.fgrey('',par))

    # datumed receiver data
    spmig.modelPW3(rdat,slo,sdat,'ref',par)
    Result(rdat,fdmod.fgrey('',par))


#    for iexp in range(par['ns']):
#        etag = "-e%03d" % iexp
#
#        Flow(sdat+etag,sdat,'window squeeze=n n4=1 f4=%d' % iexp)
#        spmig.modelPW3('tdt'+etag,slo,sdat+etag,'ref',par)
#
#    Flow(rdat,
#         map(lambda x: 'tdt-e%03d'  % x,range(par['ns'])),
#         'cat space=n axis=4 ${SOURCES[1:%d]}'%par['ns'])
#    Result(rdat,fdmod.fgrey('',par))
        
# ------------------------------------------------------------
def migration(img,cig,drv,sdat,rdat,swfl,rwfl,slo,ss,par):

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp

        # select datumed data for individual shots
        Flow(sdat+etag,sdat,'window n4=1 f4=%d squeeze=n' %iexp)
        Flow(rdat+etag,rdat,'window n4=1 f4=%d squeeze=n' %iexp)

        # window frequency slices for MVA
        Flow(sdat+etag+'-win',
             sdat+etag,'window squeeze=n n3=%(nwmva)d f3=%(fwmva)d j3=%(jwmva)d' % par)
        Flow(rdat+etag+'-win',
             rdat+etag,'window squeeze=n n3=%(nwmva)d f3=%(fwmva)d j3=%(jwmva)d' % par)

        # wavefields for individual shots
        zomig.Cwfone3(swfl+etag,sdat+etag+'-win',slo,par) #   source
        zomig.Awfone3(rwfl+etag,rdat+etag+'-win',slo,par) # receiver

        # migrate
        Flow([img+etag,cig+etag,drv+etag],
             [sdat+etag,rdat+etag,slo],
             '''
             ../Code/srmig3.x
             %s
             rwf=${SOURCES[1]}
             slo=${SOURCES[2]}
             cig=${TARGETS[1]}
             drv=${TARGETS[2]}
             ''' % spmig.param(par))
        
    # concatenate images and CIGs
    for k in ([img,drv,cig]):
        Flow(k+'-all',
             map(lambda x: k+'-e%03d'  % x,range(par['ns'])),
             'cat space=n axis=2 ${SOURCES[1:%d]}'%par['ns'])
        Flow(k,k+'-all','stack | spray axis=2 n=1 o=0 d=1')

    # plot complete images
    Plot(  img,'window | transp |'+ fdmod.cgrey('pclip=97',par))
    Result(img,[img,ss],'Overlay')
    
    for k in ([drv,cig]):
        # plot complete CIGs
        Result(k,'window n1=1 f1=%d |' % (par['nximg']/4)
               + adcig.tgrey('',par))

        # plot CIG tiles
        Result(k+'-tile',
               k,
               '''
               window j1=%(jcig)d f1=%(fcig)d |
               transp plane=12 |
               transp plane=23 |
               put n2=%(nxtile)d o2=%(oxtile)g d2=%(dxtile)g n3=1 |
               ''' % par
               + fdmod.cgrey('wantaxis2=n pclip=97',par))
    
    # clip image
    Flow(img+'-byt',
         img+'-all','byte gainpanel=a pclip=97')
    # clip CIGs
    Flow(cig+'-byt',
         cig+'-all','byte gainpanel=a pclip=97')

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp

        # plot partial images
        Plot(img+etag,img+'-byt',
             'window n2=1 f2=%d | transp |' % iexp + fdmod.cgrey('',par))
        Result(img+etag,[img+etag,ss+etag],'Overlay')

        # plot partial CIGs
        Result(cig+etag,
               cig+'-byt',               
               'window n1=1 f1=%d n2=1 f2=%d |' % (par['nximg']/4,iexp)
               + adcig.tgrey('',par))

# ------------------------------------------------------------
def tauplot(col,fat,custom,par):
    return '''
    graph pad=n transp=y yreverse=y wanttitle=n
    plotcol=%d 
    min1=%g max1=%g
    min2=%g max2=%g 
    plotfat=%d wantaxis=n screenratio=2.0 
    %s %s labelsz=4
    ''' % (col,
	   par['zmin'],par['zmax'],
           par['oht'],
           par['oht'] + par['nht']*par['dht'],
           fat,par['labelattr'],custom)


# ------------------------------------------------------------
def deltaimage(dim,cig,drv,slo,ss,par):

    # CIG envelope
    Flow(cig+'-env',
         cig,
         '''
         window |
         transp plane=12 |
         transp plane=23 |
         envelope |
         smooth rect1=50 rect2=50 repeat=1 |
         scale axis=123 |
         clip clip=0.75         
         ''')

    # pick tau
    Flow(cig+'-pck',
         cig+'-env',
         '''
         pick rect1=%(npck1)d rect2=%(npck2)d rect3=1 vel0=0.0 |
         scale rscale=%(scpck)g |
         transp plane=13
         ''' % par)
    Result(cig+'-pck',
           'window | transp |' + fdmod.cgrey('color=j',par))

    # plot picked tau on CIG 
    Flow(cig+'-sbt',
         cig,'byte gainpanel=a pclip=90')
    for ipos in range(0,par['nximg'],10):
        ptag = "-x%03d" % ipos

        Plot(cig+'-pck'+ptag,
             cig+'-pck',
             'window n1=1 f1=%d |' % ipos + tauplot(6,4,'screenratio=3',par) )
        Plot(cig+ptag,
             cig+'-sbt',
             'window n1=1 f1=%d |' % ipos + adcig.tgrey('screenratio=3',par))
        Result(cig+ptag,[cig+ptag,cig+'-pck'+ptag],'Overlay')

    # velocity for tau spreading
    Flow(cig+'-vel',
         slo,
         '''
         window squeeze=n min1=%(oximg)g n1=%(nximg)d j1=%(jximg)d |
         math output="1/input"
         ''' % par)
    # tau spreading
    Flow(cig+'-tau',[drv,cig+'-pck',cig+'-vel'],
         '''
         ../Code/tauspread.x
         fit=${SOURCES[1]}
         vel=${SOURCES[2]} zmin=2.5
         ''')

    # delta image
    Flow(dim,[drv,cig+'-tau'],'add mode=p ${SOURCES[1]}')

    # plot complete image perturbation
    Plot(dim,
         'window n4=1 min4=0 | transp |' + fdmod.cgrey('pclip=97',par))
    Result(dim,[dim,ss],'Overlay')
    
    # plot DIM tiles
    Result(dim+'-tile',
           dim,
           '''
           window j1=%(jcig)d f1=%(fcig)d |
           transp plane=12 |
           transp plane=23 |
           put n2=%(nxtile)d o2=%(oxtile)g d2=%(dxtile)g n3=1 |
           ''' % par
           + fdmod.cgrey('wantaxis2=n pclip=97',par))

# ------------------------------------------------------------
def adjoint(bsl,dim,swfl,rwfl,slo,ss,par):

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp

        Flow(bsl+etag,
             [dim,swfl+etag,rwfl+etag,slo],
             '''
             rtoc |
             ../Code/srmvatau.x
             adj=y nht=%g oht=%g dht=%g %s
             swf=${SOURCES[1]}
             rwf=${SOURCES[2]}
             slo=${SOURCES[3]} |
             real
             ''' % (par['nht'],par['oht'],par['dht'],spmig.param(par)))

    Flow(bsl+'-all',
         map(lambda x: bsl+'-e%03d'  % x,range(par['ns'])),
         'cat space=n axis=2 ${SOURCES[1:%d]}'%par['ns'])
    Flow(bsl,bsl+'-all','stack | spray axis=2 n=1 o=0 d=1')

    # plot complete slowness backprojection
    Plot(bsl,
         'window | transp |' + fdmod.cgrey('color=e pclip=97',par))
    Result(bsl,[bsl,ss],'Overlay')

    # clip slowness perturbation
    Flow(bsl+'-byt',
         bsl+'-all','byte gainpanel=a pclip=97')

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp
        
        # plot partial slowness backprojection
        Plot(bsl+etag,bsl+'-byt',
             'window n2=1 f2=%d | transp |' % iexp + fdmod.cgrey('color=e',par))
        Result(bsl+etag,[bsl+etag,ss+etag],'Overlay')


# ------------------------------------------------------------
def forward(dsl,dim,swfl,rwfl,slo,ss,par):

    for iexp in range(par['ns']):
        etag = "-e%03d" % iexp
        
        Flow(dim+etag,
             [dsl,swfl+etag,rwfl+etag,slo],
             '''
             rtoc |
             ../Code/srmvatau.x
             adj=n nht=%g oht=%g dht=%g %s
             swf=${SOURCES[1]}
             rwf=${SOURCES[2]}
             slo=${SOURCES[3]} |
             real
             ''' % (par['nht'],par['oht'],par['dht'],spmig.param(par)))

    Flow(dim+'-all',
         map(lambda x: dim+'-e%03d'  % x,range(par['ns'])),
         'cat space=n axis=2 ${SOURCES[1:%d]}'%par['ns'])
    Flow(dim,dim+'-all','stack | spray axis=2 n=1 o=0 d=1')

    # plot complete image perturbation
    Plot(dim,
         'window n4=1 min4=0 | transp |' + fdmod.cgrey('pclip=97',par))
    Result(dim,[dim,ss],'Overlay')
    
    # plot DIM tiles
    Result(dim+'-tile',
           dim,
           '''
           window j1=%(jcig)d f1=%(fcig)d |
           transp plane=12 |
           transp plane=23 |
           put n2=%(nxtile)d o2=%(oxtile)g d2=%(dxtile)g n3=1 |
           ''' % par
           + fdmod.cgrey('wantaxis2=n pclip=97',par))
    

    # clip image perturbation
#    Flow(dim+'-byt',
#         dim+'-all','byte gainpanel=a pclip=99')

#    for iexp in range(par['ns']):
#        etag = "-e%03d" % iexp
#        
#        # plot partial image perturbation
#        Plot(dim+etag,dim+'-byt',
#             'window n4=1 min4=0 n2=1 f2=%d | transp |' % iexp + fdmod.cgrey('',par))
#        Result(dim+etag,[dim+etag,ss+etag],'Overlay')
