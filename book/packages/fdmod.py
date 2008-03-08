from rsfproj import *

# default parameters
def param(par):
    if(not par.has_key('lt')):       par['lt']='t'
    if(not par.has_key('lz')):       par['lz']='z'
    if(not par.has_key('lx')):       par['lx']='x'

    if(not par.has_key('ut')):       par['ut']='s'
    if(not par.has_key('uz')):       par['uz']='m'
    if(not par.has_key('ux')):       par['ux']='m'

    if(not par.has_key('nb')):       par['nb']=0
    if(not par.has_key('nbz')):      par['nbz']=100
    if(not par.has_key('nbx')):      par['nbx']=100
    if(not par.has_key('tz')):       par['tz']=0.0035
    if(not par.has_key('tx')):       par['tx']=0.0035

    if(not par.has_key('nbell')):    par['nbell']=1

    if(not par.has_key('snap')):     par['snap']='y'
    if(not par.has_key('jsnap')):    par['jsnap']=100
    if(not par.has_key('jdata')):    par['jdata']=1

    if(not par.has_key('ompchunk')): par['ompchunk']=1
    if(not par.has_key('ompnth')):   par['ompnth']=0
    if(not par.has_key('free')):     par['free']='n'
   
    if(not par.has_key('ot')):       par['ot']=0.
    if(not par.has_key('nt')):       par['nt']=1
    if(not par.has_key('dt')):       par['dt']=1.

    if(not par.has_key('tmin')):     par['tmin']=par['ot']
    if(not par.has_key('tmax')):     par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    if(not par.has_key('xmin')):     par['xmin']=par['ox']
    if(not par.has_key('xmax')):     par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    if(not par.has_key('zmin')):     par['zmin']=par['oz']
    if(not par.has_key('zmax')):     par['zmax']=par['oz'] + (par['nz']-1) * par['dz']

    if(not par.has_key('ratio')):    par['ratio']=1.0*(par['zmax']-par['zmin'])/(par['xmax']-par['xmin'])
    if(not par.has_key('height')):   par['height']=par['ratio']*14
    if(par['height']>10): par['height']=10

    if(not par.has_key('scalebar')): par['scalebar']='n'
    if(not par.has_key('labelattr')): par['labelattr']=' labelsz=5 labelfat=3 titlefat=3 '

    if(not par.has_key('nq1')): par['nq1']=par['nz']
    if(not par.has_key('oq1')): par['oq1']=par['oz']
    if(not par.has_key('dq1')): par['dq1']=par['dz']

    if(not par.has_key('nq2')): par['nq2']=par['nx']
    if(not par.has_key('oq2')): par['oq2']=par['ox']
    if(not par.has_key('dq2')): par['dq2']=par['dx']
    
    
    par['labelattr']=' '+par['labelattr']+' '
# ------------------------------------------------------------
# plotting functions
def cgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def wgrey(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g |
    grey labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    label1=%s unit1=%s
    label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def cgraph(custom,par):
    return '''
    graph labelrot=n wantaxis=n title="" yreverse=y
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def ccont(custom,par):
    return '''
    contour labelrot=n wantaxis=n title=""
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    plotcol=2 plotfat=3
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+' '+custom)

def dgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['labelattr']+' '+custom)

def egrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
# plot wavelet
def waveplot(custom,par):
    return '''
    graph min2=-1 max2=+1 title=""
    plotfat=5 plotcol=2
    label1=%s unit1=%s label2="" unit2=""
    %s
    ''' % (par['lt'],par['ut'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
# create wavelet
def wavelet(wav,frequency,par):
    par['frequency'] = frequency
    
    Flow(wav,None,
         '''
         spike nsp=1 mag=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g k1=%(kt)d |
         pad end1=%(nt)d |
         ricker1 frequency=%(frequency)g |
         window n1=%(nt)d |
         scale axis=123 |
         put label1=t 
         ''' % par)    

# ------------------------------------------------------------
def horizontal(cc,coord,par):
    Flow(cc+'_',None,'math n1=%(nx)d d1=%(dx)g o1=%(ox)g output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g" ' % coord)
    Flow(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def vertical(cc,coord,par):
    Flow(cc+'_',None,'math n1=%(nz)d d1=%(dz)g o1=%(oz)g output=0' % par)
    Flow(cc+'_x',cc+'_','math output="%g" '% coord)
    Flow(cc+'_z',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def point(cc,xcoord,zcoord,par):
    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def point3(cc,xcoord,zcoord,magn,par):
    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc+'_r',cc+'_','math output="%g"' % magn)
    Flow(cc,[cc+'_x',cc+'_z',cc+'_r'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | transp
         ''', stdin=0)

def circle(cc,xcenter,zcenter,radius,sampling,par):
    Flow(cc+'_x',None,
         'math n1=%d d1=%g o1=%g output="%g+%g*cos(x1/180*3.14)"'
         % (sampling,360./sampling,0.,xcenter,radius) )
    Flow(cc+'_z',None,
         'math n1=%d d1=%g o1=%g output="%g+%g*sin(x1/180*3.14)"'
         % (sampling,360./sampling,0.,zcenter,radius) )
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)


def boxarray(cc,nz,oz,dz,nx,ox,dx,par):
    Flow(cc+'_',None,
         '''
         math output=1
         n1=%d d1=%g o1=%g
         n2=%d d2=%g o2=%g
         ''' % (nz,dz,oz,nx,dx,ox) )
    Flow(cc+'_z',cc+'_','math output="x1" | put n1=%d n2=1' % (nz*nx))
    Flow(cc+'_x',cc+'_','math output="x2" | put n1=%d n2=1' % (nz*nx))
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def ssplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=o plotcol=6 plotfat=10 wantaxis=n %s' % custom,par)

def rrplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=4 plotfat=5 wantaxis=n %s' % custom,par)

def qqplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=1 plotfat=3 wantaxis=n %s' % custom,par)

def qqwin(par):
    return '''
    nq1=%(nq1)d oq1=%(oq1)g dq1=%(dq1)g
    nq2=%(nq2)d oq2=%(oq2)g dq2=%(dq2)g
    ''' % par

# ------------------------------------------------------------
# rays plot
def rayplot(hwt,j1ray,j2ray,j1wft,j2wft,custom,par):

    Plot(hwt+'ray',hwt,'window j1=%d j2=%d f2=%d | transp |' %(j1ray,j2ray,j2wft)
         + cgraph('plotcol=1 wantaxis=n '+custom,par))
    Plot(hwt+'wft',hwt,'window j1=%d j2=%d f2=%d |'          %(j1wft,j2wft,j2wft)
         + cgraph('plotcol=2 wantaxis=n symbol=. '+custom,par))

    Plot  (hwt,[hwt+'ray',hwt+'wft'],'Overlay')
  
# ------------------------------------------------------------
# execute acoustic finite-differences modeling
def amodel(data,wfld,  wavl,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [data,wfld],[wavl,velo,dens,sou,rec],
          '''
          afmod ompchunk=%(ompchunk)d
          verb=y abc=y free=n dens=y
          snap=%(snap)s jsnap=%(jsnap)d
          nbz=%(nbz)d tz=%(tz)g
          nbx=%(nbx)d tx=%(tx)g
          vel=${SOURCES[1]}
          den=${SOURCES[2]}
          sou=${SOURCES[3]}
          rec=${SOURCES[4]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

def lmodel(data,wfld,ldata,lwfld,  wavl,velo,refl,sou,rec,custom,par):
    par['fdcustom'] = custom

    Flow([ data,wfld,ldata,lwfld],[wavl,velo,refl,sou,rec],
         '''
         aborn ompchunk=%(ompchunk)d
         verb=y abc=y free=n
         snap=%(snap)s jsnap=%(jsnap)d
         nbz=%(nbz)d tz=%(tz)g
         nbx=%(nbx)d tx=%(tx)g
         vel=${SOURCES[1]}
         ref=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         lid=${TARGETS[2]}
         liw=${TARGETS[3]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# isotropic stiffness tensor
def isotropic(cc,vp,vs,ro,par):
    
    # Lame parameters
    # lambda = ro * (vp^2 - 2 vs^2)
    Flow(cc+'lambda',[ro,vp,vs],
         '''
         math output="ro*(vp*vp-2*vs*vs)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         vs=${SOURCES[2]}     
         ''')
    
    #     mu = ro *           vs^2
    Flow(cc+'mu',[ro,vs],
         '''
         math output="ro*vs*vs"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')
    
    # c11 = lambda + 2 mu
    # c33 = lambda + 2 mu
    # c13 = lambda
    # c44 = mu
    Flow(cc+'11',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'13',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l"')
    Flow(cc+'33',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'44',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="m"')
    Flow(cc,[cc+'11',cc+'13',cc+'33',cc+'44'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# anisotropic stiffness tensor
def anisotropic(cc,vp,vs,ro,epsilon,delta,par):
    Flow(cc+'33',[vp,ro],
         '''
         math output="ro*vp*vp"
         vp=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')    
    Flow(cc+'44',[vs,ro],
         '''
         math output="ro*vs*vs"
         vs=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')
    Flow(cc+'11',[cc+'33',epsilon],
         '''
         math output="2*epsilon*c33+c33"
         c33=${SOURCES[0]}
         epsilon=${SOURCES[1]}
         ''')
    Flow(cc+'13',[cc+'33',cc+'44',delta],
         '''
         math output="sqrt(2*c33*(c33-c44)*delta+(c33-c44)*(c33-c44))-c44"
         c33=${SOURCES[0]}
         c44=${SOURCES[1]}
         delta=${SOURCES[2]}
         ''')
    
    Flow(cc,[cc+'11',cc+'13',cc+'33',cc+'44'],
         'cat axis=3 space=n ${SOURCES[1:4]}')
    
# ------------------------------------------------------------
# acoustic modeling
def awefd(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)
def awefd1(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    awefd(odat,owfl,idat,velo,dens,sou,rec,custom+' expl=y ',par)
    

def lwefd(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([bdat,bwfl,sdat,swfl],[idat,velo,dens,refl,sou,rec],
         '''
         lwefd
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         ref=${SOURCES[3]}
         sou=${SOURCES[4]}
         rec=${SOURCES[5]}
         wfl=${TARGETS[1]}
         lid=${TARGETS[2]}
         liw=${TARGETS[3]}
         %(fdcustom)s
         ''' % par)
def lwefd1(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom,par):
    lwefd(bdat,bwfl,sdat,swfl,idat,velo,dens,refl,sou,rec,custom+' expl=y ',par)

# ------------------------------------------------------------
# elastic modeling
def ewefd(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d nbell=%(nbell)d
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# heat modeling
def hdefd(dat,wfl,  wav,con,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [dat,wfl],[wav,con,sou,rec],
          '''
          hdefd
          verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
          con=${SOURCES[1]}
          sou=${SOURCES[2]}
          rec=${SOURCES[3]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

# ------------------------------------------------------------    
# F-D modeling from arbitrary source/receiver geometry
def awe(odat,wfld,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,wfld],[idat,velo,dens,sou,rec],
          '''
          awe ompchunk=%(ompchunk)d
          verb=y abc=y free=n dens=y
          snap=%(snap)s jsnap=%(jsnap)d
          nbz=%(nbz)d tz=%(tz)g
          nbx=%(nbx)d tx=%(tx)g
          vel=${SOURCES[1]}
          den=${SOURCES[2]}
          sou=${SOURCES[3]}
          rec=${SOURCES[4]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

# ------------------------------------------------------------
# shot-record reverse-time migration
def rtm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    awefd(sout,swfl,sdat,velo,dens,sacq,iacq,custom,par)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(tout,rwfl,tdat,velo,dens,racq,iacq,custom,par)
    Flow(rout,tout,'reverse which=2 opt=i verb=y')

    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[sout,rout],'xcor uu=${SOURCES[1]} axis=2 verb=y nbuf=500')
    
#    corr = imag+'_cor'
#    Flow(corr,[sout,rout],'paradd mode=p ${SOURCES[1]} memsize=%d' %mem)
#    Flow(imag,corr,'stack axis=2')

# ------------------------------------------------------------
# exploding-reflector reverse-time migration
def zom(imag,data,rdat,velo,dens,sacq,racq,custom,par):

    rwfl = imag+'_ur' # receiver wavefield
    rout = imag+'_dr' # receiver data (not the input rdat)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    twfl = imag+'_tur'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(data,twfl,tdat,velo,dens,sacq,racq,custom+' jsnap=%d' % par['nt'],par)

    Flow(imag,twfl,'window n3=1')

# ------------------------------------------------------------
# wavefield-over-model plot
def wom(wom,wfld,velo,vmean,par):

    if(not par.has_key('wweight')): par['wweight']=10
    if(not par.has_key('wclip')):   par['wclip']=1.0

    chop = wfld+'_chop'
    Flow(chop,wfld,
         '''
         window
         min1=%(zmin)g max1=%(zmax)g
         min2=%(xmin)g max2=%(xmax)g |
         scale axis=123 |
         clip clip=%(wclip)g
         ''' % par)

    Flow(wom,[velo,chop],
         '''
         add add=-%g |
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g |
         math w=${SOURCES[1]} output="input+%g*w"
         ''' % (vmean,
                par['nt']/par['jsnap'],
                par['ot'],
                par['dt']*par['jsnap'],
                par['wweight']))

# ------------------------------------------------------------
# (elastic) wavefield-over-model
def wem(wom,wfld,velo,vmean,par):

    if(not par.has_key('wweight')): par['wweight']=10
    if(not par.has_key('wclip')):   par['wclip']=1.0

    Flow(velo+'-spray',
         velo,
         '''
         add add=-%g |
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g
         ''' % (vmean,
                par['nt']/par['jsnap'],
                par['ot'],
                par['dt']*par['jsnap']))

    for i in range(2):
        Flow(wom+'-%d' %i ,wfld,'window n3=1 f3=%d' %i)
        Flow(wom+'-%d' %i +'chop',
             wom+'-%d' %i,
             '''
             window
             min1=%(zmin)g max1=%(zmax)g
             min2=%(xmin)g max2=%(xmax)g |
             scale axis=123 |
             clip clip=%(wclip)g
             ''' % par)

        Flow(wom+'-%d' %i +'-wom',
             [velo+'-spray',wom+'-%d' %i + 'chop'],
             'math w=${SOURCES[1]} output="input+%g*w" | transp plane=34' % par['wweight'])

    Flow(wom,[wom+'-0-wom',wom+'-1-wom'],'cat axis=3 space=n ${SOURCES[1]}')
    
        
# ------------------------------------------------------------
# image-over-model plot
def iom(iom,imag,velo,vmean,par):

    if(not par.has_key('iweight')): par['iweight']=10
    if(not par.has_key('iclip')):   par['iclip']=1.0

    chop = imag+'_chop'
    Flow(chop,imag,
         '''
         window
         min1=%(zmin)g max1=%(zmax)g
         min2=%(xmin)g max2=%(xmax)g |
         scale axis=123 |
         clip clip=%(iclip)g
         ''' % par)

    Flow(iom,[velo,chop],
         '''
         add add=-%g |
         scale axis=123 |
         math w=${SOURCES[1]} output="input+%g*w"
         ''' % (vmean,par['iweight']))

# ------------------------------------------------------------
# wavefield snapshot plots
def wframe(frame,movie,index,custom,par):

    Flow([movie+'_plt',movie+'_bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)

#    Result(frame,[movie+'_plt',movie+'_bar'],
#           'window n3=1 f3=%d bar=${SOURCES[1]} |'% index + wgrey(custom,par))
    Plot  (frame,[movie+'_plt',movie+'_bar'],
           'window n3=1 f3=%d bar=${SOURCES[1]} |'% index + wgrey(custom,par))
    
# ------------------------------------------------------------
# elastic wavefield movie
def emovie(wfld,custom,axis,par):

    # loop over wavefield components
    for i in range(2):
        Flow(wfld+str(i+1),wfld,
             '''
             window n3=1 f3=%d |
             window min1=%g max1=%g min2=%g max2=%g
             ''' % (i,par['zmin'],par['zmax'],par['xmin'],par['xmax']))

    # join component wavefields
    Flow(wfld+'all',[wfld+'1',wfld+'2'],
         'cat axis=%d space=n ${SOURCES[1:2]}' % axis)

    if(axis==1):        
        height=2*par['height']
        if(height>10): height=10
        ratio =2*par['ratio']
    if(axis==2):
        height=par['height']
        if(height>10): height=10
        ratio =0.5*par['ratio']
    
    Result(wfld,wfld+'all',
           '''
           grey title="" wantaxis=y screenratio=%f screenht=%f
           gainpanel=a pclip=99 %s
           %s
           ''' % (ratio,height,par['labelattr'],custom) )

# ------------------------------------------------------------
# elastic wavefield movie
def eframe(frame,movie,index,custom,axis,par):

    if(axis==1):        
        height=2*par['height']
        if(height>10): height=10
        ratio =2*par['ratio']
    if(axis==2):
        height=par['height']
        if(height>10): height=10
        ratio =0.5*par['ratio']

    Flow([movie+'_plt',movie+'_bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)
    
    Result(frame,[movie+'_plt',movie+'_bar'],
           'window n3=1 f3=%d bar=${SOURCES[1]} |'% index +
           '''
           grey title="" wantaxis=y screenratio=%f screenht=%f
           gainpanel=a pclip=99 %s
           %s
           ''' % (ratio,height,par['labelattr'],custom) )
    
# ------------------------------------------------------------
# plot elastic wavelet
def ewavelet(wavelet,custom,par):
    
    for i in range(2):
        Plot(wavelet+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             waveplot('%d %s'% (i,custom) ,par))
    Result(wavelet,[wavelet+'1',wavelet+'2'],'SideBySideIso')
