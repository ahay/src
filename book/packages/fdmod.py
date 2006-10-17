from rsfproj import *

# overwrite default parameters
def param(par):
    if(not par.has_key('nbz')):      par['nbz']=100
    if(not par.has_key('nbx')):      par['nbx']=100
    if(not par.has_key('tz')):       par['tz']=0.0035
    if(not par.has_key('tx')):       par['tx']=0.0035

    if(not par.has_key('snap')):     par['snap']='y'
    if(not par.has_key('jsnap')):    par['jsnap']=100

    if(not par.has_key('ompchunk')): par['ompchunk']=1

    if(not par.has_key('ot')): par['ot']=0.
    if(not par.has_key('nt')): par['nt']=1
    if(not par.has_key('dt')): par['dt']=1.

    if(not par.has_key('tmin')):     par['tmin']=par['ot']
    if(not par.has_key('tmax')):     par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    if(not par.has_key('xmin')):     par['xmin']=par['ox']
    if(not par.has_key('xmax')):     par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    if(not par.has_key('zmin')):     par['zmin']=par['oz']
    if(not par.has_key('zmax')):     par['zmax']=par['oz'] + (par['nz']-1) * par['dz']

    if(not par.has_key('ratio')):    par['ratio']=(par['zmax']-par['zmin'])/(par['xmax']-par['xmin'])
    if(not par.has_key('height')):   par['height']=par['ratio']*14

# ------------------------------------------------------------
# plotting functions
def cgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title="" wantscalebar=n
    pclip=100
    min1=%g max1=%g %s
    min2=%g max2=%g %s
    screenratio=%g screenht=%g
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],
           par['xmin'],par['xmax'],par['lx'],
           par['ratio'],par['height'],
           custom)

def wgrey(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g |
    grey labelrot=n wantaxis=y title="" wantscalebar=n
    pclip=100 gainpanel=a
    screenratio=%g screenht=%g
    %s %s %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['ratio'],par['height'],
           par['lz'],par['lx'],
           custom)

def cgraph(custom,par):
    return '''
    graph labelrot=n wantaxis=n title="" yreverse=y 
    min2=%g max2=%g %s
    min1=%g max1=%g %s
    screenratio=%g screenht=%g
    %s
    ''' % (
           par['zmin'],par['zmax'],par['lz'],
           par['xmin'],par['xmax'],par['lx'],
           par['ratio'],par['height'],
           custom)

def dgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100
    min1=%g max1=%g %s
    min2=%g max2=%g %s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],
           par['xmin'],par['xmax'],par['lx'],
           custom)

# ------------------------------------------------------------
# create wavelet
def wavelet(wav,frequency,par):
    par['frequency'] = frequency
    
    Flow(wav,None,
         '''
         spike nsp=1 mag=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g k1=%(kt)d |
         ricker1 frequency=%(frequency)g |
         scale axis=123 |
         put label1=t 
         ''' % par)    

# ------------------------------------------------------------
def horizontal(cc,coord,par):
    par['coord'] = coord
    
    Flow(cc+'_',None,'math n1=%(nx)d d1=%(dx)g o1=%(ox)g output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%(coord)g" ' % par)
    Flow(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def vertical(cc,coord,par):
    par['coord'] = coord
    
    Flow(cc+'_',None,'math n1=%(nz)d d1=%(dz)g o1=%(oz)g output=0' % par)
    Flow(cc+'_x',cc+'_','math output="%(coord)g" ' % par)
    Flow(cc+'_z',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp
         ''', stdin=0)

def point(cc,xcoord,zcoord,magnitude,par):

    Flow(cc+'_',None,'math n1=1 d1=1 o1=0 output=0' % par)
    Flow(cc+'_z',cc+'_','math output="%g"' % zcoord)
    Flow(cc+'_x',cc+'_','math output="%g"' % xcoord)
    Flow(cc+'_r',cc+'_','math output="%g"' % magnitude)    
    Flow(cc,[cc+'_x',cc+'_z',cc+'_r'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | transp
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

# F-D modeling from arbitrary source/receiver geometry
def awe(odat,wfld,  idat,velo,dens,sou,rec,custom,par):
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

# shot-record reverse-time migration
def rtm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat)
    rout = imag+'_dr' # receiver data (not the input rdat)

    # source wavefield (z,x,t)
    awe(sout,swfl,sdat,velo,dens,sacq,iacq,custom,par)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awe(tout,rwfl,tdat,velo,dens,racq,iacq,custom,par)
    Flow(rout,tout,'reverse which=2 opt=i verb=y')

    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[sout,rout],'xcor uu=${SOURCES[1]} axis=2 verb=y nbuf=500')
    
#    corr = imag+'_cor'
#    Flow(corr,[sout,rout],'paradd mode=p ${SOURCES[1]} memsize=%d' %mem)
#    Flow(imag,corr,'stack axis=2')

# exploding-reflector reverse-time migration
def zom(imag,data,rdat,velo,dens,sacq,racq,custom,par):

    rwfl = imag+'_ur' # receiver wavefield
    rout = imag+'_dr' # receiver data (not the input rdat)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    twfl = imag+'_tur'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awe(data,twfl,tdat,velo,dens,sacq,racq,custom,par)

    Flow(imag,twfl,'window n3=1 f3=%d' % (par['nt']/par['jsnap']-1) )

# wavefield-over-model plots
def wom(wom,wfld,velo,vmean,par):

    chop = wfld+'_chop'
    Flow(chop,wfld,
         '''
         window
         min1=%(zmin)g max1=%(zmax)g
         min2=%(xmin)g max2=%(xmax)g |
         scale axis=123
         ''' % par)

    Flow(wom,[velo,chop],
         '''
         add add=-%g |
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g |
         math w=${SOURCES[1]} output="0.005*input+w"
         ''' % (vmean,
                par['nt']/par['jsnap'],
                par['ot'],
                par['dt']*par['jsnap']))
