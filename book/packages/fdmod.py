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

    if(not par.has_key('tmin')):     par['tmin']=par['ot']
    if(not par.has_key('tmax')):     par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    if(not par.has_key('xmin')):     par['xmin']=par['ox']
    if(not par.has_key('xmax')):     par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    if(not par.has_key('zmin')):     par['zmin']=par['oz']
    if(not par.has_key('zmax')):     par['zmax']=par['oz'] + (par['nz']-1) * par['dz']

    if(not par.has_key('ratio')):    par['ratio']=(par['zmax']- par['zmin'])/(par['xmax']- par['xmin'])
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
def rtm(imag,sdat,rdat,velo,dens,sacq,racq,custom,mem,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat)
    rout = imag+'_dr' # receiver data (not the input rdat)

    # source wavefield (z,x,t)
    awe(sout,swfl,sdat,velo,dens,sacq,sacq,custom+' jsnap=1 ',par)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    twfl = imag+'_tur'
    tout = imag+'_tdr'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awe(tout,twfl,tdat,velo,dens,racq,racq,custom+' jsnap=1 ',par)
    Flow(rwfl,twfl,'reverse which=4 opt=i verb=y memsize=1000')
    Flow(rout,tout,'reverse which=2 opt=i verb=y')

    corr = imag+'_cor'

    # conventional (cross-correlation zero-lag) imaging condition
    Flow(corr,[swfl,rwfl],'paradd mode=p ${SOURCES[1]} memsize=%d' %mem)
    Flow(imag,corr,'stack axis=3')

# exploding-reflector reverse-time migration
def zom(imag,data,rdat,velo,dens,sacq,racq,custom,par):

    rwfl = imag+'_ur' # receiver wavefield
    rout = imag+'_dr' # receiver data (not the input rdat)

    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    twfl = imag+'_tur'

    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awe(data,twfl,tdat,velo,dens,sacq,racq,custom,par)

    Flow(imag,twfl,'window n3=1 f3=%d' % (par['nt']/par['jsnap']) )

# wavefield-over-model plots
def wom(wom,wfld,velo,par):

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
         scale axis=123 |
         spray axis=3 n=%d o=%g d=%g |
         math w=${SOURCES[1]} output="(input-0.5)+10*w"
         ''' % (par['nt']/par['jsnap'],
                par['ot'],
                par['dt']*par['jsnap']))
