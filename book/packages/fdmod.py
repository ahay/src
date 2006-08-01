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
    grey labelrot=n wantaxis=y wanttitle=n wantscalebar=n
    title="" pclip=99
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
    grey labelrot=n wantaxis=y wanttitle=y wantscalebar=n
    title="" pclip=100 gainpanel=a
    screenratio=%g screenht=%g
    %s %s %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['ratio'],par['height'],
           par['lz'],par['lx'],
           custom)

def cgraph(custom,par):
    return '''
    graph labelrot=n
    yreverse=y wantaxis=n title=" " 
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
    grey labelrot=n wantaxis=y wanttitle=y
    title="" pclip=99
    min2=%g max2=%g %s
    %s
    ''' % (par['xmin'],par['xmax'],par['lx'],
           custom)

# execute acoustic finite-differences modeling
def amodel(data,wfld,  wavl,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [data,wfld],[wavl,velo,dens,sou,rec],
          '''
          afmodP ompchunk=%(ompchunk)d
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

