try:    from rsf.cluster import *
except: from rsf.proj    import *

# default parameters
def param(par):
    
    if(not par.has_key('ot')):       par['ot']=0.
    if(not par.has_key('nt')):       par['nt']=1
    if(not par.has_key('dt')):       par['dt']=1.
    if(not par.has_key('lt')):       par['lt']='t'
    if(not par.has_key('ut')):       par['ut']='s'
            
    if(not par.has_key('ox')):       par['ox']=0.
    if(not par.has_key('nx')):       par['nx']=1
    if(not par.has_key('dx')):       par['dx']=1.
    if(not par.has_key('lx')):       par['lx']='x'
    if(not par.has_key('ux')):       par['ux']='km'
    
    if(not par.has_key('oy')):       par['oy']=0.
    if(not par.has_key('ny')):       par['ny']=1
    if(not par.has_key('dy')):       par['dy']=1.
    if(not par.has_key('ly')):       par['ly']='y'    
    if(not par.has_key('uy')):       par['uy']='km'
     
    if(not par.has_key('oz')):       par['oz']=0.
    if(not par.has_key('nz')):       par['nz']=1
    if(not par.has_key('dz')):       par['dz']=1.
    if(not par.has_key('lz')):       par['lz']='z'
    if(not par.has_key('uz')):       par['uz']='km'

    if(not par.has_key('tmin')):     par['tmin']=par['ot']
    if(not par.has_key('tmax')):     par['tmax']=par['ot'] + (par['nt']-1) * par['dt']
    if(not par.has_key('xmin')):     par['xmin']=par['ox']
    if(not par.has_key('xmax')):     par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    if(not par.has_key('ymin')):     par['ymin']=par['oy']
    if(not par.has_key('ymax')):     par['ymax']=par['oy'] + (par['ny']-1) * par['dy']
    if(not par.has_key('zmin')):     par['zmin']=par['oz']
    if(not par.has_key('zmax')):     par['zmax']=par['oz'] + (par['nz']-1) * par['dz']

    dx=par['xmax']-par['xmin'];
    dy=par['ymax']-par['ymin'];
    dz=par['zmax']-par['zmin'];
    dt=par['tmax']-par['tmin'];
   
    if(not par.has_key('iratio')):
        if(dx==0.0): par['iratio']=1.0
        else:        par['iratio']=1.0*(dz)/(dx)

    if(not par.has_key('iheight')):
        if(par['iratio']>=1): par['iheight']=10
        else:                 par['iheight']=14*par['iratio']

    if(not par.has_key('dratio')):
        par['dratio']=par['iratio']

    if(not par.has_key('dheight')):
       par['dheight']=par['iheight']

    if(not par.has_key('scalebar')): par['scalebar']='n'    
    if(not par.has_key('labelattr')): par['labelattr']=' parallel2=n labelsz=6 labelfat=3 titlesz=12 titlefat=3 xll=2 ' + ' '
    
    par['labelrot0']=' parallel2=n format1=%3.0f format2=%3.0f format3=%3.0f '
    par['labelrot1']=' parallel2=n format1=%3.1f format2=%3.1f format3=%3.1f '
    par['labelrot2']=' parallel2=n format1=%3.2f format2=%3.2f format3=%3.2f '
    
# ------------------------------------------------------------
# grey 2D image
def igrey2d(custom,par):
    return '''
    grey title="" pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['iratio'],par['iheight'],par['scalebar'],
           par['labelattr']+custom)

# grey 2D data
def dgrey2d(custom,par):
    return '''
    grey title="" pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['dratio'],par['dheight'],par['scalebar'],
           par['labelattr']+custom)

# wiggle 2D data
def dwigl2d(custom,par):
    return '''
    wiggle title="" pclip=100
    transp=y yreverse=y wherexlabel=t poly=y seamean=n
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['dratio'],par['dheight'],par['scalebar'],
           par['labelattr']+custom)

def egrey2d(custom,par):
    return '''
    grey title="" pclip=100
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['dratio'],par['dheight'],par['scalebar'],
           par['labelattr']+custom)

def ewigl2d(custom,par):
    return '''
    transp |
    wiggle title="" pclip=100
    transp=n yreverse=y wherexlabel=t poly=y seamean=n
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['dratio'],par['dheight'],par['scalebar'],
           par['labelattr']+custom)

# ------------------------------------------------------------
# plot wavelet
def waveplot(custom,par):
    return '''
    graph title=""
    min2=-1 max2=+1 
    plotfat=5 plotcol=5
    label1=%s unit1=%s
    label2="" unit2=""
    screenratio=0.5 screenht=7
    %s
    ''' % (par['lt'],par['ut'],
           par['labelattr']+custom)

# ------------------------------------------------------------
def cgraph2d(custom,par):
    return '''
    graph title=""
    labelrot=n wantaxis=n yreverse=y wherexlabel=t
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['iratio'],par['iheight'],par['scalebar'],
           par['labelattr']+custom)

def bbplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex | window |
    ''' + cgraph2d('plotcol=6 plotfat=2 %s'%custom,par)

def ssplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=6 plotfat=5 %s'%custom,par)

def rrplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=3 plotfat=5 %s'%custom,par)

def qqplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=1 plotfat=5 %s'%custom,par)

# ------------------------------------------------------------
def gainall(custom,par):
    return '''
    byte gainpanel=a pclip=100 %s
    '''%custom


