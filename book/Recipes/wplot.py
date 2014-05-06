try:    from rsf.cluster import *
except: from rsf.proj    import *
import pplot

# reset default parameters
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

    if par['nt']>1:
        par['df']=0.5/(par['ot']+(par['nt']-1)*par['dt'])
    else: 
        par['df']=1.
    par['of']=0
    par['nf']=par['nt']
    if(not par.has_key('lf')):       par['lf']='f'
    if(not par.has_key('uf')):       par['uf']='Hz'
    if(not par.has_key('fmin')):     par['fmin']=par['of']
    if(not par.has_key('fmax')):     par['fmax']=par['of'] + (par['nf']-1) * par['df']
    
    # make room to plot acquisition
    if(not par.has_key('zmin')):     par['zmin']=min(par['zmin'],-0.025*(par['zmax']-par['zmin']))

    dx=par['xmax']-par['xmin'];
    dy=par['ymax']-par['ymin'];
    dz=par['zmax']-par['zmin'];
    dt=par['tmax']-par['tmin'];
    
    if(not par.has_key('iratio2d')):
        if(dx==0.0): par['iratio2d']=1.0
        else:        par['iratio2d']=1.0*(dz)/(dx)
    if(not par.has_key('iheight2d')):
        if(par['iratio2d']>=0.8): par['iheight2d']=10
        else:                     par['iheight2d']=12*par['iratio2d']
            
    if(not par.has_key('dratio2d')):
        par['dratio2d']=par['iratio2d']

    if(not par.has_key('dheight2d')):
#       par['dheight2d']=par['iheight2d']
       par['dheight2d']=12*par['dratio2d']

    if((dx+dy) == 0.0)  : yxratio=1.0
    else                : yxratio=1.0*dx/(dx+dy)
    if((dz+dy) == 0.0)  : yzratio=1.0
    else                : yzratio=1.0*dz/(dz+dy)
    if((dt+dy) == 0.0)  : ytratio=1.0
    else                : ytratio=dt/(dt+dy);
    
    par['pointt']=ytratio;
    par['pointz']=yzratio;
    par['pointx']=yxratio;

    if((dx+dy) == 0.0): par['iratio3d']=1
    else:               par['iratio3d']=(dz+dy)/(dx+dy)

    if(par['iratio3d']>1): par['iheight3d']=12
    else:                  par['iheight3d']=12*par['iratio3d']

    if((dx+dy) == 0.0): par['dratio3d']=1
    else:               par['dratio3d']=(dt+dy)/(dx+dy)
        
    if(par['dratio3d']>1): par['dheight3d']=12
    else:                  par['dheight3d']=12*par['dratio3d']
        
    if(not par.has_key('scalebar')): par['scalebar']='n'    
    if(not par.has_key('labelattr')): par['labelattr']=' parallel2=n labelsz=8 labelfat=5 titlesz=12 titlefat=3 xll=2.5 yll=1.5 ' + ' '
    
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
           par['iratio2d'],par['iheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

# grey 2D frame of a cube
def ifrm2d(index,custom,par):
    return '''
    window n3=1 f3=%d |
    grey title=""
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (index,
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['iratio2d'],par['iheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

def igrey3d(custom,par):
    return '''
    grey3 title="" framelabel=n parallel2=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    xll=1.5 yll=1.5
    %s
    '''%(par['lz'],par['uz'],
         par['lx'],par['ux'],
         par['ly'],par['uy'],
         par['nz']/2,par['nx']/2,par['ny']/2,
         par['iratio3d'],par['iheight3d'],par['pointz'],par['pointx'],
         par['labelattr']+' '+custom)

def imovie3d(movie,byte,nfrm,custom,par):
    for ifrm in range(nfrm):
        ftag="-f%03d"%ifrm
        Plot(movie+ftag,movie+'byt',
             'window n4=1 f4=%d |'%ifrm + igrey3d(custom,par))
    Result(movie,[movie+'-f%03d'%ifrm for ifrm in range(nfrm)],'Movie')

def ifrmE2d(wfrm,wbyt,index,custom,par,xscale=0.5,yscale=0.5,shift=-11):
    Plot(wfrm+'_V',wbyt,'window n3=1 f3=0 |'+ ifrm2d(index,'',par))
    Plot(wfrm+'_H',wbyt,'window n3=1 f3=1 |'+ ifrm2d(index,'',par)) 
    pplot.p1x2(wfrm,wfrm+'_V',wfrm+'_H',xscale,yscale,shift)

def iovlE2d(out,inp,par,xscale=0.5,yscale=0.5,shift=-11):
    Plot(out+'_V',inp,'Overlay')
    Plot(out+'_H',inp,'Overlay')
    pplot.p1x2(out,out+'_V',out+'_H',xscale,yscale,shift)

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
           par['dratio2d'],par['dheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

def dgrey3d_init(vel,par):
    dx=par['xmax']-par['xmin'];
    dy=par['ymax']-par['ymin'];
    dt=par['tmax']-par['tmin'];
    dz=dt*vel;
    
    if((dx+dy) == 0.0): par['pointx']=1.0
    else              : par['pointx']=dx/(dx+dy)
    if((dt+dy) == 0.0): par['pointt']=1.0
    else              : par['pointt']=dz/(dz+dy);
    if((dx+dy) == 0.0): par['dratio3d']=1
    else:               par['dratio3d']=(dz+dy)/(dx+dy)
    if(par['dratio3d']>1): par['dheight3d']=11
    else:                  par['dheight3d']=11*par['dratio3d']

def dgrey3d(custom,par):
    return '''
    grey3 title="" framelabel=n parallel2=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    xll=1.5 yll=1.5
    %s
    '''%(par['lt'],par['ut'],
         par['lx'],par['ux'],
         par['ly'],par['uy'],
         par['nt']/2,par['nx']/2,par['ny']/2,
         par['dratio3d'],par['dheight3d'],par['pointt'],par['pointx'],
         par['labelattr']+' '+custom)

def dgreyE2d(data,dbyt,custom,par,xscale=0.5,yscale=0.5,shift=-11):
    Plot(data+'_V',dbyt,'window n2=1 f2=0 | transp |'+ dgrey2d('',par))
    Plot(data+'_H',dbyt,'window n2=1 f2=1 | transp |'+ dgrey2d('',par)) 
    pplot.p1x2(data,data+'_V',data+'_H',xscale,yscale,shift)

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
           par['dratio2d'],par['dheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

def egrey2d(custom,par):
    return '''
    grey title="" pclip=100
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['dratio2d'],par['dheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

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
           par['dratio2d'],par['dheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
def gainall(custom,par):
    return '''
    byte gainpanel=a pclip=100 %s
    '''%custom

def frame2d(frame,movie,index,custom,par):
    Flow(movie+'_p',movie,gainall(custom,par))
    Result(frame,movie+'_p',
           'window n3=1 f3=%d |'%index
           + igrey2d(custom,par))
    
# ------------------------------------------------------------
# plot wavelet
def waveplot(custom,par):
    return '''
    graph title=""
    min1=%g max1=%g
    min2=-1 max2=+1 
    plotfat=8 plotcol=5
    label1=%s unit1=%s
    label2="" unit2=""
    screenratio=0.3 screenht=4
    %s
    ''' % (par['tmin'],par['tmax'],
           par['lt'],par['ut'],
        par['labelattr']+' '+custom)

# plot spectrum
def specplot(custom,par):
    return '''
    scale axis=123 |
    graph title=""
    min1=%g max1=%g
    min2=0 max2=+1 
    plotfat=5 plotcol=5
    label1=%s unit1=%s
    label2="" unit2=""
    screenratio=0.3 screenht=4
    %s
    ''' % (par['fmin'],par['fmax'],
           par['lf'],par['uf'],
        par['labelattr']+' '+custom)

def waveplotE2d(wav,custom,par):

     Plot(wav+'_V',wav,
          'window n2=1 f2=0 | transp | window |'+
          waveplot(custom,par))
     Plot(wav+'_H',wav,
          'window n2=1 f2=1 | transp | window |'+
          waveplot(custom,par))

     pplot.p1x2(wav,wav+'_V',wav+'_H',0.5,0.5,-11.5)
     Result(wav,wav,'Overlay')

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
           par['iratio2d'],par['iheight2d'],par['scalebar'],
           par['labelattr']+' '+custom)

def bbplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex | window |
    ''' + cgraph2d('plotcol=6 plotfat=2 %s'%custom,par)

def ssplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=6 plotfat=15 %s'%custom,par)

def rrplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=3 plotfat=5 %s'%custom,par)

def qqplot2d(custom,par):
    return '''
    window n1=2 | dd type=complex |
    ''' + cgraph2d('symbol=. plotcol=1 plotfat=5 %s'%custom,par)

# ------------------------------------------------------------
# rays plot
def rays2d(plot,hwt,fray,jray,custom,par):
    Plot(plot,hwt,
         'window squeeze=n f1=%d j1=%d | transp |'%(fray,jray)
         + cgraph2d('plotcol=5 wantaxis=n'+' '+custom,par))

# wfts plot
def wfts2d(plot,hwt,fwft,jwft,custom,par):
    Plot(plot,hwt,
         'window squeeze=n f2=%d j2=%d | transp |'%(fwft,jwft)
         + cgraph2d('plotcol=2 wantaxis=n plotfat=2 symbol=.'+' '+custom,par))

# contour plot
def ccont2d(custom,par):
    return '''
    contour labelrot=n wantaxis=n title=""
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=n
    plotcol=2 plotfat=2
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['iratio2d'],par['iheight2d'],
        par['labelattr']+' '+custom)

# ------------------------------------------------------------
# ------------------------------------------------------------
# wavefield-over-model plot
def wom2d(wom,wfld,velo,vmean,nfrm,weight,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    wtmp = wfld + 'tmp'
    vtmp = wfld + 'vel'

    Flow(wom,[velo,wfld],
        '''
        %sscale < ${SOURCES[1]} axis=123 >%s datapath=%s/;
        '''%(M8R,wtmp,DPT) 
        +
        '''
        %sadd < ${SOURCES[0]} add=-%g |
        scale axis=123 |
        spray axis=3 n=%d o=%g d=%g
        >%s datapath=%s/;
        '''%(M8R,vmean,nfrm,0,1,vtmp,DPT) 
        +
        '''
        %sadd scale=1,%g <%s %s >${TARGETS[0]};
        '''%(M8R,weight,vtmp,wtmp) 
        +
        '''
        %srm %s %s
        '''%(M8R,wtmp,vtmp),
        stdin=0,
        stdout=0)


def ovl2d(wom,wfld,velo,vmean,nfrm,weight,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    wtmp = wfld + 'tmp'
    vtmp = wfld + 'vel'

    Flow(wom,[velo,wfld],
        '''
        %sscale < ${SOURCES[1]} axis=123 >%s datapath=%s/;
        '''%(M8R,wtmp,DPT) 
        +
        '''
        %sadd < ${SOURCES[0]} add=-%g |
        scale axis=123 >%s datapath=%s/;
        '''%(M8R,vmean,vtmp,DPT) 
        +
        '''
        %sadd scale=1,%g <%s %s >${TARGETS[0]};
        '''%(M8R,weight,vtmp,wtmp) 
        +
        '''
        %srm %s %s
        '''%(M8R,wtmp,vtmp),
        stdin=0,
        stdout=0)
    
# ------------------------------------------------------------
# static movie = tightly packed, side-by-side frames
def istagger2d(cube,custom,par,nfrm=2,scale=1,ratio=1,ymax=10,xmax=14):
    oy=-ymax*(1./scale-1); dy=abs(oy)/(nfrm-1)
    ox=0;                  dx=xmax*(1./scale-1)/(nfrm-1)

    Flow(cube+'byt',cube,'byte gainpanel=a pclip=100 %s'%custom)    
    for ifrm in range(nfrm):
        tag="%04d"%ifrm
        
        Plot(cube+tag,cube+'byt',
             'window n3=1 f3=%d |'%ifrm + igrey2d('titlesz=%d title=%d %s'%(4/scale,ifrm,custom),par))
        Plot(cube+tag+'_',
             cube+tag,'Overlay',vppen='yscale=%f xscale=%f xcenter=%f ycenter=%f'
             %(scale,scale,ox-ifrm*dx,oy+ifrm*dy))        
    Result(cube,[cube+"%04d_"%ifrm for ifrm in range(nfrm)],'Overlay')

# ------------------------------------------------------------
# plot 3x3 matrix of plots
def inine(cube,byte,custom,par,scale=0.3,ymax=10,ratio=1):
    nfrm=9

    dy=ymax;
    dx=ymax/ratio
        
    for ifrm in range(nfrm):
        tag="%d"%ifrm
        
        Plot(cube+tag,byte,
             'window n3=1 f3=%d |'%ifrm 
             + igrey2d('wantaxis=n titlesz=%d title=%d %s'%(4/scale,ifrm,custom),par))
        Plot(cube+tag+'_',
             cube+tag,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'
             %(scale,scale,-1+int(-2+ifrm/3)*dy,-1-(ifrm%3)*dx))        
    Result(cube,[cube+"%d_"%ifrm for ifrm in range(nfrm)],'Overlay')

# ------------------------------------------------------------
# plot 1x5 vector of plots (assumes 9 frame in input)
def ifive(cube,byte,custom,par,scale=0.2,xmax=12):
    nfrm=5
    jfrm=2
    
    dx=xmax

    #    Flow(cube+'byt',cube,'byte gainpanel=a pclip=100 %s'%custom)    
    for ifrm in range(nfrm):
        tag="%d"%ifrm
        
        Plot(cube+tag,byte,
             'window n3=1 f3=%d |'%(ifrm*jfrm)
             + igrey2d('wantaxis=n titlesz=%d title=%d %s'%(4/scale,ifrm,custom),par))
        Plot(cube+tag+'_',
             cube+tag,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'
             %(scale,scale,-1,-1-(ifrm%5)*dx))        
    Result(cube,[cube+"%d_"%ifrm for ifrm in range(nfrm)],'Overlay')

# ------------------------------------------------------------
# plot 5x1 vector of plots (assumes 9 frame in input)
def cfive(cube,byte,custom,par,scale=0.2,ymax=9):
    nfrm=5
    jfrm=2
    
    dy=ymax

    for ifrm in range(nfrm):
        tag="%d"%ifrm
        
        Plot(cube+tag,byte,
             'window n3=1 f3=%d |'%(ifrm*jfrm)
             + igrey2d('wantaxis=n titlesz=%d title=%d %s'%(4/scale,ifrm,custom),par))
        Plot(cube+tag+'_',
             cube+tag,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=%f'
             %(scale,scale,-2-4*dy+(ifrm%5)*dy,-2))        
    Result(cube,[cube+"%d_"%ifrm for ifrm in range(nfrm)],'Overlay')
