try:    from rsf.cluster import *
except: from rsf.proj    import *
import pplot,math

# ------------------------------------------------------------
def obsolete(old,new):
    print "------------------------------------------------------------------"
    print "the function",old,"is OBSOLETE; use",new
    print "------------------------------------------------------------------"

# ------------------------------------------------------------
def Temp(o,i,r):
    Flow(o,i,r+ ' datapath=%s '%os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH')))
# ------------------------------------------------------------

# default parameters
def param(par):
    if(not par.has_key('verb')):     par['verb']='n'
    
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

    if(not par.has_key('nbx')):      par['nbx']=0
    if(not par.has_key('nby')):      par['nby']=0
    if(not par.has_key('nbz')):      par['nbz']=0
    if(not par.has_key('tx')):       par['tx']=0.0035
    if(not par.has_key('ty')):       par['ty']=0.0035
    if(not par.has_key('tz')):       par['tz']=0.0035

    if(not par.has_key('verb')):     par['verb']='n'
    if(not par.has_key('nb')):       par['nb']=0
    if(not par.has_key('nbell')):    par['nbell']=5
    if(not par.has_key('snap')):     par['snap']='y'
    if(not par.has_key('jsnap')):    par['jsnap']=100
    if(not par.has_key('jdata')):    par['jdata']=1
    if(not par.has_key('dabc')):     par['dabc']='y'
    if(not par.has_key('ompchunk')): par['ompchunk']=1
    if(not par.has_key('ompnth')):   par['ompnth']=0
    if(not par.has_key('free')):     par['free']='n'
    if(not par.has_key('fsrf')):     par['fsrf']='n'
   
    if(not par.has_key('ratio')):
        if(dx==0.0):
            par['ratio']=1.0
        else:
            par['ratio']=1.0*(dz)/(dx)

    if(not par.has_key('height')):
        if(par['ratio']>=1):
            par['height']=10
        else:
            par['height']=13.625*par['ratio']

    if(not par.has_key('dratio')):
        if(dx==0.0):
            par['dratio']=1.0
        else:
            par['dratio']=1.0*(dt)/(dx)

    if(not par.has_key('dheight')):
        if(par['dratio']>1):
            par['dheight']=10
        else:
            par['dheight']=13*par['dratio']

    if((dx+dy)   == 0.0)  : yxratio=1.0
    else                  : yxratio=1.0*dx/(dx+dy)
    if((dz+dy)   == 0.0)  : yzratio=1.0
    else                  : yzratio=1.0*dz/(dz+dy)
    if((2*dt+dy) == 0.0)  : ytratio=1.0
    else                  : ytratio=2*dt/(2*dt+dy);
    
    par['pointt']=ytratio;
    par['pointz']=yzratio;
    par['pointx']=yxratio;
    
    if((dx+dy) == 0.0):
        par['ratio3d']=1
    else:
        par['ratio3d']=(dz+dy)/(dx+dy)
        
    if(par['ratio3d']>1):
        par['height3d']=10
    else:
        par['height3d']=14*par['ratio3d']
        
    if((dx+dy) == 0.0):
        par['tratio3d']=1
    else:
        par['tratio3d']=(2*dt+dy)/(dx+dy)
        
    if(par['tratio3d']>1):
        par['theight3d']=10
    else:
        par['theight3d']=11*par['tratio3d']

    if(not par.has_key('scalebar')): par['scalebar']='n'
    if(not par.has_key('labelattr')): par['labelattr']=' parallel2=n labelsz=6 labelfat=3 titlesz=12 titlefat=3 '
    
    if(not par.has_key('nqz')): par['nqz']=par['nz']
    if(not par.has_key('oqz')): par['oqz']=par['oz']
    if(not par.has_key('dqz')): par['dqz']=par['dz']

    if(not par.has_key('nqx')): par['nqx']=par['nx']
    if(not par.has_key('oqx')): par['oqx']=par['ox']
    if(not par.has_key('dqx')): par['dqx']=par['dx']

    par['xratio']=2
    par['yratio']=2
    par['tratio']=2
    par['aratio']=2

    par['labelrot'] =' parallel2=n '
    par['labelrot0']=' parallel2=n format1=%3.0f format2=%3.0f format3=%3.0f '
    par['labelrot1']=' parallel2=n format1=%3.1f format2=%3.1f format3=%3.1f '
    par['labelrot2']=' parallel2=n format1=%3.2f format2=%3.2f format3=%3.2f '


# ------------------------------------------------------------
def modpar(par):
    default(par)

def default(par):
    if(not par.has_key('nbx')):      par['nbx']=0
    if(not par.has_key('nby')):      par['nby']=0
    if(not par.has_key('nbz')):      par['nbz']=0
    if(not par.has_key('tx')):       par['tx']=0.0035
    if(not par.has_key('ty')):       par['ty']=0.0035
    if(not par.has_key('tz')):       par['tz']=0.0035

    if(not par.has_key('nb')):       par['nb']=0
    if(not par.has_key('nbell')):    par['nbell']=5
    if(not par.has_key('snap')):     par['snap']='y'
    if(not par.has_key('jsnap')):    par['jsnap']=100
    if(not par.has_key('jdata')):    par['jdata']=1
    if(not par.has_key('dabc')):     par['dabc']='y'
    if(not par.has_key('ompchunk')): par['ompchunk']=1
    if(not par.has_key('ompnth')):   par['ompnth']=0
    if(not par.has_key('free')):     par['free']='n'

    if(not par.has_key('nqz')):      par['nqz']=par['nz']
    if(not par.has_key('oqz')):      par['oqz']=par['oz']
    if(not par.has_key('dqz')):      par['dqz']=par['dz']

    if(not par.has_key('nqx')):      par['nqx']=par['nx']
    if(not par.has_key('oqx')):      par['oqx']=par['ox']
    if(not par.has_key('dqx')):      par['dqx']=par['dx']
    
# ------------------------------------------------------------
# plotting functions
def cgrey(custom,par):
    return '''
    grey 
    title=""
    pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g wantscalebar=%s
    %s
    ''' % (par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],par['scalebar'],
           par['labelattr']+custom)

def ccut3d(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g min3=%g max3=%g %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['ymin'],par['ymax'],
	   custom)

#    byte gainpanel=a pclip=100 %s |
def cgrey3d(custom,par):
    return '''
    grey3 title="" framelabel=n parallel2=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    xll=1.5 yll=1.5
    %s
    ''' % (
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['nz']/2,par['nx']/2,par['ny']/2,
           par['ratio3d'],par['height3d'],par['pointz'],par['pointx'],
           par['labelattr']+' '+custom)

def wgrey(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g |
    grey parallel2=n labelrot=n wantaxis=y title=""
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
    graph labelrot=n wantaxis=n title="" yreverse=y wherexlabel=t
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
    grey parallel2=n labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=0.5 screenht=7
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['labelattr']+' '+custom)

def dwigl(custom,par):
    return '''
    wiggle parallel2=n labelrot=n wantaxis=y title=""
    transp=y yreverse=y wherexlabel=t poly=y seamean=n
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=0.5 screenht=7
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['labelattr']+' '+custom)

def dgrey3d(custom,par):
    return '''
    window min2=%g max2=%g min3=%g max3=%g |
    byte gainpanel=a pclip=100 %s |
    grey3 title="" framelabel=n parallel2=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    xll=1.5 yll=1.5
    %s
    ''' % (
           par['xmin'],par['xmax'],
           par['ymin'],par['ymax'],
           custom,
           par['lt'],par['ut'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['nt']/2,par['nx']/2,par['ny']/2,
           par['tratio3d'],par['theight3d'],par['pointt'],par['pointx'],
           par['labelattr']+' '+custom)

def egrey(custom,par):
    return '''
    grey parallel2=n labelrot=n wantaxis=y title=""
    pclip=100
    min2=%g max2=%g label2=%s unit2=%s
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=%g screenht=%g
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['ratio'],par['height'],
           par['labelattr']+' '+custom)

def ewigl(custom,par):
    return '''
    transp |
    wiggle parallel2=n labelrot=n wantaxis=y title=""
    transp=n yreverse=y wherexlabel=t poly=y seamean=n
    pclip=100
    min1=%g max1=%g label1=%s unit1=%s
    min2=%g max2=%g label2=%s unit2=%s
    screenratio=%g screenht=%g
    %s
    ''' % (par['tmin'],par['tmax'],par['lt'],par['ut'],
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['ratio'],par['height'],
           par['labelattr']+' '+custom)


def fgrey(custom,par):
    return '''
    window | real | transp |
    grey parallel2=n labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    min2=%g max2=%g label2=%s unit2=%s
    label1="f" unit1="Hz"
    screenratio=%g screenht=%g
    %s
    ''' % (par['xmin'],par['xmax'],par['lx'],par['ux'],
           par['ratio'],par['height'],
           par['labelattr']+' '+custom)

# ------------------------------------------------------------
def center3d(x,y,z,par):
    return '''
    frame1=%d frame2=%d frame3=%d
    ''' % ( (float(z)-par['oz'])/par['dz']+1,
            (float(x)-par['ox'])/par['dx']+1,
            (float(y)-par['oy'])/par['dy']+1 )

# ------------------------------------------------------------
# plot wavelet
def waveplot(custom,par):
    return '''
    graph title="" min1=%g min2=-1 max2=+1
    plotfat=5 plotcol=5
    label1=%s unit1=%s
    label2="" unit2=""
    parallel2=n screenratio=0.5 screenht=7
    %s
    ''' % (par['ot'],par['lt'],par['ut'],
           par['labelattr']+' '+custom)

def spectrum(custom,par):
    return '''
    spectra |
    scale axis=123 |
    graph title="" plotfat=5 plotcol=5
    label1="f" unit1="Hz"
    min2=0 max2=1 label2=""
    %s
    ''' %(par['labelattr']+' '+custom)


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
def oldhorizontal(cc,coord,par):
    Temp(cc+'_',None,'math n1=%(nx)d d1=%(dx)g o1=%(ox)g output=0' % par)
    Temp(cc+'_z',cc+'_','math output="%g" ' % coord)
    Temp(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n ${SOURCES[1]} |
         transp |
         put label1="" unit1="" label2="" unit2=""
         ''')

def horizontal(cc,coord,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))
    
    cco=cc+'o'
    ccz=cc+'z'
    ccx=cc+'x'
    
    Flow(cc,None,
         '''
         %smath output=0 n1=%d o1=%g d1=%g >%s datapath=%s/;
         '''%(M8R,par['nx'],par['ox'],par['dx'],cco,DPT) +
         '''
         %smath <%s output="%g" >%s datapath=%s/;
         '''%(M8R,cco,coord,ccz,DPT) +
         '''
         %smath <%s output="x1" >%s datapath=%s/;
         '''%(M8R,cco,ccx,DPT) +
         '''
         %scat axis=2 space=n %s %s | transp | put label1="" unit1="" label2="" unit2="">${TARGETS[0]};
         '''%(M8R,ccx,ccz) +
         '''     
         %srm %s %s %s
         '''%(M8R,cco,ccx,ccz),
              stdin=0,
              stdout=0)


def horizontalupercent(cc,coord,par):
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))
    
    cco=cc+'o'
    ccz=cc+'z'
    ccx=cc+'x'
    
    Flow(cc,None,
         '''  
         math output=0 n1=%d o1=%g d1=%g >%s datapath=%s/ &&
         '''%(par['nx'],par['ox'],par['dx'],cco,DPT) +
         '''  
         math <%s output="%g" >%s datapath=%s/ &&
         '''%(cco,coord,ccz,DPT) +
         '''  
         math <%s output="x1" >%s datapath=%s/ &&
         '''%(cco,ccx,DPT) +
         '''  
         cat axis=2 space=n %s %s | transp | put label1="" unit1="" label2="" unit2="">${TARGETS[0]} &&
         '''%(ccx,ccz) +
         '''  
         rm %s %s %s
         '''%(cco,ccx,ccz),
              stdin=0,
              stdout=0)




# ------------------------------------------------------------
def cable2d(cc,zrec,orec,nrec,drec,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))
    
    cco=cc+'o'
    ccz=cc+'z'
    ccx=cc+'x'
    
    Flow(cc,None,
         '''
         %smath output=0 n1=%d o1=%g d1=%g >%s datapath=%s/;
         '''%(M8R,nrec,orec,drec,cco,DPT) +
         '''
         %smath <%s output="%g" >%s datapath=%s/;
         '''%(M8R,cco,zrec,ccz,DPT) +
         '''
         %smath <%s output="x1" >%s datapath=%s/;
         '''%(M8R,cco,ccx,DPT) +
         '''
         %scat axis=2 space=n %s %s | transp | put label1="" unit1="" label2="" unit2="">${TARGETS[0]};
         '''%(M8R,ccx,ccz) +
         '''     
         %srm %s %s %s
         '''%(M8R,cco,ccx,ccz),
              stdin=0,
              stdout=0)


# ------------------------------------------------------------
def horizontal3d(cc,coord,par):
    Temp(cc+'_',None,
         'math n1=%(nx)d d1=%(dx)g o1=%(ox)g n2=%(ny)d d2=%(dy)g o2=%(oy)g output=0' % par)
    Temp(cc+'_z',cc+'_','math output="%g" | put n1=%d n2=%d n3=1' % (coord,par['nx'],par['ny']) )

    if(par['nx']==1):
        Temp(cc+'_x',cc+'_',
             'math output="%g" | put n1=%d n2=%d n3=1' % (par['ox'],par['nx'],par['ny']) )
    else:
        Temp(cc+'_x',cc+'_',
             'math output="x1" | put n1=%d n2=%d n3=1' % (      par['nx'],par['ny']) )

    if(par['ny']==1):
        Temp(cc+'_y',cc+'_',
             'math output="%g" | put n1=%d n2=%d n3=1' % (par['oy'],par['nx'],par['ny']) )
    else:
        Temp(cc+'_y',cc+'_',
             'math output="x2" | put n1=%d n2=%d n3=1' % (          par['nx'],par['ny']) )

    Flow(cc,[cc+'_x',cc+'_y',cc+'_z'],
         '''
         cat axis=3 space=n ${SOURCES[1:3]} | 
         transp plane=13 | transp plane=23 |
	     put label1="" unit1="" label2="" unit2=""
         ''')

def oldvertical(cc,coord,par):
    Temp(cc+'_',None,'math n1=%(nz)d d1=%(dz)g o1=%(oz)g output=0' % par)
    Temp(cc+'_x',cc+'_','math output="%g" '% coord)
    Temp(cc+'_z',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n ${SOURCES[1]} | 
         transp |
         put label1="" unit1="" label2="" unit2=""
         ''')
    
def vertical(cc,coord,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))
    
    cco=cc+'o'
    ccz=cc+'z'
    ccx=cc+'x'
    
    Flow(cc,None,
         '''
         %smath output=0 n1=%d o1=%g d1=%g >%s datapath=%s/;
         '''%(M8R,par['nz'],par['oz'],par['dz'],cco,DPT) +
         '''
         %smath <%s output="x1" >%s datapath=%s/;
         '''%(M8R,cco,ccz,DPT) +
         '''
         %smath <%s output="%g" >%s datapath=%s/;
         '''%(M8R,cco,coord,ccx,DPT) +
         '''
         %scat axis=2 space=n %s %s | transp | put label1="" unit1="" label2="" unit2="">${TARGETS[0]};
         '''%(M8R,ccx,ccz) +
         '''     
         %srm %s %s %s
         '''%(M8R,cco,ccx,ccz),
              stdin=0,
              stdout=0)
    

def vertical3d(cc,coordx,coordy,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    cco=cc+'o'
    ccx=cc+'x'
    ccy=cc+'y'
    ccz=cc+'z'

    Flow(cc,None,
         '''
         %smath n1=%d d1=%g o1=%g output=0 >%s datapath=%s/;
         '''%(M8R,par['nz'],par['dz'],par['oz'],cco,DPT) +
         '''
         %smath <%s output="%g" >%s datapath=%s/;
         '''%(M8R,cco,coordx,ccx,DPT) +
         '''
         %smath <%s output="%g" >%s datapath=%s/;
         '''%(M8R,cco,coordy,ccy,DPT) +
         '''
         %smath <%s output="x1" >%s datapath=%s/;
         '''%(M8R,cco,ccz,DPT) +
         '''
         %scat axis=2 space=n %s %s %s |
         transp |
	 put label1="" unit1="" label2="" unit2="" >${TARGETS[0]};
         '''%(M8R,ccx,ccy,ccz) +
         '''
         %srm %s %s %s %s
         '''%(M8R,cco,ccx,ccy,ccz),
         stdin=0,
         stdout=0
        )
    
def point(cc,xcoord,zcoord,par):
    Flow(cc,None,
         'spike nsp=2 mag=%g,%g n1=2 k1=1,2'%(xcoord,zcoord))
    
def point3d(cc,xcoord,ycoord,zcoord,par):
    Flow(cc,None,
         'spike nsp=3 mag=%g,%g,%g n1=3 k1=1,2,3'%(xcoord,ycoord,zcoord))
    
def point3(cc,xcoord,zcoord,magn,par):
    Flow(cc,None,
         'spike nsp=3 mag=%g,%g,%g n1=3 k1=1,2,3'%(xcoord,zcoord,magn))

# ------------------------------------------------------------
def circle(cc,xcenter,zcenter,radius,sampling,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))
    
    ccx=cc+'x'
    ccz=cc+'z'

    Flow(cc,None,
         '''
         %smath n1=%d d1=%g o1=%g output="%g+%g*cos(%g*x1/180.)" >%s datapath=%s/;
         '''%(M8R,sampling,360./sampling,0.,xcenter,radius,math.pi,ccx,DPT) +
         '''
         %smath n1=%d d1=%g o1=%g output="%g-%g*sin(%g*x1/180)" >%s datapath=%s/;
         '''%(M8R,sampling,360./sampling,0.,zcenter,radius,math.pi,ccz,DPT) +
         '''
         %scat axis=2 space=n %s %s |
         transp |
         put label1="" unit1="" label2="" unit2="" >${TARGETS[0]};
         '''%(M8R,ccx,ccz) +
         '''
         %srm %s %s
         '''%(M8R,ccx,ccz),
         stdin=0,
         stdout=0
        )
    
# ------------------------------------------------------------
def ellipse(cc,xcenter,zcenter,semiA,semiB,sampling,par):
    Temp(cc+'_r',None,
         '''
         math n1=%d d1=%g o1=%g
         output="((%g)*(%g))/sqrt( ((%g)*cos(%g*x1/180))^2 + ((%g)*sin(%g*x1/180))^2 )"
         '''% (sampling,360./sampling,0.,
                semiA,semiB,semiB,math.pi,semiA,math.pi) )
    Temp(cc+'_x',cc+'_r',
         'math output="%g+input*cos(%g*x1/180)"'
         % (xcenter,math.pi) )
    Temp(cc+'_z',cc+'_r',
         'math output="%g+input*sin(%g*x1/180)"'
         % (zcenter,math.pi) )

    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp |
         put label1="" unit1="" label2="" unit2=""
         ''', stdin=0)
    
def dipping(cc,intercept,slope,par):
    Temp(cc+'_',None,'math n1=%(nx)d d1=%(dx)g o1=%(ox)g output=0' % par)
    Temp(cc+'_z',cc+'_','math output="%g+x1*%g" '%(intercept,slope))
    Temp(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp |
	 put label1="" unit1="" label2="" unit2=""
         ''', stdin=0)

def boxarray(cc,nz,oz,dz,nx,ox,dx,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    cco=cc+'o'
    ccz=cc+'z'
    ccx=cc+'x'

    Flow(cc,None,
         '''
         %smath output=1 n1=%d o1=%g d1=%g n2=%d o2=%g d2=%g >%s datapath=%s/;
         '''%(M8R,nz,oz,dz,nx,ox,dx,cco,DPT) +
         '''
         %smath <%s output="x1" | put n1=%d n2=1 >%s datapath=%s/;
         '''%(M8R,cco,nz*nx,ccz,DPT) +
         '''
         %smath <%s output="x2" | put n1=%d n2=1 >%s datapath=%s/;
         '''%(M8R,cco,nz*nx,ccx,DPT) +
         '''
         %scat axis=2 space=n %s %s | transp >${TARGETS[0]};
         '''%(M8R,ccx,ccz) +
         '''     
         %srm %s %s %s
         '''%(M8R,cco,ccx,ccz),
              stdin=0,
              stdout=0)

def boxarray3d(cc,nz,oz,dz,nx,ox,dx,ny,oy,dy,par):
    Temp(cc+'_',None,
         '''
         math output=1
         n1=%d d1=%g o1=%g
         n2=%d d2=%g o2=%g
	 n3=%d d3=%g o3=%g
         ''' % (nz,dz,oz,
	 	nx,dx,ox,
		ny,dy,oy) )
    Temp(cc+'_z',cc+'_','math output="x1" | put n1=%d n2=1 n3=1' % (nz*nx*ny))
    Temp(cc+'_x',cc+'_','math output="x2" | put n1=%d n2=1 n3=1' % (nz*nx*ny))
    Temp(cc+'_y',cc+'_','math output="x3" | put n1=%d n2=1 n3=1' % (nz*nx*ny))
    Flow(cc,[cc+'_x',cc+'_y',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | transp |
         put label1="" unit1="" label2="" unit2="" label3="" unit3=""
         ''',stdin=0)

#
# distribute random locations within a defined box 
#
#       _________________  <-- zmin
#      | . .  .  .   .   | 
#      |    .   .   .  . |
#      | .         .   . |
#      |     .   .   .   |
#       _________________  <-- zmax
#      ^                 ^
#      |                 |
#     xmin              xmax
#    
# ns: number of random points
#
# seed: seed that start the pseudo-random number sequence. 
#       It is hardcoded on purpose, so the results can be 
#       reproduced   
#

def rand_boxarray(rsrc,xmin,xmax,zmin,zmax,ns,seed=33121):

    m8r='$RSFROOT/bin/sf'
    xrand='xrand'+rsrc 
    zrand='zrand'+rsrc 

    custom1=' seed=%d'%seed
    custom2=' seed=%d'%(seed*4+1)

    Flow(rsrc,None,
        '''
        spike n1=%d | noise type=n rep=y range=1 mean=0 %s |
	    add add=1 | scale axis=123 | add scale=%f |
        add add=%g > %s ;
        '''%(ns,custom1,(xmax-xmin),xmin,xrand)+
        '''
        %sspike n1=%d | noise type=n rep=y range=1 mean=0 %s |
	    add add=1 | scale axis=123 | add scale=%f |
        add add=%g > %s ;
        '''%(m8r,ns,custom2,(zmax-zmin),zmin,zrand)+
        '''
        %scat axis=2 %s %s |transp >${TARGETS[0]} ; 
        '''%(m8r,xrand,zrand)+
        '''
        %srm %s %s
        '''%(m8r,xrand,zrand),stdout=0)

def makebox(box,zmin,zmax,xmin,xmax,par):
    Temp(box+'_z',None,
         '''
         spike nsp=5 mag=%g,%g,%g,%g,%g
         n1=5 o1=0 d1=1 k1=1,2,3,4,5 |
         transp
         '''%(zmin,zmax,zmax,zmin,zmin))
    Temp(box+'_x',None,
         '''
         spike nsp=5 mag=%g,%g,%g,%g,%g
         n1=5 o1=0 d1=1 k1=1,2,3,4,5 |
         transp
         '''%(xmin,xmin,xmax,xmax,xmin))
    Flow(box,[box+'_x',box+'_z'],'cat axis=1 space=n ${SOURCES[1]}')

# ------------------------------------------------------------
def makeline(line,zmin,zmax,xmin,xmax,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    linez=line+'z'
    linex=line+'x'

    Flow(line,None,
         '''
         %sspike nsp=2 mag=%g,%g n1=2 o1=0 d1=1 k1=1,2 |
         transp >%s datapath=%s/;
         '''%(M8R,zmin,zmax,linez,DPT) +
         '''
         %sspike nsp=2 mag=%g,%g n1=2 o1=0 d1=1 k1=1,2 |
         transp >%s datapath=%s/;
         '''%(M8R,xmin,xmax,linex,DPT) +
         '''
         %scat axis=1 space=n %s %s >${TARGETS[0]};
         '''%(M8R,linex,linez) +
         '''     
         %srm %s %s
         '''%(M8R,linex,linez),
              stdin=0,
              stdout=0)

# ------------------------------------------------------------
def hline(cc,sx,ex,coord,par):
    
    nx=(ex-sx)/par['dx']+1
    dx=par['dx']

    Temp(cc+'_',None,'math n1=%d d1=%g o1=%g output=0' % (nx,dx,sx))
    Temp(cc+'_z',cc+'_','math output="%g" ' % coord)
    Temp(cc+'_x',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp |
	 put label1="" unit1="" label2="" unit2=""
         ''', stdin=0)

def vline(cc,sz,ez,coord,par):

    nz=(ez-sz)/par['dz']+1
    dz=par['dz']
    
    Temp(cc+'_',None,'math n1=%d d1=%g o1=%g output=0' % (nz,dz,sz))
    Temp(cc+'_x',cc+'_','math output="%g" ' % coord)
    Temp(cc+'_z',cc+'_','math output="x1" ')
    Flow(cc,[cc+'_x',cc+'_z'],
         '''
         cat axis=2 space=n
         ${SOURCES[0]} ${SOURCES[1]} | transp |
	 put label1="" unit1="" label2="" unit2=""
         ''', stdin=0)

def box(cc,sx,ex,sz,ez,par):
   
    hline(cc+'h1',sx,ex,sz,par)
    hline(cc+'h2',sx,ex,ez,par)

    vline(cc+'v1',sz,ez,sx,par)
    vline(cc+'v2',sz,ez,ex,par)

    Flow(cc,[cc+'h1',cc+'h2',cc+'v1',cc+'v2'],
         'cat ${SOURCES[1:4]} space=n axis=2 | put label1="" unit1="" label2="" unit2=""')

# ------------------------------------------------------------
def bbplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    window |
    ''' + cgraph('plotcol=6 plotfat=2 wantaxis=n %s' % custom,par)

def ssplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=6 plotfat=10 wantaxis=n %s' % custom,par)

def rrplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=3 plotfat=5 wantaxis=n %s' % custom,par)

def qqplot(custom,par):
    return '''
    window n1=2 |
    dd type=complex |
    ''' + cgraph('symbol=. plotcol=1 plotfat=5 wantaxis=n %s' % custom,par)

def qqbox2d(par):
    return '''
    nqz=%(nqz)d oqz=%(oqz)g
    nqx=%(nqx)d oqx=%(oqx)g
    ''' % par

def qqbox3d(par):
    return '''
    nqz=%(nqz)d oqz=%(oqz)g
    nqx=%(nqx)d oqx=%(oqx)g
    nqy=%(nqy)d oqy=%(oqy)g
    ''' % par

def iwindow(par):
    win = ' ' + \
          '''
          nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g 
          nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g
          ''' % par + ' '

    return win

# ------------------------------------------------------------
# rays plot
def rayplot(hwt,j1ray,j2ray,j1wft,j2wft,custom,par):

    Plot(hwt+'ray',hwt,'window squeeze=n j1=%d j2=%d f2=%d | transp |' %(j1ray,j2ray,j2wft)
         + cgraph('plotcol=5 wantaxis=n '+custom,par))
    Plot(hwt+'wft',hwt,'window j1=%d j2=%d f2=%d |'          %(j1wft,j2wft,j2wft)
         + cgraph('plotcol=2 squeeze=n wantaxis=n symbol=. '+custom,par))

    Plot  (hwt,[hwt+'ray',hwt+'wft'],'Overlay')

# ------------------------------------------------------------
# acoustic modeling
def awefd(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# acoustic modeling that only generates data 
# ------------------------------------------------------------
def awefd_data(odat,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat],[idat,velo,dens,sou,rec],
         '''
         awefd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}'''%par+'''
         wfl=%stmp  '''%odat+ '''
         %s >${TARGETS[0]}; 
         $RSFROOT/bin/sfrm %stmp
         ''' %(par['fdcustom'],odat),stdout=0)

# ------------------------------------------------------------
# acoustic modeling that only generates source wavefield 
# ------------------------------------------------------------
def awefd_swfl(owfl,idat,velo,dens,sou,custom,par):
    par['fdcustom'] = custom
    
    Flow([owfl],[idat,velo,dens,sou],
         '''
         awefd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         wfl=${TARGETS[0]} '''%par+'''
         %s >%stmp ; $RSFROOT/bin/sfrm %stmp
         ''' %(par['fdcustom'],owfl,owfl),stdout=0)
    
# ------------------------------------------------------------
# acoustic modeling that only generates receiver wavefield:
# the input data is not reversed, all the 
# secondary steps are taken in the rule, and all the 
# intermediate files are deleted
# ------------------------------------------------------------
def awefd_rwfl(owfl,idat,velo,dens,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([owfl],[idat,velo,dens,rec],
         '''
         $RSFROOT/bin/sfreverse opt=i which=2 |
         awefd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}'''%par+'''
         wfl=%s_junk
         %s >%stmp  ; $RSFROOT/bin/sfrm %stmp ;
         $RSFROOT/bin/sfreverse opt=i which=4 <%s_junk >${TARGETS[0]};
         $RSFROOT/bin/sfrm %s_junk
         ''' %(owfl,par['fdcustom'],owfl,owfl,owfl,owfl),stdout=0)

def awefd1(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    awefd(odat,owfl,idat,velo,dens,sou,rec,custom+' expl=y ',par)

# constant-density acoustic FD modeling
def cdafd(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd2d cden=y
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         sou=${SOURCES[2]}
         rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         '''%par+custom)
def cdafd1(odat,owfl,idat,velo,sou,rec,custom,par):
    cdafd(odat,owfl,idat,velo,sou,rec,custom+' expl=y ',par)

def cdafd3d(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd3d cden=y
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         dabc=%(dabc)s nb=%(nb)d
         vel=${SOURCES[1]}
         sou=${SOURCES[2]}
         rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         '''%par+custom)
    
# ------------------------------------------------------------
# Born modeling
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
# anisotropic stiffness tensor
#def anisotropic(cc,vp,vs,ro,epsilon,delta,par):
#    Flow(cc+'33',[vp,ro],
#         '''
#         math output="ro*vp*vp"
#         vp=${SOURCES[0]}
#         ro=${SOURCES[1]}
#         ''')    
#    Flow(cc+'44',[vs,ro],
#         '''
#         math output="ro*vs*vs"
#         vs=${SOURCES[0]}
#         ro=${SOURCES[1]}
#         ''')
#    Flow(cc+'11',[cc+'33',epsilon],
#         '''
#         math output="2*epsilon*c33+c33"
#         c33=${SOURCES[0]}
#         epsilon=${SOURCES[1]}
#         ''')
#    Flow(cc+'13',[cc+'33',cc+'44',delta],
#         '''
#         math output="sqrt(2*c33*(c33-c44)*delta+(c33-c44)*(c33-c44))-c44"
#         c33=${SOURCES[0]}
#         c44=${SOURCES[1]}
#         delta=${SOURCES[2]}
#         ''')
#    
#    Flow(cc,[cc+'11',cc+'13',cc+'33',cc+'44'],
#         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def animodel(v,vv,eta,delta,theta):
    Flow(v+'_vv',[vv],           'window')
    Flow(v+'_vn',[v+'_vv',delta],'math del=${SOURCES[1]} output="input*sqrt(1+2*del)"')
    Flow(v+'_vh',[v+'_vn',eta  ],'math eta=${SOURCES[1]} output="input*sqrt(1+2*eta)"')
    Flow(v,[v+'_vv',v+'_vn',v+'_vh',theta],'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# acoustic anisotropic modeling
def anifd2d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['anifdcustom'] = custom

    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         ttifd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d dabc=%(dabc)s
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(anifdcustom)s
         ''' % par)

# ------------------------------------------------------------
# elastic modeling
def ewefd(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd2d
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d jdata=%(jdata)d
	 nb=%(nb)d nbell=%(nbell)d
	 ssou=%(ssou)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd2(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd2dtti
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
# exploding-reflector reverse-time migration
def zom(imag,data,velo,dens,rcoo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    obsolete("fdmod/zom","rtm/zofmig")
    
    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom

    rdat = imag+'rdat'
    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo,dens],
         '''
         %sreverse < ${SOURCES[0]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdat,DPT) +
         '''
         %sawefd2d < %s cden=n %s verb=n
         vel=${SOURCES[2]}
         den=${SOURCES[3]}
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdat,awepar+' jsnap=%d'%(par['nt']-1),rwfl,DPT) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s %s
         '''%(M8R,rdat,rwfl),
              stdin=0,
              stdout=0)
    
def cdzom(imag,data,velo,rcoo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    obsolete('fdmod/cdzom','rtm/zofmigCD')

    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom

    rdat = imag+'rdat'
    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sreverse < ${SOURCES[0]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdat,DPT) +
         '''
         %sawefd2d < %s cden=y %s verb=n
         vel=${SOURCES[2]}
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdat,awepar+' jsnap=%d'%(par['nt']-1),rwfl,DPT) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s %s
         '''%(M8R,rdat,rwfl),
              stdin=0,
              stdout=0)
    
def cdzom3d(imag,data,velo,rcoo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    obsolete('fdmod/cdzom','rtm/zofmigCD')

    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom

    rdat = imag+'rdat'
    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sreverse < ${SOURCES[0]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdat,DPT) +
         '''
         %sawefd3d < %s cden=y %s verb=n
         vel=${SOURCES[2]}
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdat,awepar+' jsnap=%d'%(par['nt']-1),rwfl,DPT) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s %s
         '''%(M8R,rdat,rwfl),
              stdin=0,
              stdout=0)
    
# ------------------------------------------------------------
# wavefield-over-model plot
def wom(wom,wfld,velo,vmean,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    if(not par.has_key('wweight')): par['wweight']=1

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
        '''%(M8R,vmean,
             (par['nt']-1)/par['jsnap']+1,
             par['ot'],
             par['dt']*par['jsnap'],
             vtmp,DPT) 
	+
        '''
        %sadd scale=1,%g <%s %s >${TARGETS[0]};
        '''%(M8R,par['wweight'],vtmp,wtmp) 
	+
        '''
        %srm %s %s
        '''%(M8R,wtmp,vtmp),
        stdin=0,
        stdout=0)

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
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    if(not par.has_key('iweight')): par['iweight']=1

    itmp = imag + 'tmp'
    vtmp = imag + 'vel'

    Flow(iom,[velo,imag],
        '''
        %sscale < ${SOURCES[1]} axis=123 >%s datapath=%s/;
        '''%(M8R,itmp,DPT) 
	+
        '''
        %sadd < ${SOURCES[0]} add=-%g |
        scale axis=123
        >%s datapath=%s/;
        '''%(M8R,vmean,
             vtmp,DPT) 
	+
        '''
        %sadd scale=1,%g <%s %s >${TARGETS[0]};
        '''%(M8R,par['iweight'],vtmp,itmp) 
	+
        '''
        %srm %s %s
        '''%(M8R,itmp,vtmp),
        stdin=0,
        stdout=0)
 
# ------------------------------------------------------------
# wavefield snapshot plots
def wframe(frame,movie,index,custom,par):

    Flow([movie+'_plt',movie+'_bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s'%custom)
    Plot  (frame,[movie+'_plt',movie+'_bar'],
           'window n3=1 f3=%d |'%index
           + cgrey(custom,par))
    
# ------------------------------------------------------------
# elastic wavefield movie frames
def eframe(frame,movie,index,custom,axis,par,xscale=0.75,yscale=0.75,shift=-8.25):

    Flow([movie+'-plt',movie+'-bar'],movie,
         'byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)

    for i in range(2):
        Plot(frame+'-'+str(i),movie+'-plt',
             'window n3=1 f3=%d n4=1 f4=%d |' % (i,index)
             + cgrey('',par))
#        Result(frame+'-'+str(i),movie+'-plt',
#             'window n3=1 f3=%d n4=1 f4=%d |' % (i,index)
#             + cgrey('',par))

    if(axis==1): pplot.p2x1(frame,frame+'-1',frame+'-0',yscale,xscale,shift)
    else:        pplot.p1x2(frame,frame+'-0',frame+'-1',yscale,xscale,shift)

# ------------------------------------------------------------
# elastic wavefield movie
def emovie(movie,wfld,nframes,custom,axis,par,xscale=0.75,yscale=0.75,shift=-8.25):
    
    for iframe in range(nframes):
        tag = "-%02d" % iframe
        eframe(movie+tag,wfld,iframe,custom,axis,par,xscale,yscale,shift)
        
    allframes = map(lambda x: movie+'-%02d'  % x,range(nframes))
    Result(movie,allframes,'Movie')

# ------------------------------------------------------------
# elastic data
def edata(plot,data,custom,par):

    Flow([plot+'_plt',plot+'_bar'],data,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)
    
    for i in range(2):
	ctag = "-%d"%i
	
        Plot(  plot+ctag,[plot+'_plt',plot+'_bar'],
               'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
               + dgrey('screenratio=1.75 screenht=20 %s' % custom,par))
    pplot.p1x2(plot,plot+'-0',plot+'-1',0.5,0.5,-11)

    #    Result(plot+str(i+1),[plot+'_plt',plot+'_bar'],
    #           'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
    #           + dgrey('%s' % custom,par))

# ------------------------------------------------------------
# elastic image
def eimage(plot,imag,custom,par):

    Flow([plot+'_plt',plot+'_bar'],imag,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s ' % custom)        

    for i in range(4):        
        Plot  (plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
               + cgrey('%s'% custom,par))
        Result(plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
               + cgrey('%s'% custom,par))
        
# ------------------------------------------------------------
# plot elastic wavelet
def ewavelet(wavelet,custom,par):
    
    for i in range(2):
        Plot(wavelet+'-'+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             waveplot('%d %s'% (i,custom) ,par))
    pplot.p1x2(wavelet,wavelet+'-1',wavelet+'-2',0.5,0.5,-11.5)
def ewavelet3d(wavelet,custom,par):
    
    for i in range(3):
        Plot(wavelet+'-'+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             waveplot('%d %s'% (i,custom) ,par))
    pplot.p1x3(wavelet,wavelet+'-1',wavelet+'-2',wavelet+'-3',0.35,0.35,-11)
    
# ------------------------------------------------------------
# acoustic RTM
def artm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    if(not par.has_key('nbuf')): par['nbuf']=100
    
    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)
    
    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    
    # source wavefield (z,x,t)
    awefd(sout,swfl,sdat,velo,dens,sacq,iacq,iwindow+custom,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    awefd(rout,twfl,tdat,velo,dens,racq,iacq,iwindow+custom,par)
    Flow(rwfl,twfl,'reverse which=4 opt=i verb=y')
    
    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[swfl,rwfl],
         'xcor2d uu=${SOURCES[1]} axis=3 verb=y nbuf=%(nbuf)d ompnth=%(ompnth)d' % par)

# ------------------------------------------------------------
# elastic RTM
def ertm(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):

    if(not par.has_key('nbuf')): par['nbuf']=100
    
    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)
    
    iwindow = ' ' + \
              '''
              nqz=%(nqz)d oqz=%(oqz)g
              nqx=%(nqx)d oqx=%(oqx)g
              jsnap=%(jdata)d jdata=%(jdata)d
              ''' % par + ' '
    
    # source wavefield (z,x,t)
    ewefd2(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y opot=n" + custom+iwindow,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    ewefd2(rout,twfl,tdat,cccc,dens,racq,iacq," ssou=n opot=n" + custom+iwindow,par)

    # ------------------------------------------------------------
    Flow(swfl+'1',swfl,'window n3=1 f3=0')
    Flow(swfl+'2',swfl,'window n3=1 f3=1')

    Flow(rwfl+'1',twfl,'window n3=1 f3=0 | reverse which=4 opt=i verb=y')
    Flow(rwfl+'2',twfl,'window n3=1 f3=1 | reverse which=4 opt=i verb=y')

    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(imag+i+j,[swfl+i,rwfl+j],
                 'xcor uu=${SOURCES[1]} axis=3 verb=y nbuf=%(nbuf)d ompnth=%(ompnth)d' % par)
            
    Flow(imag,[imag+'11',imag+'12',imag+'21',imag+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def awefd2d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd2d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d
         nb=%(nb)d dabc=%(dabc)s
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)
    
def awefd3d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd3d
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d jdata=%(jdata)d
         nb=%(nb)d dabc=%(dabc)s
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd2d(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd2d
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nbell=%(nbell)d
         nb=%(nb)d  dabc=%(dabc)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

def ewefd3d(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         ewefd3d
         ompchunk=%(ompchunk)d  ompnth=%(ompnth)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nbell=%(nbell)d
         nb=%(nb)d  dabc=%(dabc)s
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
def gauss1x(gaus,xcen,xsig,par):
    Flow(gaus,None,
         '''
         math output="exp(-((x1-%g)*(x1-%g))/(2*%g))"
         ''' % (xcen,xcen,xsig*xsig) +
         '''
         n1=%(nx)d d1=%(dx)g o1=%(ox)g |
         scale axis=123
         ''' % par)
def gauss1z(gaus,zcen,zsig,par):
    Flow(gaus,None,
         '''
         math output="exp(-((x1-%g)*(x1-%g))/(2*%g))"
         ''' % (zcen,zcen,zsig*zsig) +
         '''
         n1=%(nz)d d1=%(dz)g o1=%(oz)g |
         scale axis=123
         ''' % par)
def gauss2d(gaus,xcen,zcen,xsig,zsig,par):
    Flow(gaus,None,
         '''
         math output="exp(-((x1-%g)*(x1-%g))/(2*%g)-((x2-%g)*(x2-%g))/(2*%g))"
	 ''' % (zcen,zcen,zsig*zsig,xcen,xcen,xsig*xsig) +
         '''
         n1=%(nz)d d1=%(dz)g o1=%(oz)g
         n2=%(nx)d d2=%(dx)g o2=%(ox)g |
         scale axis=123 
         ''' % par)

def gauss3d(gaus,xcen,ycen,zcen,xsig,ysig,zsig,par):
    Flow(gaus,None,
         '''
         math output="exp(-((x1-%g)*(x1-%g))/(2*%g)-((x2-%g)*(x2-%g))/(2*%g)-((x3-%g)*(x3-%g))/(2*%g))"
         ''' % (zcen,zcen,zsig*zsig,xcen,xcen,xsig*xsig,ycen,ycen,ysig*ysig) +
         '''
         n1=%(nz)d d1=%(dz)g o1=%(oz)g
         n2=%(nx)d d2=%(dx)g o2=%(ox)g 
	 n3=%(ny)d d3=%(dy)g o3=%(oy)g |
         scale axis=123
         ''' % par)

# ------------------------------------------------------------
def rbox2d(rbox,xlow,xhig,zlow,zhig,par):

    kx=int((xlow-par['ox'])/par['dx'])
    lx=int((xhig-par['ox'])/par['dx'])
    kz=int((zlow-par['oz'])/par['dz'])
    lz=int((zhig-par['oz'])/par['dz'])
    
    Flow(rbox,None,
         '''
         spike nsp=1 mag=1
         n1=%(nz)d d1=%(dz)g o1=%(oz)g
         n2=%(nx)d d2=%(dx)g o2=%(ox)g
         ''' % par +
         '''
         k1=%d l1=%d k2=%d l2=%d |
         scale axis=123
         '''%(kz,lz,kx,lx))


def quiver(vect,custom,par):

	Plot(vect+'o',vect,'window n1=1|' +
     	cgraph('squeeze=n plotcol=0 plotfat=3 symbol=. '+custom,par))
	
	Plot(vect+'h',vect,
	cgraph('squeeze=n plotcol=0 plotfat=1 '+custom,par))	

	Plot(vect,[vect+'h',vect+'o'],'Overlay')

# ------------------------------------------------------------
def hic2d(hic,wfs,wfr,cc,custom,par):
    Flow(hic,[wfs,wfr,cc],
         '''
         laps2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         ''' %par + ' %s '%custom)

def hic3d(hic,wfs,wfr,cc,custom,par):
    Flow(hic,[wfs,wfr,cc],
         '''
         laps3d verb=y
         nhx=%(nhx)d nhy=%(nhy)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         ''' %par + ' %s '%custom)
