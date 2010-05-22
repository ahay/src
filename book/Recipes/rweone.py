from rsf.proj import *

def ccgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y wanttitle=n
    title="" pclip=99 label1="z" unit1=m label2="x" unit2=m %s
    min1=%g max1=%g min2=%g max2=%g
    ''' % (custom,par['zmin'],par['zmax'],par['xmin'],par['xmax'])

def rcgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y wanttitle=n
    title="" pclip=99 label1="t" unit1=s label2="x" unit2=m %s
    ''' % (custom)

def ccgraph(custom,par):
    return '''
    graph labelrot=n %s
    yreverse=y wantaxis=n title=" " 
    min1=%g max1=%g min2=%g max2=%g
    ''' % (custom,par['xmin'],par['xmax'],par['zmin'],par['zmax'])

# plot coordinate system
def cos(cos,par):
    wft = cos + '-wft'
    ray = cos + '-ray'

    Plot(ray,cos,'window j1=%(jray)d | transp |' % par
         + ccgraph('plotcol=1',par))
    Plot(wft,cos,'window j2=%(jwft)d |' % par
         + ccgraph('plotcol=2',par))
    
    Plot(cos,[ray,wft],'Overlay')

# plot CC grid
def ccplot(dat,opt,par):
    Plot(dat,dat,'transp |' +ccgrey(opt,par))

# plot RC grid
def rcplot(dat,opt,par):
    Plot(dat,dat,'transp |' +rcgrey(opt,par))

#def cosplot(cos,vel,par):
#    wft = cos + '-wft'
#    ray = cos + '-ray'
#
#    Plot(vel,ccgrey('allpos=y bias=1000',par))
#
#    Plot(ray,cos,'window j1=%(jray)d | transp |' % par
#         + ccgraph('plotcol=1',par))
#    Plot(wft,cos,'window j2=%(jwft)d |' % par
#         + ccgraph('plotcol=2',par))
#    
#    Plot(cos,[ray,wft],'Overlay')
    
# compute RC coefficients
def abm(abm,abr,slo,cos,par):

    _slo = '_' + slo
    Flow(_slo,[slo,cos],
         '''
         c2r rays=${SOURCES[1]} adj=n |
         put label1=t label2=g
         ''')
    Flow([abm,abr],[cos,_slo],
         '''
         rweab naref=1 nbref=1
         slo=${SOURCES[1]} abr=${TARGETS[1]}
         ''')

# transform data to frequency
def freq(frq,dat,nw,fw,jw):
    Flow(frq,dat,
         '''
         fft1 inv=n opt=n |
         window squeeze=n min1=1 j1=%d f1=%d n1=%d |
         transp |
         put label1=g label2=t label3=w
         ''' % (jw,fw,nw) )

# plot 2x3 images in RC
def imrc(img,frq,abm,abr,method,par):

    if(method=='F15'): method='method=0 c1=0.50   c2=0.00'
    if(method=='F45'): method='method=0 c1=0.50   c2=0.25'
    if(method=='F60'): method='method=0 c1=0.4761 c2=0.3767'
    if(method=='SSF'): method='method=1'
    if(method=='FFD'): method='method=2 c1=0.50   c2=0.00'
    if(method=='PSC'): method='method=3 c1=0.50   c2=0.00'
    
    Flow(img,[frq,abm,abr],
         '''
         rwezomig ntap=100 %s
         abm=${SOURCES[1]}
         abr=${SOURCES[2]} |
         put label1=g label2=t
         ''' % method)

# plot 2x3 images in CC
def imcc(imgCC,imgRC,cos,method,par):
    
    Flow(imgCC,[imgRC,cos],
         '''
         c2r rays=${SOURCES[1]} adj=y linear=n 
         nsz=5 nsx=5
         a2n=%(nz)d a2o=%(oz)g a2d=%(dz)g
         a1n=%(nx)d a1o=%(ox)g a1d=%(dx)g |
         put label1=z label2=x
         ''' % par)

def complot(plot,p0,p1,p2,p3,p4,p5):

    j0 = '_' + p0
    j1 = '_' + p1
    j2 = '_' + p2
    j3 = '_' + p3
    j4 = '_' + p4
    j5 = '_' + p5
    
    Plot(j0,p0,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-10.0 xcenter=-0.0')
    Plot(j1,p1,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-10.0 xcenter=-12.0')
    Plot(j2,p2,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-10.0 xcenter=-24.0')
    Plot(j3,p3,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-00.0 xcenter=-0.0')
    Plot(j4,p4,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-00.0 xcenter=-12.0')
    Plot(j5,p5,'Overlay',vppen='yscale=.35 xscale=.35 ycenter=-00.0 xcenter=-24.0')

    Result(plot,[j0,j1,j2,j3,j4,j5],'Overlay')
