try:    from rsf.cluster import *
except: from rsf.proj    import *
from rsf.recipes import spmig,sgmig,zomig,fdmod
import random
import functools, operator

random.seed(1004)
def add(x,y): return x+y
def myid(n): return '_'+functools.reduce(operator.add,['%d'%random.randint(0,9) for i in range(n)])

def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if('verb' in par):
        p = p + ' verb='  +     par['verb']
    if('nrmax' in par):
        p = p + ' nrmax=' + str(par['nrmax'])
    if('dtmax' in par):
        p = p + ' dtmax=' + str(par['dtmax'])
    if('eps' in par):
        p = p + ' eps='   + str(par['eps'])
    if('tmx' in par):
        p = p + ' tmx='   + str(par['tmx'])
    if('tmy' in par):
        p = p + ' tmy='   + str(par['tmy'])
    if('pmx' in par):
        p = p + ' pmx='   + str(par['pmx'])
    if('pmy' in par):
        p = p + ' pmy='   + str(par['pmy'])
    if('misc' in par):
        p = p + ' '       +     par['misc']
    p = p + ' '
    return p

# ------------------------------------------------------------
def wempar(par):
    if('verb' not in par):    par['verb']='y'
    if('eps' not in par):     par['eps']=0.1

    if('nrmax' not in par):   par['nrmax']=1
    if('dtmax' not in par):   par['dtmax']=0.00005

    if('tmx' not in par):     par['tmx']=16
    if('tmy' not in par):     par['tmy']=16
    
    if('incore' not in par):  par['incore']='y'

# ------------------------------------------------------------
def slowness2d(slow,velo,par):
    Flow(slow,velo,
         '''
         math "output=1/input" |
         transp plane=34 | 
         transp plane=12 | 	
	 transp plane=23 |
         put label2=y o2=0 d2=1
         ''')

def slowness(slow,velo,par):
    slowness2d(slow,velo,par)
    
# ------------------------------------------------------------
# WEM
# ------------------------------------------------------------

# WR: forward in time
def fWRwem(data,wfld,slow,par):
#    Flow(wfld,[data,slow],
#         '''
#         wex causal=y %s slo=${SOURCES[1]} |
#         window |
#         transp
#         ''' % param(par))
    Flow(wfld+'_tmp',[data,slow],
         'wexwfl causal=y %s slo=${SOURCES[1]}' % param(par))
    Flow(wfld,wfld+'_tmp','window | transp')

# WR: backward in time
def bWRwem(data,wfld,slow,par):
#    Flow(wfld,[data,slow],
#         '''
#         wex causal=n %s slo=${SOURCES[1]} |
#         window |
#         transp
#         ''' % param(par))
    Flow(wfld+'_tmp',[data,slow],
         'wexwfl causal=n %s slo=${SOURCES[1]}' % param(par))
    Flow(wfld,wfld+'_tmp','window | transp')


def wemWR(data,wfld,slow,causal,par):
    Flow(wfld,[data,slow],
         'wexwfl %s slo=${SOURCES[1]}' % param(par) + ' causal=%s '%causal) 

def iwindow(par):
    win = ' ' + \
          '''
          nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g 
          nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g
          jsnap=%(jdata)d jdata=%(jdata)d
          ''' % par + ' '

    return win
    
# ------------------------------------------------------------
# RTM
# ------------------------------------------------------------

# WR: forward in time
def fWRrtm(data,wfld,velo,dens,coor,custom,par):
    fdmod.anifd2d(wfld+'_',wfld,
                  data,velo,dens,coor,coor,iwindow(par)+custom,par)
# WR: backward in time
def bWRrtm(data,wfld,velo,dens,coor,custom,par):
    Flow(data+'_R',data,'reverse which=2 opt=i verb=y')
    fdmod.anifd2d(wfld+'_',wfld+'_R',
                  data+'_R',velo,dens,coor,coor,iwindow(par)+custom,par)

# WR: forward in time
def fWRawe(data,wfld,velo,dens,coor,custom,par):
    fdmod.awefd2d(wfld+'_',wfld,
                  data,velo,dens,coor,coor,iwindow(par)+custom,par)
# WR: backward in time
def bWRawe(data,wfld,velo,dens,coor,custom,par):
    Flow(data+'_R',data,'reverse which=2 opt=i verb=y')
    fdmod.awefd2d(wfld+'_',wfld+'_R',
                  data+'_R',velo,dens,coor,coor,iwindow(par)+custom,par)
    Flow(wfld,wfld+'_R','reverse which=2 opt=i verb=y')

# WR: forward in time
def fWRcda(data,wfld,velo,coor,custom,par):
    fdmod.cdafd(wfld+'_',wfld,
                data,velo,coor,coor,iwindow(par)+custom,par)
# WR: backward in time
def bWRcda(data,wfld,velo,coor,custom,par):
    Flow(data+'_R',data,'reverse which=2 opt=i verb=y')
    fdmod.cdafd(wfld+'_',wfld,
                data+'_R',velo,coor,coor,iwindow(par)+custom,par)

def reverse(wfldO,wfldI,par):
    Flow(wfldO,wfldI,'reverse which=4 opt=i verb=y')
    
# ------------------------------------------------------------
# variable-density shot-record migration
def rtmcic(imag,velo,dens,
           sdat,scoo,
           rdat,rcoo,
           custom,par):

    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom

    swfl=imag+'swfl'+myid(16)
    rdrv=imag+'rdrv'+myid(16)
    rwfl=imag+'rwfl'+myid(16)

    Flow(imag,[sdat,scoo,rdat,rcoo,velo],
         '''
         %sawefd2d < ${SOURCES[0]} cden=y %s verb=n
         vel=${SOURCES[4]}
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+' '+awepar,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=y %s verb=n
         vel=${SOURCES[4]}
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+' '+awepar,rwfl,DPT) +
         '''
         %scicold2d <%s isreversed=0 ur=%s axis=3 verb=n %s >${TARGETS[0]};
         '''%(M8R,swfl,rwfl,custom) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)

# constant-density shot-record migration
def cdrtm(imag,velo,
          sdat,scoo,
          rdat,rcoo,
          custom,par):

    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom

    swfl=imag+'swfl'+myid(16)
    rdrv=imag+'rdrv'+myid(16)
    rwfl=imag+'rwfl'+myid(16)

    Flow(imag,[sdat,scoo,rdat,rcoo,velo],
         '''
         %sawefd2d < ${SOURCES[0]} cden=y %s verb=n
         vel=${SOURCES[4]}
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+' '+awepar,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=y %s verb=n
         vel=${SOURCES[4]}
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+' '+awepar,rwfl,DPT) +
         '''
         %scicold2d <%s isreversed=0 ur=%s axis=3 verb=n %s >${TARGETS[0]};
         '''%(M8R,swfl,rwfl,custom) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)

#def cdrtmOBSOLETE(imag,velo,
#          sdat,scoo,
#          rdat,rcoo,
#          custom,par):
#
#    M8R='$RSFROOT/bin/sf'
#    DPT=os.environ.get('TMPDATAPATH')
#
#    awewin = 'nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g'%par
#    awepar = 'ompchunk=%(ompchunk)d ompnth=%(ompnth)d verb=y free=n snap=%(snap)s jsnap=%(jdata)d jdata=%(jdata)d dabc=%(dabc)s nb=%(nb)d'%par + ' ' + custom
#
#    swfl=imag+'swfl'
#    rdrv=imag+'rdrv'
#    rwrv=imag+'rwrv'
#    rwfl=imag+'rwfl'
#
#    Flow(imag,[sdat,scoo,rdat,rcoo,velo],
#         '''
#         %sawefd2d < ${SOURCES[0]} cden=y %s
#         vel=${SOURCES[4]}
#         sou=${SOURCES[1]}
#         rec=${SOURCES[1]}
#         wfl=%s datapath=%s/
#         >/dev/null;
#         '''%(M8R,awewin+' '+awepar,swfl,DPT) +
#         '''
#         %sreverse < ${SOURCES[2]} which=2 opt=i verb=y >%s datapath=%s/;
#         '''%(M8R,rdrv,DPT) +
#         '''
#         %sawefd2d < %s cden=y %s
#         vel=${SOURCES[4]}
#         sou=${SOURCES[3]}
#         rec=${SOURCES[3]}
#         wfl=%s datapath=%s/
#         >/dev/null;
#         '''%(M8R,rdrv,awewin+' '+awepar,rwrv,DPT) +
#         '''
#         %sreverse < %s which=4 opt=i verb=y >%s datapath=%s/;
#         '''%(M8R,rwrv,rwfl,DPT) +
#         '''
#         %sxcor2d <%s uu=%s axis=3 verb=y %s >${TARGETS[0]};
#         '''%(M8R,swfl,rwfl,custom) +
#         '''
#         %srm %s %s %s %s
#         '''%(M8R,swfl,rdrv,rwrv,rwfl),
#              stdin=0,
#              stdout=0)

# ------------------------------------------------------------
# IC
# ------------------------------------------------------------

# CIC: cross-correlation
def cic(imag,swfl,rwfl,custom,par,isreversed=0):
    Flow(imag,[swfl,rwfl],
         '''
         cicold2d verb=y
         ur=${SOURCES[1]} 
         '''%par+ 'isreversed=%d'%isreversed + ' ' + custom )
    
# EIC
def eic(cip,swfl,rwfl,cc,custom,par,isreversed=0):    
    Flow(cip,[swfl,rwfl,cc],
         '''
         eicold2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         '''%par+ 'isreversed=%d'%isreversed + ' '+ custom )
    
# CIC: deconvolution
def dic(imag,swfl,rwfl,eps,custom,par):
    par['diccustom'] = custom
    
    Flow(imag,[swfl,rwfl],
         '''
         math s=${SOURCES[0]} r=${SOURCES[1]} 
         output="-(conj(s)*r)/(conj(s)*s+%g)" |
         transp plane=23 | stack | real
         ''' %eps)

# ------------------------------------------------------------
def wem(imag,sdat,rdat,slow,custom,par):
    
    fWRwem(sdat,swfl,slow,par)
    bWRwem(rdat,rwfl,slow,par)

    cic(imag,swfl,rwfl,custom,par)

# ------------------------------------------------------------
def rtm(imag,sdat,rdat,velo,dens,custom,par):
    
    fWRrtm(sdat,swfl,velo,dens,par)
    bWRrtm(rdat,rwfl,velo,dens,par)

    cic(imag,swfl,rwfl,custom,par)

# ------------------------------------------------------------
# SURVEY-SINKING migration
# ------------------------------------------------------------
#def sinking(par):
#    # surface data
#    sgmig.wflds('d0','cmps',par)
#
#    # datuned data
#    sgmig.datum('d1','sd','d0',par)
#    
#    for k in ('0','1'):
#        s = 's' + k # slowness
#        d = 'd' + k # prestack data
#        i = 'i' + k # prestack image
#        e = 'e' + k # inner offset data
#        z = 'z' + k # inner offset image
#
#        # prestack migration
#        sgmig.image(i,s,d,par)
#        
#        # near offset migration
#        Flow(e,d,'window squeeze=n n3=8')
#        sgmig.image(z,s,e,par)


# ------------------------------------------------------------
# SHOT-PROFILE MIGRATION
# ------------------------------------------------------------
#def profile(par):
#    # surface wavefields
#    spmig.wflds('d0s','d0r','wave','shot',par)
#    
#    # datumed wavefields
#    spmig.datum('d1s','d1r','sd','d0s','d0r',par)
#    
#    for k in ('0','1'):
#        s = 's' + k       # slowness
#        j = 'j' + k       # image
#        ds= 'd' + k + 's' # source   wavefield
#        dr= 'd' + k + 'r' # receiver wavefield
#        
#        spmig.image(j,s,ds,dr,par)
#        
# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
#def igrey(custom,par):
#    return '''
#    grey labelrot=n title="" %s min2=%g max2=%g min1=%g max1=%g
#    ''' % (custom,par['xmin'],par['xmax'],par['zmin'],par['zmax'])

#def result(par):
#
#    for k in ('0','1','2','d','o'):
#        s = 's' + k
#        Result(s,s,'window      | transp |'
#               +igrey('title=s pclip=100 color=j allpos=y',par))
#
#        for k in ('0','1'):
#            z = 'z' + k
#            i = 'i' + k
#            j = 'j' + k
#            
#            Result(z,z,'window n3=1 | transp |'
#                   +igrey('title=z pclip=99',par))
#            Result(i,i,'window n3=1 | transp |'
#                   +igrey('title=i pclip=99',par))
#            Result(j,j,'window      | transp |'
#                   +igrey('title=j pclip=99',par))
#
# ------------------------------------------------------------
