from rsfproj import *
import spmig, sgmig, zomig,fdmod

def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if(par.has_key('verb')):
        p = p + ' verb='  +     par['verb']
    if(par.has_key('nrmax')):
        p = p + ' nrmax=' + str(par['nrmax'])
    if(par.has_key('dtmax')):
        p = p + ' dtmax=' + str(par['dtmax'])
    if(par.has_key('eps')):
        p = p + ' eps='   + str(par['eps'])
    if(par.has_key('tmx')):
        p = p + ' tmx='   + str(par['tmx'])
    if(par.has_key('tmy')):
        p = p + ' tmy='   + str(par['tmy'])
    if(par.has_key('pmx')):
        p = p + ' pmx='   + str(par['pmx'])
    if(par.has_key('pmy')):
        p = p + ' pmy='   + str(par['pmy'])
    if(par.has_key('misc')):
        p = p + ' '       +     par['misc']
    p = p + ' '
    return p

# ------------------------------------------------------------
def wempar(par):
    if(not par.has_key('verb')):    par['verb']='y'
    if(not par.has_key('eps')):     par['eps']=0.1

    if(not par.has_key('nrmax')):   par['nrmax']=1
    if(not par.has_key('dtmax')):   par['dtmax']=0.00005

    if(not par.has_key('tmx')):     par['tmx']=16
    if(not par.has_key('tmy')):     par['tmy']=16
    
    if(not par.has_key('incore')):  par['incore']='y'

# ------------------------------------------------------------
def slowness(slow,velo,par):
    Flow(slow,velo,
         '''
         window |
         math "output=1/input" |
         transp |
         spray axis=2 n=1 o=0 d=1 |
         put label2=y unit2=""
         ''')

# ------------------------------------------------------------
# WEM
# ------------------------------------------------------------

# WR: forward in time
def fWRwem(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wex causal=y %s slo=${SOURCES[1]} |
         window |
         transp
         ''' % param(par))

# WR: backward in time
def bWRwem(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wex causal=n %s slo=${SOURCES[1]} |
         window |
         transp
         ''' % param(par))

# ------------------------------------------------------------
# RTM
# ------------------------------------------------------------

# WR: forward in time
def fWRrtm(data,wfld,velo,dens,coor,custom,par):
    iwindow = ' ' + \
        '''
        nqz=%(nqz)d oqz=%(oqz)g
        nqx=%(nqx)d oqx=%(oqx)g
        jsnap=%(jdata)d jdata=%(jdata)d
        ''' % par + ' '

    fdmod.anifd2d(wfld+'_out',wfld,
                  data,velo,dens,coor,coor,iwindow+custom,par)

# WR: backward in time
def bWRrtm(data,wfld,velo,dens,coor,custom,par):
    iwindow = ' ' + \
        '''
        nqz=%(nqz)d oqz=%(oqz)g
        nqx=%(nqx)d oqx=%(oqx)g
        jsnap=%(jdata)d jdata=%(jdata)d
        ''' % par + ' '

    Flow(data+'_rev',data,'reverse which=2 opt=i verb=y')
    fdmod.anifd2d(wfld+'_out',wfld+'_rev',
                  data+'_rev',velo,dens,coor,coor,iwindow+custom,par)
    Flow(wfld,wfld+'_rev','reverse which=4 opt=i verb=y')

# ------------------------------------------------------------
# IC
# ------------------------------------------------------------

# CIC
def cic(imag,swfl,rwfl,custom,par):
    par['ciccustom'] = custom

    Flow(imag,[swfl,rwfl],
         '''
         xcor2d 
         uu=${SOURCES[1]} axis=3 verb=y ompnth=%(ompnth)d 
         %(ciccustom)s
         ''' % par)

# EIC
def eic(cip,swfl,rwfl,cc,custom,par):
    par['eiccustom'] = custom
    
    Flow(cip,[swfl,rwfl,cc],
         '''
         laps2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         %(eiccustom)s
         ''' %par)

# ------------------------------------------------------------
def wem(imag,sdat,rdat,slow,custom,par):
    
    fWRwem(sdat,swfl,slow,par)
    bWRwem(rdat,rwfl,slow,par)

    cic(imag,swfl,rwfl,custom,par)

# ------------------------------------------------------------
def rtm(imag,sdat,rdat,slow,custom,par):
    
    fWRrtm(sdat,swfl,slow,par)
    bWRrtm(rdat,rwfl,slow,par)

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
