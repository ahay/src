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
# WEXIMG
# ------------------------------------------------------------

# Wavefield Reconstruction: forward in time
def fWRwex(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl causal=y %s slo=${SOURCES[1]} 
         ''' % param(par))

# Wavefield Reconstruction: backward in time
def bWRwex(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl causal=n %s slo=${SOURCES[1]} 
         ''' % param(par))

# Causal datuming (forward in time, causal=y)
# (useful for datuming source wavefields)

def Cdtwex(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl inv=0 datum=1 causal=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# Anti-causal datuming (backward in time, causal=n)
# (useful for datuming receiver wavefields)

def Adtwex(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl inv=0 datum=1 causal=n %s
         slo=${SOURCES[1]}
         ''' % param(par))

# Modeling
def wexMOD(ref,data,slow,wfls,par)
    Flow(data,[ref,slow,wfls],
         '''
         rtoc |
         sfweximg
         adj=0 save=0 feic=0 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         ''' %param(par))

# Migration
def wexCIMG(img,data,slow,wfls,par)
    Flow(img,[data,slow,wfls],
         '''
         weximg
         adj=1 save=0 feic=0 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]} | 
         real
         ''' % param(par))

def wexEIMG(cip,data,slow,wfls,cc,par)
    Flow(cip,[data,slow,wfls,cc],
         '''
         weximg %s
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         adj=1 save=0 feic=1 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         cc=${SOURCES[3]}
         ''' % param(par),par)

# ------------------------------------------------------------
# WEXMVA
# ------------------------------------------------------------

# forward (C.I.C.)
def wexCFMVA(dimg,dslo,slow,wfls,wflr,par)
    Flow(dimg,[dslo,wfls,wflr,slow],
         '''
         wexmva
         adj=0 feic=0 verb=y %s
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))

# adjoint (C.I.C.)
def wexCAMVA(dslo,dimg,slow,wfls,wflr,par)
    Flow(dslo,[dimg,wfls,wflr,slow],
         '''
         wexmva
         adj=1 feic=0 verb=y %s
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))

# forward (E.I.C.)
def wexEFMVA(dimg,dslo,slow,wfls,wflr,cc,par)
    Flow(dimg,[dslo,wfls,wflr,slow,cc],
         '''
         wexmva %s
         adj=0 feic=1 verb=y
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % param(par),par)

# adjoint (E.I.C.)
def wexEFMVA(dslo,dimg,slow,wfls,wflr,cc,par)
    Flow(dslo,[dimg,wfls,wflr,slow,cc],
         '''
         wexmva %s
         adj=1 feic=1 verb=y
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % param(par),par)



