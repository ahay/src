try:    from rsf.cluster import *
except: from rsf.proj    import *

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

def eicpar(par):
    p = ' '
    if('nhx' in par):
        p = p + ' nhx='   + str(par['nhx'])
    if('nhy' in par):
        p = p + ' nhy='   + str(par['nhy'])
    if('nhz' in par):
        p = p + ' nhz='   + str(par['nhz'])
    if('nht' in par):
        p = p + ' nht='   + str(par['nht'])
    if('oht' in par):
        p = p + ' oht='   + str(par['oht'])
    if('dht' in par):
        p = p + ' dht='   + str(par['dht'])
    p = p + ' '
    return(p)


# ------------------------------------------------------------
# prepare slowness
def slowness(slow,velo,par):
    Flow(slow,velo,
         '''
         math "output=1/input" |
         transp plane=12 | 	
	 transp plane=23
         ''')

# ------------------------------------------------------------
# WEXIMG
# ------------------------------------------------------------

# Wavefield Reconstruction: forward in time
def fWRwex(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl verb=y causal=y %s slo=${SOURCES[1]} 
         ''' % param(par))

# Wavefield Reconstruction: backward in time
def bWRwex(data,wfld,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl verb=y causal=n %s slo=${SOURCES[1]} 
         ''' % param(par))

# Wavefield Reconstruction
def wfl(data,wfld,slow,custom,par):
    Flow(wfld,[data,slow],
         '''
         wex verb=y irun=wfl causal=n
         slo=${SOURCES[1]}
         %s
         ''' % (param(par)+custom))

# Wavefield Datuming
def dtm(ddtm,data,slow,custom,par):
    Flow(ddtm,[ddat,slow],
         '''
         wex verb=y irun=dtm causal=n
         slo=${SOURCES[1]}
         %s
         ''' % (param(par)+custom))

# migration with Conventional Imaging Condition
def cic(icic,sdat,rdat,slow,custom,par):
    Flow(icic,[sdat,rdat,slow],
         '''
         wex verb=y irun=cic
         dat=${SOURCES[1]}
         slo=${SOURCES[2]}
         %s
         ''' % (param(par)+custom))

# Causal datuming (forward in time, causal=y)
# (useful for datuming source wavefields)

def Cdtwex(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl verb=y inv=0 datum=1 causal=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# Anti-causal datuming (backward in time, causal=n)
# (useful for datuming receiver wavefields)

def Adtwex(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         wexwfl verb=y inv=0 datum=1 causal=n %s
         slo=${SOURCES[1]}
         ''' % param(par))


def datum(swf1,rwf1,slow,swf0,rwf0,par):
    Cdtwex(swf1,swf0,slow,par) #      causal
    Adtwex(rwf1,rwf0,slow,par) # anti-causal


# Modeling
def wexMOD(ref,data,slow,wfls,par):
    Flow(data,[ref,slow,wfls],
         '''
         rtoc |
         weximg %s
         adj=0 save=0 feic=0 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         ''' %param(par))

# Migration
def wexCIMG(img,data,slow,wfls,par):
    Flow(img,[data,slow,wfls],
         '''
         weximg %s
         adj=1 save=0 feic=0 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]} | 
         real
         ''' %param(par))

def wexEIMG(cip,data,slow,wfls,cc,par):
    par['temp'] = param(par)
    Flow(cip,[data,slow,wfls,cc],
         '''
         weximg %(temp)s
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         adj=1 save=0 feic=1 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         cc=${SOURCES[3]}
         ''' % par)
    par.pop('temp')

def wexXIMG(cig,data,slow,wfls,gg,par):
    par['temp'] = param(par)
    Flow(cig,[data,slow,wfls,gg],
         '''
         weximg %(temp)s
         nhx=%(nhx)d nhz=0 nhy=%(nhy)d nht=0 dht=%(dht)g
         adj=1 save=0 feic=1 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         cc=${SOURCES[3]}
         ''' % par)
    par.pop('temp')

def wexTIMG(cig,data,slow,wfls,gg,par):
    par['temp'] = param(par)
    Flow(cig,[data,slow,wfls,gg],
         '''
         weximg %(temp)s
         nhx=0 nhz=0 nhy=0 nht=%(nht)d dht=%(dht)g
         adj=1 save=0 feic=1 verb=y
         slo=${SOURCES[1]}
         swfl=${SOURCES[2]}
         cc=${SOURCES[3]}
         ''' % par)
    par.pop('temp')

def migcic(imag,sdat,rdat,slow,custom,par):
    Flow(imag,[rdat,sdat,slow],
         '''
         wexmig
         adj=1 save=0 feic=0 verb=y
         swfl=${SOURCES[1]}
         slo=${SOURCES[2]}
	 %s 
         ''' % (param(par)+custom))

def migeic(imag,cips,sdat,rdat,slow,cc,custom,par):
    Flow([imag,cips],[rdat,sdat,slow,cc],
         '''
         wexmig
         adj=1 save=0 feic=1 verb=y 
         swfl=${SOURCES[1]}
         slo=${SOURCES[2]}
	 cc=${SOURCES[3]}
         cip=${TARGETS[1]}
         %s
         ''' % (param(par)+eicpar(par)+custom))


# ------------------------------------------------------------
# WEXMVA
# ------------------------------------------------------------

# forward (C.I.C.)
def wexCFMVA(dimg,dslo,slow,wfls,wflr,par):
    Flow(dimg,[dslo,wfls,wflr,slow],
         '''
         wexmva
         adj=0 feic=0 verb=y %s
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))

# adjoint (C.I.C.)
def wexCAMVA(dslo,dimg,slow,wfls,wflr,par):
    Flow(dslo,[dimg,wfls,wflr,slow],
         '''
         wexmva
         adj=1 feic=0 verb=y %s
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))

# forward (E.I.C.)
def wexEFMVA(dimg,dslo,slow,wfls,wflr,cc,par):
    par['temp'] = param(par)

    Flow(dimg,[dslo,wfls,wflr,slow,cc],
         '''
         wexmva %(temp)s
         adj=0 feic=1 verb=y
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % par)
    par.pop('temp')

# adjoint (E.I.C.)
def wexEAMVA(dslo,dimg,slow,wfls,wflr,cc,par):
    par['temp'] = param(par)

    Flow(dslo,[dimg,wfls,wflr,slow,cc],
         '''
         wexmva %(temp)s
         adj=1 feic=1 verb=y
         nhx=%(nhx)d nhz=%(nhz)d nhy=%(nhy)d nht=%(nht)d dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % par)
    par.pop('temp')

# forward (X.I.C.)
def wexXFMVA(dimg,dslo,slow,wfls,wflr,cc,par):
    par['temp'] = param(par)

    Flow(dimg,[dslo,wfls,wflr,slow,cc],
         '''
         wexmva %(temp)s
         adj=0 feic=1 verb=y
         nhx=%(nhx)d nhz=0 nhy=%(nhy)d nht=0 dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % par)
    par.pop('temp')

# adjoint (E.I.C.)
def wexXAMVA(dslo,dimg,slow,wfls,wflr,cc,par):
    par['temp'] = param(par)

    Flow(dslo,[dimg,wfls,wflr,slow,cc],
         '''
         wexmva %(temp)s
         adj=1 feic=1 verb=y
         nhx=%(nhx)d nhz=0 nhy=%(nhy)d nht=0 dht=%(dht)g
         swfl=${SOURCES[1]}
         rwfl=${SOURCES[2]}
         slo=${SOURCES[3]}
         cc=${SOURCES[4]}
         ''' % par)
    par.pop('temp')



