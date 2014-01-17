try:    from rsf.cluster import *
except: from rsf.proj    import *

import fdmod

def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if(par.has_key('verb')):  p = p + ' verb='  +     par['verb']
    if(par.has_key('nrmax')): p = p + ' nrmax=' + str(par['nrmax'])
    if(par.has_key('dtmax')): p = p + ' dtmax=' + str(par['dtmax'])
    if(par.has_key('eps')):   p = p + ' eps='   + str(par['eps'])
    if(par.has_key('tmx')):   p = p + ' tmx='   + str(par['tmx'])
    if(par.has_key('tmy')):   p = p + ' tmy='   + str(par['tmy'])
    if(par.has_key('pmx')):   p = p + ' pmx='   + str(par['pmx'])
    if(par.has_key('pmy')):   p = p + ' pmy='   + str(par['pmy'])
    if(par.has_key('misc')):  p = p + ' '       +     par['misc']
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

# ------------------------------------------------------------
def eicpar(par):
    p = ' '
    if(par.has_key('nhx')):
        p = p + ' nhx='   + str(par['nhx'])
    if(par.has_key('nhy')):
        p = p + ' nhy='   + str(par['nhy'])
    if(par.has_key('nhz')):
        p = p + ' nhz='   + str(par['nhz'])
    if(par.has_key('nht')):
        p = p + ' nht='   + str(par['nht'])
    if(par.has_key('oht')):
        p = p + ' oht='   + str(par['oht'])
    if(par.has_key('dht')):
        p = p + ' dht='   + str(par['dht'])
    p = p + ' '
    return(p)

# ------------------------------------------------------------
# prepare slowness
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
# Wavefield Reconstruction
def wfl(wfld,data,slow,custom,par):
    Flow(wfld,[data,slow],
         '''
         wei verb=y irun=wfl causal=n
         slo=${SOURCES[1]}
         %s
         ''' % (param(par)+custom))

# ------------------------------------------------------------
# Wavefield Datuming
def dtm(ddat,data,slow,custom,par):
    Flow(ddat,[data,slow],
         '''
         wei verb=y irun=dtm causal=n
         slo=${SOURCES[1]}
         %s
         ''' % (param(par)+custom))

# ------------------------------------------------------------
# migration with Conventional Imaging Condition
def cicmig(icic,sdat,rdat,slow,custom,par):
    Flow(icic,[sdat,rdat,slow],
         '''
         wei verb=y irun=cic
         dat=${SOURCES[1]}
         slo=${SOURCES[2]}
         %s
         ''' % (param(par)+custom))

# ------------------------------------------------------------
# migration with Extended Imaging Condition
def eicmig(icic,ieic,sdat,rdat,slow,ccoo,custom,par):
    Flow([icic,ieic],[sdat,rdat,slow,ccoo],
         '''
         wei verb=y irun=eic
         dat=${SOURCES[1]}
         slo=${SOURCES[2]}
         coo=${SOURCES[3]}
         cip=${TARGETS[1]}
         %s
         ''' % (param(par)+eicpar(par)+custom))

# ------------------------------------------------------------
# migration with Extended Imaging Condition
def hicmig(icic,ieic,sdat,rdat,slow,ccoo,custom,par):
    Flow([icic,ieic],[sdat,rdat,slow,ccoo],
         '''
         wei verb=y irun=hic
         dat=${SOURCES[1]}
         slo=${SOURCES[2]}
         coo=${SOURCES[3]}
         cip=${TARGETS[1]}
         %s
         ''' % (param(par)+eicpar(par)+custom))
    
# ------------------------------------------------------------
# wavefields from arbitrary sources
def genwfl(wfl,sou,coo,slo,down,causal,custom,par):
    Flow(wfl,[sou,slo,coo],
         '''
         weigwf verb=y
         slo=${SOURCES[1]}
         coo=${SOURCES[2]}
         down=%s causal=%s 
         %s
         ''' %(down,causal,param(par)+custom))

# ------------------------------------------------------------
# 3D imaging condition
#def cic3d(imag,swfl,rwfl,par):
#    Flow(imag,[swfl,rwfl],
#         '''
#         cic uu=${SOURCES[1]} axis=3 verb=y ompnth=%(ompnth)d
#         ''' % par)
    
#def eic3d(cip,swfl,rwfl,cc,custom,par):
#    par['eiccustom'] = custom
#   
#    Flow(cip,[swfl,rwfl,cc],
#         '''
#         eic verb=y
#         nhx=%(nhx)d nhy=%(nhy)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
#         ur=${SOURCES[1]}
#         cc=${SOURCES[2]}
#         %(eiccustom)s
#         ''' %par)

# ------------------------------------------------------------
# FWI kernel (Z extrapolation)
def fwikerZ(ker,dws,ss,dwr,rr,slo,pad,custom,par):

     padL=min(par['nx'],int(0.5*(pad-par['nx'])))
     Flow(ker+'_padX',slo,
          '''
          remap1 n1=%d o1=%g d1=%g
          '''%(pad,par['ox']-padL*par['dx'],par['dx']))
          
     genwfl(ker+'_SW',dws,ss,ker+'_padX','y','y','',par)
     genwfl(ker+'_RW',dwr,rr,ker+'_padX','y','n','',par)

     Flow(ker,[ker+'_SW',ker+'_RW'],
          '''
          math output="conj(us)*ur"
          us=${SOURCES[0]} ur=${SOURCES[1]} |
          window min1=%g n1=%g |
          stack axis=3 | real | transp
          '''%(par['ox'],par['nx']),stdin=0)

# ------------------------------------------------------------
# FWI kernel (X extrapolation)
def fwikerX(ker,dws,ss,dwr,rr,slo,pad,custom,par):
     Flow(ss+'_T',ss,'reverse which=1 opt=i')
     Flow(rr+'_T',rr,'reverse which=1 opt=i')

     padT=min(par['nz'],int(0.5*(pad-par['nz'])))
     Flow(ker+'_padZ',slo,
          '''
          window | transp | 
          remap1 n1=%d o1=%g d1=%g |
          transp plane=23 
          '''%(pad,par['oz']-padT*par['dz'],par['dz']))
     
     genwfl(ker+'_SW',dws,ss+'_T',ker+'_padZ','y','y','',par)
     genwfl(ker+'_RW',dwr,rr+'_T',ker+'_padZ','n','n','',par)

     Flow(ker,[ker+'_SW',ker+'_RW'],
          '''
          math output="conj(us)*ur"
          us=${SOURCES[0]} ur=${SOURCES[1]} |
          window min1=%g n1=%g |
          stack axis=3 | real 
          '''%(par['oz'],par['nz']),stdin=0)


