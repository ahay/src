try:    from rsf.cluster import *
except: from rsf.proj    import *

import fdmod

def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if 'verb' in par:  p = p + ' verb='  +     par['verb']
    if 'nrmax' in par: p = p + ' nrmax=' + str(par['nrmax'])
    if 'dtmax' in par: p = p + ' dtmax=' + str(par['dtmax'])
    if 'eps' in par:   p = p + ' eps='   + str(par['eps'])
    if 'tmx' in par:   p = p + ' tmx='   + str(par['tmx'])
    if 'tmy' in par:   p = p + ' tmy='   + str(par['tmy'])
    if 'pmx' in par:   p = p + ' pmx='   + str(par['pmx'])
    if 'pmy' in par:   p = p + ' pmy='   + str(par['pmy'])
    if 'misc' in par:  p = p + ' '       +     par['misc']
    p = p + ' '
    return p

# ------------------------------------------------------------
def wempar(par):
    if 'verb' not in par:    par['verb']='y'
    if 'eps' not in par:     par['eps']=0.1

    if 'nrmax' not in par:   par['nrmax']=1
    if 'dtmax' not in par:   par['dtmax']=0.00005

    if 'tmx' not in par:     par['tmx']=16
    if 'tmy' not in par:     par['tmy']=16

# ------------------------------------------------------------
def eicpar(par):
    p = ' '
    if 'nhx' in par:
        p = p + ' nhx='   + str(par['nhx'])
    if 'nhy' in par:
        p = p + ' nhy='   + str(par['nhy'])
    if 'nhz' in par:
        p = p + ' nhz='   + str(par['nhz'])
    if 'nht' in par:
        p = p + ' nht='   + str(par['nht'])
    if 'oht' in par:
        p = p + ' oht='   + str(par['oht'])
    if 'dht' in par:
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


