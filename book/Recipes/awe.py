try:    from rsf.cluster import *
except: from rsf.proj    import *

# ------------------------------------------------------------
def param(par):
    if(not par.has_key('nb')):       par['nb']=0
    if(not par.has_key('nbell')):    par['nbell']=5

    if(not par.has_key('snap')):     par['snap']='y'
    if(not par.has_key('jsnap')):    par['jsnap']=100
    if(not par.has_key('jdata')):    par['jdata']=1
    if(not par.has_key('dabc')):     par['dabc']='y'
    if(not par.has_key('ompchunk')): par['ompchunk']=1
    if(not par.has_key('ompnth')):   par['ompnth']=0
    if(not par.has_key('fsrf')):     par['fsrf']='n'
    if(not par.has_key('verb')):     par['verb']='n'

    if(not par.has_key('gaus')):     par['gaus']='y'
    
# ------------------------------------------------------------
def awepar(par):
    awe = ' ' + \
          '''
          ompchunk=%(ompchunk)d ompnth=%(ompnth)d
          verb=%(verb)s fsrf=%(fsrf)s
          snap=%(snap)s jsnap=%(jsnap)d jdata=%(jdata)d
          dabc=%(dabc)s nb=%(nb)d
          '''%par + ' '
    return awe

def eicpar(par):
    eic = ' ' + \
          '''
          nhx=%(nhx)d nhy=%(nhy)d nhz=%(nhz)d nht=%(nht)d
          gaus=%(gaus)s
          '''%par + ' '
    return eic

# ------------------------------------------------------------
# wavelet
def wavelet(wav,frq,custom,par):    
    Flow(wav,None,
         '''
         spike nsp=1 mag=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g k1=%(kt)d |
         pad end1=%(nt)d |
         '''%par + 
         '''
         ricker1 frequency=%g |
         '''%frq + 
         '''      
         window n1=%(nt)d |
         scale axis=123 |
         put label1=%(lt)s unit1=%(ut)s |
         transp
         '''%par)   

# ------------------------------------------------------------
# variable-density acoustic FD modeling
def awefd2d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd2d cden=n
         vel=${SOURCES[1]} den=${SOURCES[2]}
         sou=${SOURCES[3]} rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         ''' + ' ' + awepar(par) + ' ' + custom)

def awefd3d(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awefd3d cden=n
         vel=${SOURCES[1]} den=${SOURCES[2]}
         sou=${SOURCES[3]} rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         ''' + ' ' + awepar(par) + ' ' + custom)
         
# ------------------------------------------------------------
# constant-density acoustic FD modeling
def cdafd2d(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd2d cden=y
         vel=${SOURCES[1]}
         sou=${SOURCES[2]} rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         ''' + ' ' + awepar(par) + ' ' + custom)

def cdafd3d(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd3d cden=y
         vel=${SOURCES[1]}
         sou=${SOURCES[2]} rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         ''' + ' ' + awepar(par) + ' ' + custom)
         
# ------------------------------------------------------------
# AWE modeling operator
def aweop2d(wfl,sou,vel,custom,par):
    Flow(wfl,[sou,vel],
         '''
         aweop2d
         vel=${SOURCES[1]}
         '''+ ' ' + awepar(par) + ' ' + custom)

def aweop3d(wfl,sou,vel,custom,par):
    Flow(wfl,[sou,vel],
         '''
         aweop3d
         vel=${SOURCES[1]}
         '''+ ' ' + awepar(par) + ' ' + custom)
         
# ------------------------------------------------------------
# CIC forward: opr * wfl = img
def cic2dF(img,wfl,opr,custom,par):
    Flow(img,[wfl,opr],
         '''
         cicop2d adj=n wflcausal=y oprcausal=n 
         opr=${SOURCES[1]} 
         ''' + ' ' + custom)
def cic3dF(img,wfl,opr,custom,par):
    Flow(img,[wfl,opr],
         '''
         cicop3d adj=n wflcausal=y oprcausal=n 
         opr=${SOURCES[1]} 
         ''' + ' ' + custom)

# CIC adjoint: opr * img = wfl
def cic2dA(wfl,img,opr,custom,par):
    Flow(wfl,[img,opr],
         '''
         cicop2d adj=y wflcausal=n oprcausal=n 
         opr=${SOURCES[1]} 
         ''' + ' ' + custom)
def cic3dA(wfl,img,opr,custom,par):
    Flow(wfl,[img,opr],
         '''
         cicop3d adj=y wflcausal=n oprcausal=n 
         opr=${SOURCES[1]} 
         ''' + ' ' + custom)
         
# ------------------------------------------------------------         
# EIC forward: opr * wfl = img
def eic2dF(img,wfl,opr,cip,custom,par):    
    Flow(img,[wfl,opr,cip],
         '''
         eicop2d adj=n wflcausal=y oprcausal=n 
         opr=${SOURCES[1]} cip=${SOURCES[2]}
         ''' + eicpar(par) + ' ' + custom )
def eic3dF(img,wfl,opr,cip,custom,par):    
    Flow(img,[wfl,opr,cip],
         '''
         eicop3d adj=n wflcausal=y oprcausal=n 
         opr=${SOURCES[1]} cip=${SOURCES[2]}
         ''' + eicpar(par) + ' ' + custom )
         
# EIC adjoint: opr * img = wfl
def eic2dA(wfl,img,opr,cip,custom,par):    
    Flow(wfl,[img,opr,cip],
         '''
         eicop2d adj=y wflcausal=n oprcausal=n 
         opr=${SOURCES[1]} cip=${SOURCES[2]}
         ''' + ' ' + custom )
def eic3dA(wfl,img,opr,cip,custom,par):    
    Flow(wfl,[img,opr,cip],
         '''
         eicop3d adj=y wflcausal=n oprcausal=n 
         opr=${SOURCES[1]} cip=${SOURCES[2]}
         ''' + ' ' + custom )
         
# ------------------------------------------------------------
# EIC 'source' kernel
def kerS2d(ker,uS,uR,res,cip,vel,custom,par):
    eic2dA(ker+'_gS',res,uR,cip,'gaus=y wflcausal=y oprcausal=n'+' '+custom,par)
    aweop2d(ker+'_aS',ker+'_gS',vel,'adj=y'+' '+custom,par)
    cic2dF(ker,uS,ker+'_aS','wflcausal=y oprcausal=y'+' '+custom,par)

def kerS3d(ker,uS,uR,res,cip,vel,custom,par):
    eic3dA(ker+'_gS',res,uR,cip,'gaus=y wflcausal=y oprcausal=n'+' '+custom,par)
    aweop3d(ker+'_aS',ker+'_gS',vel,'adj=y'+' '+custom,par)
    cic3dF(ker,uS,ker+'_aS','wflcausal=y oprcausal=y'+' '+custom,par)
    
# ------------------------------------------------------------
# EIC 'receiver' kernel
def kerR2d(ker,uS,uR,res,cip,vel,custom,par):
    Flow(  res+'_',res,'reverse which=7 opt=i')
    eic2dA(ker+'_gR',res+'_',uS,cip,'gaus=y wflcausal=y oprcausal=y'+' '+custom,par)
    aweop2d(ker+'_aR',ker+'_gR',vel,'adj=n'+' '+custom,par)
    cic2dF(ker,uR,ker+'_aR','wflcausal=y oprcausal=n'+' '+custom,par)

def kerR3d(ker,uS,uR,res,cip,vel,custom,par):
    Flow(  res+'_',res,'reverse which=15 opt=i')
    eic3dA(ker+'_gR',res+'_',uS,cip,'gaus=y wflcausal=y oprcausal=y'+' '+custom,par)
    aweop3d(ker+'_aR',ker+'_gR',vel,'adj=n'+' '+custom,par)
    cic3dF(ker,uR,ker+'_aR','wflcausal=y oprcausal=n'+' '+custom,par)
    
# ------------------------------------------------------------
# EIC total kernel
def kerT2d(ker,uS,uR,res,cip,vel,custom,par):
    kerS2d(ker+'_S',uS,uR,res,cip,vel,custom,par)
    kerR2d(ker+'_R',uS,uR,res,cip,vel,custom,par)
    Flow(ker,[ker+'_S',ker+'_R'],'add ${SOURCES[1]}')

def kerT3d(ker,uS,uR,res,cip,vel,custom,par):
    kerS3d(ker+'_S',uS,uR,res,cip,vel,custom,par)
    kerR3d(ker+'_R',uS,uR,res,cip,vel,custom,par)
    Flow(ker,[ker+'_S',ker+'_R'],'add ${SOURCES[1]}')

    
