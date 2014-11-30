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
    if(not par.has_key('fdorder')):    par['fdorder']=8
    if(not par.has_key('optfd')):    par['optfd']='y'
    if(not par.has_key('hybridbc')):    par['hybridbc']='y'

    if(not par.has_key('gaus')):     par['gaus']='y'

# ------------------------------------------------------------
def awepar(par):
    awe = ' ' + \
          '''
          ompchunk=%(ompchunk)d ompnth=%(ompnth)d
          verb=%(verb)s fsrf=%(fsrf)s
          dabc=%(dabc)s nb=%(nb)d
          snap=%(snap)s jsnap=%(jsnap)d
          '''%par + ' '
    return awe

def eicpar(par):
    eic = ' ' + \
          '''
          nhx=%(nhx)d nhy=%(nhy)d nhz=%(nhz)d nht=%(nht)d
          gaus=%(gaus)s verb=%(verb)s
          '''%par + ' '
    return eic

def iwindow(par):
    win = ' ' + \
          '''
          nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g 
          nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g
          ''' % par + ' '
    return win

def aweoptpar(par):
    aweopt = ' ' + \
          '''
          verb=%(verb)s
          fdorder=%(fdorder)d optfd=%(optfd)s
          dabc=%(dabc)s nb=%(nb)d hybridbc=%(hybridbc)s 
          snap=%(snap)s jsnap=%(jsnap)d
          '''%par + ' '
    return aweopt

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
# constant-density acoustic FD modeling with optimized fd and hybrid bc
def cdafd2dopt(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd2dopt
         vel=${SOURCES[1]}
         sou=${SOURCES[2]} rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         ''' + ' ' + aweoptpar(par) + ' ' + custom)
         
def cdafd3dopt(odat,owfl,idat,velo,sou,rec,custom,par):    
    Flow([odat,owfl],[idat,velo,sou,rec],
         '''
         awefd3dopt
         vel=${SOURCES[1]}
         sou=${SOURCES[2]} rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         ''' + ' ' + aweoptpar(par) + ' ' + custom)

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

# ------------------------------------------------------------
# variable-density RTM w/ CIC
def cicmig(icic,
           sdat,scoo,
           rdat,rcoo,
           velo,dens,
           custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=icic+'swfl'
    rdrv=icic+'rdrv'
    rwfl=icic+'rwfl'

    Flow(icic,[sdat,scoo,rdat,rcoo,velo,dens],
         '''
         %sawefd2d < ${SOURCES[0]} cden=n %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[4]}
         den=${SOURCES[5]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+awepar(par)+custom,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=n %s verb=n
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         vel=${SOURCES[4]}
         den=${SOURCES[5]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+awepar(par)+custom,rwfl,DPT) +
         '''
         %scicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)

# ------------------------------------------------------------
# constant-density RTM w/ CIC
def cicmigCD(icic,
             sdat,scoo,
             rdat,rcoo,
             velo,
             custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=icic+'swfl'
    rdrv=icic+'rdrv'
    rwfl=icic+'rwfl'

    Flow(icic,[sdat,scoo,rdat,rcoo,velo],
         '''
         %sawefd2d < ${SOURCES[0]} cden=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[4]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+awepar(par)+custom,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=y %s verb=n
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         vel=${SOURCES[4]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+awepar(par)+custom,rwfl,DPT) +
         '''
         %scicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)
    
# ------------------------------------------------------------
# veriable-density RTM w/ CIC and EIC
def eicmig(icic,
           ieic,icoo,
           sdat,scoo,
           rdat,rcoo,
           velo,dens,
           custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=ieic+'swfl'
    rdrv=ieic+'rdrv'
    rwfl=ieic+'rwfl'

    Flow([icic,ieic],[sdat,scoo,rdat,rcoo,icoo,velo,dens],
         '''
         %sawefd2d < ${SOURCES[0]} cden=n %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[5]}
         den=${SOURCES[6]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+awepar(par)+custom,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=n %s verb=n
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         vel=${SOURCES[5]}
         den=${SOURCES[6]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+awepar(par)+custom,rwfl,DPT) +
         '''
         %scicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %seicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s cip=${SOURCES[4]} 
         >${TARGETS[1]};
         '''%(M8R,swfl,eicpar(par)+custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)

# ------------------------------------------------------------
# constant-density RTM w/ CIC and EIC
def eicmigCD(icic,
             ieic,icoo,
             sdat,scoo,
             rdat,rcoo,
             velo,
             custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=ieic+'swfl'
    rdrv=ieic+'rdrv'
    rwfl=ieic+'rwfl'

    Flow([icic,ieic],[sdat,scoo,rdat,rcoo,icoo,velo],
         '''
         %sawefd2d < ${SOURCES[0]} cden=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[5]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,iwindow(par)+awepar(par)+custom,swfl,DPT) +
         '''
         %sreverse < ${SOURCES[2]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdrv,DPT) +
         '''
         %sawefd2d < %s cden=y %s verb=n
         sou=${SOURCES[3]}
         rec=${SOURCES[3]}
         vel=${SOURCES[5]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdrv,iwindow(par)+awepar(par)+custom,rwfl,DPT) +
         '''
         %scicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s 
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %seicop2d <%s adj=n wflcausal=y oprcausal=n %s
         opr=%s cip=${SOURCES[4]} 
         >${TARGETS[1]};
         '''%(M8R,swfl,eicpar(par)+custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)
    
# ------------------------------------------------------------
# zero-offset RTM - variable density
def awertm2d(imag,data,rcoo,velo,dens,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo,dens],
         '''
         %sawefd2d < ${SOURCES[0]} adj=y cden=n %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         den=${SOURCES[3]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)

def awertm3d(imag,data,rcoo,velo,dens,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo,dens],
         '''
         %sawefd3d < ${SOURCES[0]} adj=y cden=n %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         den=${SOURCES[3]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n4=1 f4=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)
    
# zero-offset RTM - constant density
def cdartm2d(imag,data,rcoo,velo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sawefd2d < ${SOURCES[0]} adj=y cden=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)

def cdartm3d(imag,data,rcoo,velo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sawefd3d < ${SOURCES[0]} adj=y cden=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n4=1 f4=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)
 # ------------------------------------------------------------
# zero-offset RTM with awefd2dopt/awefd3dopt - constant density
def cdaoptrtm2d(imag,data,rcoo,velo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sawefd2dopt < ${SOURCES[0]} adj=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)

def cdaoptrtm3d(imag,data,rcoo,velo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo],
         '''
         %sawefd3dopt < ${SOURCES[0]} adj=y %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         wfl=%s datapath=%s/ %s
         >/dev/null;
         '''%(M8R,awepar(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT,custom) +
         '''
         %swindow < %s n4=1 f4=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s
         '''%(M8R,rwfl),
              stdin=0,
              stdout=0)
# ------------------------------------------------------------
def dPAD2d(wfld,trac,ix,iz,par):
    Flow(wfld,trac,
        '''
        transp plane=23 |
        pad beg1=%d n1out=%d beg2=%d n2out=%d |
        put o1=%g d1=%g o2=%g d2=%g
        '''%(iz,par['nz'],ix,par['nx'],
             par['oz'],par['dz'],
             par['ox'],par['dx']))

def dWIN2d(trac,wfld,ix,iz,par):
    Flow(trac,wfld,
         '''
         window n1=1 f1=%d n2=1 f2=%d |
         transp
         '''%(iz,ix))

         
