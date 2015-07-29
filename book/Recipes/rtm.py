try:    from rsf.cluster import *
except: from rsf.proj    import *
import fdmod

# ------------------------------------------------------------
def awepar(par):
    awe = ' ' + \
          '''
          ompchunk=%(ompchunk)d ompnth=%(ompnth)d
          verb=y free=n
          snap=%(snap)s
          jsnap=%(jdata)d jdata=%(jdata)d
          dabc=%(dabc)s nb=%(nb)d
          '''%par + ' '
    return awe
    
def iwindow(par):
    win = ' ' + \
          '''
          nqz=%(nqz)d oqz=%(oqz)g dqz=%(dqz)g 
          nqx=%(nqx)d oqx=%(oqx)g dqx=%(dqx)g
          ''' % par + ' '
    return win

def eicpar(par):
    eic = ' ' + 'nhx=%(nhx)d nhy=%(nhy)d nhz=%(nhz)d nht=%(nht)d'%par + ' '
    return eic
    
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s 
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %seicold2d <%s isreversed=0 verb=n %s
         ur=%s cc=${SOURCES[4]} 
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s 
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %seicold2d <%s isreversed=0 verb=n %s
         ur=%s cc=${SOURCES[4]} 
         >${TARGETS[1]};
         '''%(M8R,swfl,eicpar(par)+custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)


# ------------------------------------------------------------
# CIC: cross-correlation
def cic(icic,swfl,rwfl,custom,par,isreversed=0):
    Flow(icic,[swfl,rwfl],
         '''
         cicold2d verb=y
         ur=${SOURCES[1]} 
         '''%par+ 'isreversed=%d'%isreversed + ' ' + custom )

# EIC
def eic(ieic,swfl,rwfl,cc,custom,par,isreversed=0):    
    Flow(ieic,[swfl,rwfl,cc],
         '''
         eicold2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         '''%par+ 'isreversed=%d'%isreversed + ' '+ custom )

# ------------------------------------------------------------
def zofmig(imag,data,rcoo,velo,dens,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

    rdat = imag+'rdat'
    rwfl = imag+'wwfl'

    Flow(imag,[data,rcoo,velo,dens],
         '''
         %sreverse < ${SOURCES[0]} which=2 opt=i verb=n >%s datapath=%s/;
         '''%(M8R,rdat,DPT) +
         '''
         %sawefd2d < %s cden=n %s verb=n
         sou=${SOURCES[1]}
         rec=${SOURCES[1]}
         vel=${SOURCES[2]}
         den=${SOURCES[3]}
         wfl=%s datapath=%s/
         >/dev/null;
         '''%(M8R,rdat,awepar(par)+iwindow(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s %s
         '''%(M8R,rdat,rwfl),
              stdin=0,
              stdout=0)
    
# ------------------------------------------------------------
def zofmigCD(imag,data,rcoo,velo,custom,par):
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH'))

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
         '''%(M8R,rdat,awepar(par)+iwindow(par)+' jsnap=%d'%(par['nt']-1),rwfl,DPT) +
         '''
         %swindow < %s n3=1 f3=1 >${TARGETS[0]};
         '''%(M8R,rwfl) +
         '''
         %srm %s %s
         '''%(M8R,rdat,rwfl),
              stdin=0,
              stdout=0)
    
# ------------------------------------------------------------
# FWI kernel
def fwiker(ker,dts,ss,dtr,rr,vel,den,custom,par):
    
    fdmod.awefd(ker+'_SD',ker+'_SW',dts,vel,den,ss,rr,custom+iwindow(par),par)

    Flow(dtr+'_R',dtr,'reverse which=2 opt=i verb=y')
    fdmod.awefd(ker+'_RD',ker+'_RW',dtr+'_R',vel,den,rr,ss,custom+iwindow(par),par)
    
    cic(ker,ker+'_SW',ker+'_RW','',par,isreversed=0)

def fwikerCDold(ker,dts,ss,dtr,rr,vel,custom,par):
    
    fdmod.cdafd(ker+'_SD',ker+'_SW',dts,vel,ss,rr,custom+iwindow(par),par)

    Flow(dtr+'_R',dtr,'reverse which=2 opt=i verb=y')
    fdmod.cdafd(ker+'_RD',ker+'_RW',dtr+'_R',vel,rr,ss,custom+iwindow(par),par)
    
    cic(ker,ker+'_SW',ker+'_RW','',par,isreversed=0)
    
# ------------------------------------------------------------
# FWI kernel (same as variable-density RTM w/ CIC)
def fwiker(kern,
           sdat,scoo,
           rdat,rcoo,
           velo,dens,
           custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=kern+'swfl'
    rdrv=kern+'rdrv'
    rwfl=kern+'rwfl'

    Flow(kern,[sdat,scoo,rdat,rcoo,velo,dens],
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)

# ------------------------------------------------------------
# FWI kernel (same as constant-density RTM w/ CIC)
def fwikerCD(kern,
             sdat,scoo,
             rdat,rcoo,
             velo,
             custom,par):
    
    M8R='$RSFROOT/bin/sf'
    DPT=os.environ.get('TMPDATAPATH')

    swfl=kern+'swfl'
    rdrv=kern+'rdrv'
    rwfl=kern+'rwfl'

    Flow(kern,[sdat,scoo,rdat,rcoo,velo],
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
         %scicold2d <%s isreversed=0 verb=n %s
         ur=%s
         >${TARGETS[0]};
         '''%(M8R,swfl,custom,rwfl) +
         '''
         %srm %s %s %s
         '''%(M8R,swfl,rdrv,rwfl),
              stdin=0,
              stdout=0)
