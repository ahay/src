from rsf.proj import *
import fdmod,pplot,pot,wemig,adcig

# ------------------------------------------------------------
# acoustic wave-equation modeling
def awefd(odat,owfl,idat,velo,dens,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [odat,owfl],[idat,velo,dens,sou,rec],
         '''
         sfawefd2d
         ompchunk=%(ompchunk)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# elastic wave-equation modeling
def ewefd(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
    par['fdcustom'] = custom

    print par['fdcustom']
    
    Flow( [odat,owfl],[idat,cccc,dens,sou,rec],
         '''
         sfewefd2dtti
         ompchunk=%(ompchunk)d 
         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d nbell=%(nbell)d
         ccc=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         %(fdcustom)s
         ''' % par)

# ------------------------------------------------------------
# heat equation modeling
def hdefd(dat,wfl,  wav,con,sou,rec,custom,par):
    par['fdcustom'] = custom
    
    Flow( [dat,wfl],[wav,con,sou,rec,'../code/HDEfd.x'],
          '''
          ../code/HDEfd.x
          verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
          con=${SOURCES[1]}
          sou=${SOURCES[2]}
          rec=${SOURCES[3]}
          wfl=${TARGETS[1]}
          %(fdcustom)s
          ''' % par)

# ------------------------------------------------------------
# acoustic reverse-time migration
def artm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    fdmod.awefd(sout,swfl,sdat,velo,dens,sacq,iacq,custom,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'
    temp = imag+'_tmp'
    
    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    fdmod.awefd(tout,rwfl,tdat,velo,dens,racq,iacq,custom,par)
    Flow(rout,tout,'reverse which=2 opt=i verb=y')
    
    # conventional (cross-correlation zero-lag) imaging condition
    Flow(temp,[sout,rout],'add ${SOURCES[1]} mode=p|stack axis=2')
    Flow(imag,temp,
         '''
         put
         n1=%d o1=%g d1=%g
         n2=%d o2=%g d2=%g
         ''' %(par['nqz'],par['oqz'],par['dqz'],
               par['nqx'],par['oqx'],par['dqx']) )

# ------------------------------------------------------------
def xzcoord(par):
    return '''
    put
    n2=%d o2=%g d2=%g
    n3=%d o3=%g d3=%g
    ''' %(par['nqz'],par['oqz'],par['dqz'],
          par['nqx'],par['oqx'],par['dqx'])
def xzcoord1(par):
    return '''
    put
    n1=%d o1=%g d1=%g
    n2=%d o2=%g d2=%g
    ''' %(par['nqz'],par['oqz'],par['dqz'],
          par['nqx'],par['oqx'],par['dqx'])
# ------------------------------------------------------------
# elastic reverse-time migration
def ertm(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    ewefd(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y " + custom,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'
    temp = imag+'_tmp'
    
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    ewefd(tout,rwfl,tdat,cccc,dens,racq,iacq," ssou=n " + custom,par)
    Flow(rout,tout,'reverse which=4 opt=i verb=y')

    # ------------------------------------------------------------
    Flow(sout+'1',sout,'window n2=1 f2=0')
    Flow(sout+'2',sout,'window n2=1 f2=1')

    Flow(rout+'1',rout,'window n2=1 f2=0')
    Flow(rout+'2',rout,'window n2=1 f2=1')    
    
    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(temp+i+j,[sout+i,rout+j],
                 'add ${SOURCES[1]} mode=p | stack axis=2')

            Flow(imag+i+j,temp+i+j,
                 '''
                 put
                 n1=%d o1=%g d1=%g
                 n2=%d o2=%g d2=%g
                 ''' %(par['nqz'],par['oqz'],par['dqz'],
                       par['nqx'],par['oqx'],par['dqx']) )

    Flow(imag,[imag+'11',imag+'12',imag+'21',imag+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def Aertm(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    ewefd(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y " + custom,par)
   
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'
    temp = imag+'_tmp'
    
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    ewefd(tout,rwfl,tdat,cccc,dens,racq,iacq," ssou=n " + custom,par)
    # separate wave modes for receiver wavefield
    
    Flow(rout,tout,'reverse which=4 opt=i verb=y')

    # ------------------------------------------------------------
    Flow(sout+'1',sout,'''window n2=1 f2=0 | put
                 n1=%d o1=%g d1=%g
                 n2=%d o2=%g d2=%g
                 n3=%d
                 ''' %(par['nqz'],par['oqz'],par['dqz'],
                       par['nqx'],par['oqx'],par['dqx'],par['nt']/par['jdata']) )
    Flow(sout+'2',sout,'''window n2=1 f2=1| put
                 n1=%d o1=%g d1=%g
                 n2=%d o2=%g d2=%g
                 n3=%d
                 ''' %(par['nqz'],par['oqz'],par['dqz'],
                       par['nqx'],par['oqx'],par['dqx'],par['nt']/par['jdata']) )
    # separate wave modes for source wavefield
    # use sfconvolveframes2, modify pot.potentials if necessary
    pot.potentials(   'pA','uAz','uAx','dz','dx','n','','q',par)
    
    Flow(rout+'1',rout,'''window n2=1 f2=0| put
                 n1=%d o1=%g d1=%g
                 n2=%d o2=%g d2=%g
                 n3=%d
                 ''' %(par['nqz'],par['oqz'],par['dqz'],
                       par['nqx'],par['oqx'],par['dqx'],par['nt']/par['jdata']) )
    Flow(rout+'2',rout,'''window n2=1 f2=1| put
                 n1=%d o1=%g d1=%g
                 n2=%d o2=%g d2=%g
                 n3=%d
                 ''' %(par['nqz'],par['oqz'],par['dqz'],
                       par['nqx'],par['oqx'],par['dqx'],par['nt']/par['jdata']) ) 
    
    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(imag+i+j,[sout+i,rout+j],
                 'add ${SOURCES[1]} mode=p | stack axis=3')

           

    Flow(imag,[imag+'11',imag+'12',imag+'21',imag+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')
# ------------------------------------------------------------
# stereographic imaging condition
def esic(isic,imag,custom,par):
    par['sic']=custom

    qs   = isic + '_qs'
    qr   = isic + '_qr'

    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)
    
    Flow(qs+'1',sout,'window n2=1 f2=0 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
    Flow(qs+'2',sout,'window n2=1 f2=1 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
    Flow(qr+'1',rout,'window n2=1 f2=0 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
    Flow(qr+'2',rout,'window n2=1 f2=1 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
 
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(isic+i+j,[qs+i,qr+j],
                 '''
                 sic3d ur=${SOURCES[1]} nbuf=100 verb=y
                 oa=%(oa)g na=%(na)d da=%(da)g
                 %(sic)s
                 nl=%(nl)d dl=%(dl)g
                 sig=%(sig)g
                 ''' % par)
    Flow(isic,[isic+'11',isic+'12',isic+'21',isic+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# plot elastic wavefield
def emovie(wfld,custom,axis,par):

    # loop over wavefield components
    for i in range(2):
        Flow(wfld+str(i+1),wfld,
             '''
             window n3=1 f3=%d
             '''%i)
              #|
             #window min1=%g max1=%g min2=%g max2=%g
             #''' % (i,par['zmin'],par['zmax'],par['xmin'],par['xmax']))

    # join component wavefields
    Flow(wfld+'all',[wfld+'1',wfld+'2'],
         'cat axis=%d space=n ${SOURCES[1:2]}' % axis)

    if(axis==1):        
        height=2*par['height']
        if(height>10): height=10
        ratio =2*par['ratio']
    if(axis==2):
        height=par['height']
        if(height>10): height=10
        ratio =0.5*par['ratio']
    
    Plot(wfld,wfld+'all',
         '''
         grey title="" wantaxis=y screenratio=%f screenht=%f
         gainpanel=a pclip=99 %s
         %s
         ''' % (ratio,height,par['labelattr'],custom),view=1)

# ------------------------------------------------------------
# plot elastic image
def eimage(plot,imag,custom,par):
#    print clip
    title=['pp','ps','sp','ss']
    for i in range(4):
        if(i!=0):
            flag=' wantaxis2=n'
        else:
            flag=' wantaxis2=y'
        
        Flow([plot+'_plt'+str(i),plot+'_bar'+str(i)],imag,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a')        

        Plot(plot+str(i),[plot+'_plt'+str(i),plot+'_bar'+str(i)],
             'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
             + fdmod.cgrey(flag,par))
        Result(plot+str(i),[plot+'_plt'+str(i),plot+'_bar'+str(i)],
             'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
             + fdmod.cgrey(flag+custom+' wantaxis2=y title=%s'% (title[i]),par))

# ------------------------------------------------------------
# plot elastic data
def edata(plot,data,custom,par):

    Flow([plot+'_plt',plot+'_bar'],data,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=100 %s' % custom)

    for i in range(2):
        Plot(  plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
               +fdmod.dgrey('pclip=98 %s' %custom,par))        
        Result(plot+str(i+1),[plot+'_plt',plot+'_bar'],
               'window n2=1 f2=%d j3=10 bar=${SOURCES[1]} | transp |' % i
               +fdmod.dgrey('pclip=98 %s' %custom,par))   

# ------------------------------------------------------------
# plot elastic wavelet
def ewavelet(wavelet,custom,par):

    for i in range(2):
        Plot(wavelet+str(i+1),wavelet,
             'window n2=1 f2=%d | transp | window |'%i +
             fdmod.waveplot('%d'%i,par))
    Result(wavelet,[wavelet+'1',wavelet+'2'],'Movie')
        
# ------------------------------------------------------------
def wcompare(plot,us,ur,iz,par):

    for u in ([us,ur]):
        for i in range(2):
            Flow(plot+'_'+u+str(i+1),u,
                 'window n2=1 f2=%d j3=10 | transp |' % i
                 + xzcoord(par)
                 + '| window n2=1 f2=%d' % iz)
            
        Flow(u+'all',[plot+'_'+u+'1',plot+'_'+u+'2'],
             'cat axis=2 space=n ${SOURCES[0:2]}', stdin=0)
        Plot(u+'all',fdmod.dgrey('max2=%g' % (2*par['xmax']) ,par))

    Result(plot,[us+'all',ur+'all'],'SideBySideAniso')

# ------------------------------------------------------------

# ------------------------------------------------------------

def hgrey(custom):
    return '''
    grey pclip=100 title="cig" screenratio=2
    label1=z unit1=m label2=l unit2=m labelsz=6 %s
    ''' % custom
def agrey(custom):
    return '''
    grey pclip=100 title="cag" screenratio=2 grid=y
    label1=z unit1=m label2="\F10 q\F3 " unit2="\^o\_" labelsz=6 %s
    ''' % custom

# CIG to slant-stack
def tcig2ssk(par):
    return '''
    slant adj=y p0=1000 np=600 dp=10 |
    put label2=v
    '''%par
def hcig2ssk(par):
    return '''
    slant adj=y p0=-0.5 np=300 dp=0.01 |
    put label2=tan unit2=""
    '''%par
# slant-stack to angle
def tssk2ang(par):
    return '''
    pp2pstsic a0=%(oa)g na=%(na)d da=%(da)g
    velocity=${SOURCES[1]}
    vpvs=${SOURCES[2]}
    dip=${SOURCES[3]} |
    put label2=ang
    ''' % par
def hssk2ang_pos(par):
    return '''
    pp2psang2
    vpvs=${SOURCES[1]}
    dip=${SOURCES[2]} |
    tan2ang a0=%(oap)d na=%(nap)d da=%(dap)g |
    put label2=ang
    ''' % par
def hssk2ang(par):
    return '''
    pp2psang2
    vpvs=${SOURCES[1]}
    dip=${SOURCES[2]} |
    tan2ang a0=%(oa)d na=%(na)d da=%(da)g |
    put label2=ang
    ''' % par
def pssk2ang(par):
    return '''
    tan2ang a0=-80 na=401 da=0.4
    ''' %par
# ------------------------------------------------------------
def alaps(ii,us,ur,nhz,nhx,nht,n,par):

    # reshape inputs us and ur to z-x-t
    Flow(ii+'_us',us,'transp plane=23 |' + xzcoord1(par))
    Flow(ii+'_ur',ur,'transp plane=23 |' + xzcoord1(par))

    # E.I.C. (space-lags)
    Flow(ii,['ca_us','ca_ur'],
         '''
         laps verb=y 
         ur=${SOURCES[1]}
         nhz=%d nhx=%d nht=%d |
         put label1=z label2=x unit1=km unit2=km unit3=km unit4=km unit5=s
         ''' % (nhz,nhx,nht))
    
    
    # angle-gather (one location under shot)
    for i in range(n):
      x=par['xsou']   
      xm=int(x*1000)
      if (x<par['dx']*(par['nx']-1)+par['ox']):
          Result(ii+str(xm),ii,'window  min2=%d n2=1 |'%x + hcig2ssk() + '|' + hssk2ang(par) + '|' + agrey(''))
  

    # angle-gather (horizontal space stack = shots at all positions,  this uses the trick of summing lags gathers at all spatial positions which only works in v(z) and with horizontal reflectors)
    Flow( ii+'_cig',ii,'window | stack')
    Result(ii+'_cig',hgrey(''))
    
    Flow(  ii+'_ang',ii+'_cig',
           hcig2ssk() + '|' + hssk2ang(par) +'put label2=theta unit2=degree')
    Result(ii+'_ang',agrey(''))

# ------------------------------------------------------------
def elaps(ii,us,ur,nhx,nhz,nht,n,dipa,par):
    
    # velocity ratio
    Flow('vratio11',None,'math v1=vp.rsf v2=vp.rsf output=v1/v2')
    Flow('vratio12',None,'math v1=vp.rsf v2=vs.rsf output=v1/v2')
    Flow('vratio21',None,'math v1=vs.rsf v2=vp.rsf output=v1/v2')
    Flow('vratio22',None,'math v1=vs.rsf v2=vs.rsf output=v1/v2')  

    
    title1=['P','S'];title2=['P','S'];
    for i in (1,1):
       for j in (1,2):  
       	# reshape inputs us and ur to z-x-t
   	  Flow(ii+'_us'+str(i),us+str(i),' transp plane=23 |' + xzcoord1(par))
   	  Flow(ii+'_ur'+str(j),ur+str(j),' transp plane=23 |' + xzcoord1(par))

    	# E.I.C. (space-lags)
    	  Flow(ii+str(i)+str(j),[ii+'_us'+str(i),ii+'_ur'+str(j)],
         	'''
         	laps verb=y 
         	ur=${SOURCES[1]}
        	nhz=%d nhx=%d nht=%d |
        	put label1=z label2=x 
        	    unit1=km unit2=km unit3=km unit4=km unit5=s
        	''' % (nhz,nhx,nht))
	# angle gather
   	  #x=par['xsou']
   	  x=1.4
    	  xm=int(x*1000)
#    	  print x,xm
    	  if (x<par['dx']*(par['nx']-1)+par['ox']):
              Flow(ii+str(i)+str(j)+'cig_'+str(xm),ii+str(i)+str(j),'window  min2=%g n2=1 '%x) 
              Flow(ii+str(i)+str(j)+'ssk_'+str(xm),ii+str(i)+str(j)+'cig_'+str(xm),hcig2ssk(par) )
              #Flow(ii+str(i)+str(j)+'ssk_'+str(xm),ii+str(i)+str(j),'window  min2=%g n2=1 '%x + hcig2ssk() )               		
  	      Flow(ii+str(i)+str(j)+'ang_'+str(xm),[ii+str(i)+str(j)+'ssk_'+str(xm),'vratio'+str(i)+str(j),'dipa'], hssk2ang_pos(par) )
  	      Result( ii+str(i)+str(j)+'ang_'+str(xm),agrey('title="%s%s"'%(title1[i-1],title2[j-1])))
          Flow( ii+str(i)+str(j)+'_cig',ii+str(i)+str(j),'window | stack') 
          Result(ii+str(i)+str(j)+'_cig',hgrey(''))
  
    
      
          # angle-gather (horizontal space stack = shots at all positions,  this uses the trick of summing lags gathers at all spatial positions which only works in v(z) and with horizontal reflectors)          
          Flow(ii+str(i)+str(j)+'_ssk',ii+str(i)+str(j)+'_cig',hcig2ssk(par) )            		
  	  Flow(ii+str(i)+str(j)+'_ang',[ii+str(i)+str(j)+'_ssk','vratio'+str(i)+str(j),'dipa'], hssk2ang(par) )
          Result(ii+str(i)+str(j)+'_ang',agrey('grid=n title="%s%s"'%(title1[i-1],title2[j-1]))) 
          

   
def laps(ii,us,ur,cc,nhz,nhx,nht,custom,par):
    Flow(ii,[us,ur,cc],
         '''
         laps2d verb=y 
         ur=${SOURCES[1]} cc=${SOURCES[2]}
         nhz=%d nhx=%d nht=%d    %s |
         put n1=%d n2=%d n3=%d n4=1 
             label1=hx label2=ccx label3=ccz
             d2=%g d3=%g o2=%g o3=%g
         ''' % (nhz,nhx,nht,custom,
                nhx*2+1,
                par['nqx'],par['nqz'],
                par['dqx'],par['dqz'],
                par['oqx'],par['oqz']))
    

# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------
# ------------------------------------------------------------
# elastic reverse-time migration
def ewfld(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # source wavefield (z,x,t)
    ewefd(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y " + custom,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    tout = imag+'_tdr'
    temp = imag+'_tmp'
    
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    ewefd(tout,rwfl,tdat,cccc,dens,racq,iacq," ssou=n " + custom,par)
    Flow(rout,tout,'reverse which=4 opt=i verb=y')

    # ------------------------------------------------------------
#    Flow(sout+'1',sout,'window n2=1 f2=0')
#    Flow(sout+'2',sout,'window n2=1 f2=1')
#
#    Flow(rout+'1',rout,'window n2=1 f2=0')
#    Flow(rout+'2',rout,'window n2=1 f2=1')    
    
    par['ntnew']=par['nt']/par['jdata']
    print par['ntnew'],par['nt'],par['jsnap']
#   re-organize source wavefield for the following imaging conditions
    for i in range(1,3,1):
	print "comp",i
	par['comp']=i-1
	tag="%01d"%i
	print sout+tag
        Flow(sout+tag,sout,
 	    '''window n2=1 f2=%(comp)d | 
               put n1=%(nqz)d n2=%(nqx)d n3=%(ntnew)d 
                   d1=%(dqz)g d2=%(dqx)g d3=1
	           o1=%(oqz)g o2=%(oqx)g o3=0'''%par)
#   re-organize receiver wavefield for the following imaging conditions
    for i in range(1,3,1):
        par['comp']=i-1
	tag="%01d"%i
        Flow(rout+tag,rout,
            '''window n2=1 f2=%(comp)d | 
               put n1=%(nqz)d n2=%(nqx)d n3=%(ntnew)d 
                   d1=%(nqz)d d2=%(dqx)d d3=1
                   o1=%(oqz)g o2=%(oqx)g o3=0'''%par)

    Flow('vratio11','zero','math output=1')
    Flow('vratio12','zero','math output=2')
    Flow('vratio21','zero','math output=.5')
    Flow('vratio22','zero','math output=1')


def ecic(imag,ss,rr,cc,custom,par):
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
	      wemig.cic(imag+i+j,sout+i,rout+j,'nbuf=100',par)
	      Plot(imag+i+j,fdmod.cgrey('pclip=99.8'+custom,par))
              Result(imag+i+j,[imag+i+j,'ss','rr'],'Overlay')

    Flow(imag,[imag+'11', imag+'12',imag+'21',imag+'22'],
         'cat ${SOURCES[1:4]} space=n axis=3')

def eeic(imag,ss,rr,cc,xcig,custom,par):
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1'):
        for j in ('1','2'):
            wemig.eic(imag+'Ecip'+i+j,sout+i,rout+j,'cc','',par)	
            Flow(imag+'Elaps'+i+j,imag+'Ecip'+i+j,
                 '''put n4=%(ncz)d n5=%(ncx)d 
                o4=%(ocz)g o5=%(ocx)g 
                d4=%(dcz)g d5=%(dcx)g 
                label4=ccz label5=ccx |
                window 
               '''%par)
            Flow(imag+'Ecig'+i+j,imag+'Elaps'+i+j,'window n3=1 min3=%g | transp'%xcig )
            Flow(imag+'Eang'+i+j,[imag+'Ecig'+i+j,'zero','vratio'+i+j],
                 adcig.cig2ssk(300,-1.5,0.01) + '|' +
                 adcig.xsk2ang(320, -80,0.50))
            Result(imag+'Eang'+i+j,adcig.agrey(' grid=y  ',par))
            
#            Flow('Ecigall'+i+j,'Elaps'+i+j,'window j3=3 | transp ')
#            Flow('Eangall'+i+j,['Ecigall'+i+j,'zero','vratio'+i+j],
#                 adcig.cig2ssk(300,-1.5,0.01) + '|' +
#                 adcig.xsk2ang(320, -80,0.50) )
#            Result('Eangall'+i+j,'stack axis=3|'+adcig.agrey(' grid=y  ',par))
