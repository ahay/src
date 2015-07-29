try:    from rsf.cluster import *
except: from rsf.proj    import *
import fdmod,pplot

# ------------------------------------------------------------
# acoustic wave-equation modeling
#def awefd(odat,owfl,idat,velo,dens,sou,rec,custom,par):
#    par['fdcustom'] = custom
#    
#    Flow( [odat,owfl],[idat,velo,dens,sou,rec,'../code/AWEfd.x'],
#         '''
#         ../code/AWEfd.x
#         ompchunk=%(ompchunk)d 
#         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
#         vel=${SOURCES[1]}
#         den=${SOURCES[2]}
#         sou=${SOURCES[3]}
#         rec=${SOURCES[4]}
#         wfl=${TARGETS[1]}
#         %(fdcustom)s
#         ''' % par)

# ------------------------------------------------------------
# elastic wave-equation modeling
#def ewefd(odat,owfl,idat,cccc,dens,sou,rec,custom,par):
#    par['fdcustom'] = custom
#
#    print par['fdcustom']
#    
#    Flow( [odat,owfl],[idat,cccc,dens,sou,rec,'../code/EWEfd.x'],
#         '''
#         ../code/EWEfd.x
#         ompchunk=%(ompchunk)d 
#         verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d nbell=%(nbell)d
#         ccc=${SOURCES[1]}
#         den=${SOURCES[2]}
#         sou=${SOURCES[3]}
#         rec=${SOURCES[4]}
#         wfl=${TARGETS[1]}
#         %(fdcustom)s
#         ''' % par)

# ------------------------------------------------------------
# heat equation modeling
#def hdefd(dat,wfl,  wav,con,sou,rec,custom,par):
#    par['fdcustom'] = custom
#    
#    Flow( [dat,wfl],[wav,con,sou,rec,'../code/HDEfd.x'],
#          '''
#          ../code/HDEfd.x
#          verb=y free=n snap=%(snap)s jsnap=%(jsnap)d nb=%(nb)d
#          con=${SOURCES[1]}
#          sou=${SOURCES[2]}
#          rec=${SOURCES[3]}
#          wfl=${TARGETS[1]}
#          %(fdcustom)s
#          ''' % par)

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
# acoustic reverse-time migration
def artm(imag,sdat,rdat,velo,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    iwindow = ' ' + \
    '''
    nqz=%(nqz)d oqz=%(oqz)g
    nqx=%(nqx)d oqx=%(oqx)g
    jsnap=%(jdata)d jdata=%(jdata)d
    ''' % par + ' '
    
    # source wavefield (z,x,t)
    fdmod.awefd(sout,swfl,sdat,velo,dens,sacq,iacq,custom+iwindow,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=2 opt=i verb=y')
    fdmod.awefd(rout,twfl,tdat,velo,dens,racq,iacq,custom+iwindow,par)
    Flow(rwfl,twfl,'reverse which=4 opt=i verb=y')
    
    # conventional (cross-correlation zero-lag) imaging condition
    Flow(imag,[swfl,rwfl],'xcor2d uu=${SOURCES[1]} axis=3 verb=y nbuf=100')

# ------------------------------------------------------------
# elastic reverse-time migration
def ertm(imag,sdat,rdat,cccc,dens,sacq,racq,iacq,custom,par):

    swfl = imag+'_us' #   source wavefield
    rwfl = imag+'_ur' # receiver wavefield
    twfl = imag+'_ut' #     temp wavefield
    sout = imag+'_ds' #   source data (not the input sdat!)
    rout = imag+'_dr' # receiver data (not the input rdat!)

    iwindow = ' ' + \
    '''
    nqz=%(nqz)d oqz=%(oqz)g
    nqx=%(nqx)d oqx=%(oqx)g
    jsnap=%(jdata)d jdata=%(jdata)d
    ''' % par + ' '

    # source wavefield (z,x,t)
    fdmod.ewefd(sout,swfl,sdat,cccc,dens,sacq,iacq," ssou=y opot=n" + custom+iwindow,par)
    
    # receiver wavefield (z,x,t)
    tdat = imag+'_tds'
    Flow(tdat,rdat,'reverse which=4 opt=i verb=y')
    fdmod.ewefd(rout,twfl,tdat,cccc,dens,racq,iacq," ssou=n opot=n" + custom+iwindow,par)
#    Flow(rwfl,twfl,'reverse which=8 opt=i verb=y')

    # ------------------------------------------------------------
    Flow(swfl+'1',swfl,'window n2=1 f2=0')
    Flow(swfl+'2',swfl,'window n2=1 f2=1')

    Flow(rwfl+'1',twfl,'window n2=1 f2=0 | reverse which=4 opt=i verb=y')
    Flow(rwfl+'2',twfl,'window n2=1 f2=1 | reverse which=4 opt=i verb=y')    
    
    # conventional (cross-correlation zero-lag) imaging condition
    for     i in ('1','2'):
        for j in ('1','2'):
            Flow(imag+i+j,[swfl+i,rwfl+j],
                 'xcor2d uu=${SOURCES[1]} axis=3 verb=y nbuf=100')

    Flow(imag,[imag+'11',imag+'12',imag+'21',imag+'22'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# stereographic imaging condition
#def esic(isic,imag,custom,par):
#    par['sic']=custom
#
#    qs   = isic + '_qs'
#    qr   = isic + '_qr'
#
#    sout = imag+'_ds' #   source data (not the input sdat!)
#    rout = imag+'_dr' # receiver data (not the input rdat!)
#    
#    Flow(qs+'1',sout,'window n2=1 f2=0 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
#    Flow(qs+'2',sout,'window n2=1 f2=1 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
#    Flow(qr+'1',rout,'window n2=1 f2=0 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
#    Flow(qr+'2',rout,'window n2=1 f2=1 min3=0.1 max3=0.4 | transp plane=12 memsize=500 |'+ xzcoord(par) + '| transp plane=23 memsize=500')
# 
#    for     i in ('1','2'):
#        for j in ('1','2'):
#            Flow(isic+i+j,[qs+i,qr+j],
#                 '''
#                 sic3d ur=${SOURCES[1]} nbuf=100 verb=y
#                 oa=%(oa)g na=%(na)d da=%(da)g
#                 %(sic)s
#                 nl=%(nl)d dl=%(dl)g
#                 sig=%(sig)g
#                 ''' % par)
#    Flow(isic,[isic+'11',isic+'12',isic+'21',isic+'22'],
#         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
# plot elastic wavefield
def emovie(wfld,custom,axis,par):

    # loop over wavefield components
    for i in range(2):
        Flow(wfld+str(i+1),wfld,
             '''
             window n3=1 f3=%d |
             window min1=%g max1=%g min2=%g max2=%g
             ''' % (i,par['zmin'],par['zmax'],par['xmin'],par['xmax']))

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
    
    Result(wfld,wfld+'all',
           '''
           grey title="" wantaxis=y screenratio=%f screenht=%f
           gainpanel=a pclip=99 %s
           %s
           ''' % (ratio,height,par['labelattr'],custom) )

# ------------------------------------------------------------
# plot elastic image
def eimage(plot,imag,clip,par):
    print clip
    title=['pp','ps','sp','ss']
    for i in range(4):
        if(i!=0):
            flag=' wantaxis2=n'
        else:
            flag=' wantaxis2=y'
        
        Flow([plot+'_plt'+str(i),plot+'_bar'+str(i)],imag,
         'scale axis=123 | byte bar=${TARGETS[1]} gainpanel=a pclip=%g'%clip[i])        

        Plot(plot+str(i),[plot+'_plt'+str(i),plot+'_bar'+str(i)],
             'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
             + fdmod.cgrey(flag,par))
        Result(plot+str(i),[plot+'_plt'+str(i),plot+'_bar'+str(i)],
             'window n3=1 f3=%d bar=${SOURCES[1]} |'% i
             + fdmod.cgrey(flag+' wantaxis2=y title=%s'% (title[i]),par))

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
               'window n2=1 f2=%d bar=${SOURCES[1]} | transp |' % i
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
    grey pclip=100 title="" screenratio=2
    label1=z unit1=m label2="\F10 l\F3 " unit2="m" labelsz=6 %s
    ''' % custom
def agrey(custom):
    return '''
    grey pclip=100 title="" screenratio=2
    label1=z unit1=a label2="\F10 q\F3 " unit2="\^o\_" labelsz=6 %s
    ''' % custom
def hcig2ssk():
    return '''
    slant adj=y p0=-4 np=400 dp=0.02
    '''
def hssk2ang():
    return '''
    tan2ang a0=-80 na=401 da=0.4
    '''
# ------------------------------------------------------------
# E.I.C. (space-lags)
def laps(ii,us,ur,nhz,nhx,nht,custom,par):
    Flow(ii,[us,ur],
         '''
         laps verb=y 
         ur=${SOURCES[1]}
         nhz=%d nhx=%d nht=%d
         %s
         ''' % (nhz,nhx,nht,custom))
    
