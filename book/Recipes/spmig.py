try:    from rsf.cluster import *
except: from rsf.proj    import *
import zomig, os

# ------------------------------------------------------------
def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if(par.has_key('verb')):
        p = p + ' verb=' +     par['verb']
    if(par.has_key('nrmax')):
        p = p + ' nrmax='+ str(par['nrmax'])
    if(par.has_key('dtmax')):
        p = p + ' dtmax='+ str(par['dtmax'])
    if(par.has_key('tmx')):
        p = p + ' tmx='  + str(par['tmx'])
    if(par.has_key('tmy')):
        p = p + ' tmy='  + str(par['tmy'])
    if(par.has_key('pmx')):
        p = p + ' pmx='  + str(par['pmx'])
    if(par.has_key('pmy')):
        p = p + ' pmy='  + str(par['pmy'])
    if(par.has_key('ompnth')):
        p = p + ' ompnth='  + str(par['ompnth'])
    if(par.has_key('misc')):
        p = p + ' '      +     par['misc']
    p = p + ' '
    return p

# ------------------------------------------------------------
def freqs(par):
    f  = ' '
    if(par.has_key('nw')): f = f + ' nw=' + str(par['nw'])
    if(par.has_key('ow')): f = f + ' ow=' + str(par['ow'])
    if(par.has_key('dw')): f = f + ' dw=' + str(par['dw'])
    f = f + ' '
    return f

def migpar(par):
    if(not par.has_key('verb')):    par['verb']='y'
    if(not par.has_key('eps')):     par['eps']=0.1

    if(not par.has_key('nrmax')):   par['nrmax']=1
    if(not par.has_key('dtmax')):   par['dtmax']=0.00005

    if(not par.has_key('tmx')):     par['tmx']=16
    if(not par.has_key('tmy')):     par['tmy']=16
    
    if(not par.has_key('incore')):  par['incore']='y'

# ------------------------------------------------------------
# create surface wavefield files
def wflds(swfl,rwfl,wave,shot,par):
    #
    _wave = swfl + '_' + wave
    _shot = rwfl + '_' + shot
    #
    _ssss = '_' + swfl
    _rrrr = '_' + rwfl
    
    # _wave(nw)
    Flow(_wave,wave,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         put label1=w
         ''' % par )

    Flow(_shot,shot,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         spray axis=3 n=1 o=0 d=1 |
         spray axis=5 n=1 o=0 d=1 |
         put label1=w label2=rx label3=ry label4=sx label5=sy
         ''' % par )

    # _rrrr(nw,nx,ny,ne)
    # _ssss(nw,nx,nt,ne)
    Flow([_rrrr,_ssss],[_shot,_wave],
         '''
         srsyn verb=y
         nx=%(nx)d ox=%(ox)g dx=%(dx)g
         wav=${SOURCES[1]}
         swf=${TARGETS[1]}
         ''' % par )

    # swfl(nx,ny,nw,ne)
    # rwfl(nx,ny,nw,ne)
    Flow(swfl,_ssss,
         '''
         transp plane=12 |
         transp plane=23 |
         put label5=
         ''')
    Flow(rwfl,_rrrr,
         '''
         transp plane=12 |
         transp plane=23 |
         put label5=
         ''')
    
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
# datum surface wavefields
def datum(swf1,rwf1,slow,swf0,rwf0,par):
    
    zomig.Cdtone(swf1,swf0,slow,par) #      causal
    zomig.Adtone(rwf1,rwf0,slow,par) # anti-causal

def datum3(swf1,rwf1,slow,swf0,rwf0,par):
    
    zomig.Cdtone3(swf1,swf0,slow,par) #      causal
    zomig.Adtone3(rwf1,rwf0,slow,par) # anti-causal

# ------------------------------------------------------------
# migrate
def image(imag,slow,swlf,rwfl,par):
    Flow(imag,[swlf,slow,rwfl],
         '''
         srmig3 %s
         slo=${SOURCES[1]}
         rwf=${SOURCES[2]}
         ''' % param(par))
def image_cw(imag,sslo,rslo,swlf,rwfl,par):
    Flow(imag,[swlf,sslo,rslo,rwfl],
         '''
         srmig3 %s
         slo=${SOURCES[1]}
         sls=${SOURCES[2]}
         rwf=${SOURCES[3]}
         ''' % param(par))

# model
def modelPW3(data,slow,wfld,refl,par):
    Flow(    data,[    wfld,refl,slow],
          '''
          srmod3 %s
          ref=${SOURCES[1]}
          slo=${SOURCES[2]}
          ''' % param(par))
    
def modelCW3(data,sslo,rslo,wfld,refl,par):
    Flow(    data,[         wfld,refl,sslo,rslo],
          '''
          srmod3 %s
          ref=${SOURCES[1]}
          slo=${SOURCES[2]}
          sls=${SOURCES[3]}
          ''' % param(par))

# migrate w/ CIGS
def imagePW(imag,cigs,slow,swlf,rwfl,par):
    Flow([imag,cigs],[swlf,slow,rwfl],
         '''
         srmig2 %s
         slo=${SOURCES[1]}
         rwf=${SOURCES[2]}
         cig=${TARGETS[1]}
         ''' % param(par))

def imagePW3(imag,cigs,slow,swlf,rwfl,par,custom=' '):
    Flow(   [imag,cigs],   [swlf,rwfl,slow],
         '''
         srmig3 %s
         rwf=${SOURCES[1]}
         slo=${SOURCES[2]}
         cig=${TARGETS[1]}
         %s
         ''' % (param(par),custom))

def imageCW(imag,cigs,sslo,rslo,swlf,rwfl,par):
    Flow([imag,cigs],[swlf,sslo,rslo,rwfl],
         '''
         srmig2 %s
         slo=${SOURCES[1]}
         sls=${SOURCES[2]}
         rwf=${SOURCES[3]}
         cig=${TARGETS[1]}
         ''' % param(par))

def imageCW3(imag,cigs,sslo,rslo,swlf,rwfl,par):
    Flow([imag,cigs],[swlf,sslo,rslo,rwfl],
         '''
         srmig3 %s
         slo=${SOURCES[1]}
         sls=${SOURCES[2]}
         rwf=${SOURCES[3]}
         cig=${TARGETS[1]}
         ''' % param(par))
    
# migrate (cluster call)
#def cimage(imag,slow,swlf,rwfl,par):
#    Flow(imag,[swlf,slow,rwfl],
#         '''
#         srmig %s
#         slo=${SOURCES[1]}
#         <   ${SOURCES[0]}
#         rwf=${SOURCES[2]}
#         ''' % param(par), stdin=0)

# shot-profile modeling
def model(data,slow,wfld,refl,par):
    Flow(data,[wfld,slow,refl],
          '''
          srmod3 %s
          slo=${SOURCES[1]}
          ref=${SOURCES[2]}
          ''' % param(par))

# shot-profile modeling (converted waves)
def model_cw(data,sslo,rslo,wfld,refl,par):
    Flow(data,[wfld,sslo,rslo,refl],
          '''
          srmod3 %s
          slo=${SOURCES[1]}
          sls=${SOURCES[2]}
          ref=${SOURCES[3]}
          ''' % param(par))

# ------------------------------------------------------------

# generate cluster execution commands
def script(EDIR,job,imag,cigs,slow,swfl,rwfl,par,ngroup,nshots):
    mycom = 'rm ' + EDIR + '/' + job
    os.system(mycom)

    mycom = 'cp ' + bindir + '/sfsrmig2' + ' '+ EDIR
    os.system(mycom)
    
    allk = map(lambda x: '%03d' % x,range(ngroup))
    for k in allk:
        _j = '_' + imag + '.' + k + '.rsf'
        _c = '_' + cigs + '.' + k + '.rsf'
        _s = '_' + swfl + '.' + k + '.rsf'
        _r = '_' + rwfl + '.' + k + '.rsf'
        _v = '_' + slow + '.rsf'

        mycom  = 'sfsrmig2'
        mycom = mycom + param(par)
        mycom = mycom + ' <'    + _s
        mycom = mycom + ' rwf=' + _r
        mycom = mycom + ' slo=' + _v
        mycom = mycom + ' cig=' + _c
        mycom = mycom + ' >'    + _j
        mycom = 'echo "' + mycom + '" >>' + EDIR + '/' + job
        os.system(mycom)

def execute(EDIR,JOB,ngroup,nshots,imag,cigs,slow,swfl,rwfl,par):    
    script = EDIR + '/' + JOB

    f = open(script,'w')

    allk = map(lambda x: '%03d' % x,range(ngroup))
    for k in allk:
        _j = '_' + imag + '.' + k + '.rsf'
        _c = '_' + cigs + '.' + k + '.rsf'
        _s = '_' + swfl + '.' + k + '.rsf'
        _r = '_' + rwfl + '.' + k + '.rsf'
        _v = '_' + slow + '.rsf'

        mycom  = 'sfsrmig2'
        mycom = mycom + param(par)
        mycom = mycom + ' <'    + _s
        mycom = mycom + ' rwf=' + _r
        mycom = mycom + ' slo=' + _v
        mycom = mycom + ' cig=' + _c
        mycom = mycom + ' >'    + _j
        f.write(mycom+'\n')

    f.close()

#bsub = WhereIs('bsub')
        
def run(img,cig,swf,rwf,slo,imc,par,clspar,cigpar):
    if(imc=='o'): par['misc']='itype=o'
    if(imc=='x'): par['misc']='itype=x jcx=%(jmx)d nhx=%(nhx)d hsym=y                 ' % cigpar
    if(imc=='t'): par['misc']='itype=t jcx=%(jmx)d nht=%(nht)d oht=%(oht)g dht=%(dht)g' % cigpar

    if bsub:
        cluster.sprd(swf,clspar['EDIR'],4,clspar['sgroup'],clspar['selem'],clspar['sjump'])
        cluster.sprd(rwf,clspar['EDIR'],4,clspar['sgroup'],clspar['selem'],clspar['sjump'])
        cluster.copy(slo,clspar['EDIR'])

        JOB = 'SRJOB' + '-' + slo + '-' + imc

        script(clspar['EDIR'],JOB,img,cig,slo,swf,rwf,par,
               clspar['sgroup'],
               clspar['selem'])

        cluster.launch(clspar['EDIR'],JOB,
                       clspar['sgroup'],
                       clspar['time'],
                       clspar['queue'],
                       clspar['project'])

#        Command(        clspar['EDIR']+'/'+JOB, None,
#                execute(clspar['EDIR']  ,  JOB,
#                        clspar['sgroup'],
#                        clspar['selem'],
#                        img,cig,slo,swf,rwf,par),stdout=0)
#        Command(               clspar['EDIR']+'/'+JOB,None,
#                cluster.submit(clspar['EDIR']  ,  JOB,
#                               clspar['sgroup'],
#                               clspar['time'],
#                               clspar['queue'],
#                               clspar['project']))

        # submit to cluster using 'bsub < $JOB'

        cluster.summ(img,clspar['EDIR'],clspar['sgroup'])
        cluster.summ(cig,clspar['EDIR'],clspar['sgroup'])
    else:
        imagePW(img,cig,slo,swf,rwf,par)

# ------------------------------------------------------------

# first-order scattering (slowness to image)
def s2i(dslow,dimag,swfld,rwfld,bslow,par):
    Flow(dimag,[dslow,swfld,rwfld,bslow],
         '''
         srmva adj=n %s nsc=0
         swf=${SOURCES[1]}
         rwf=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))
    
# first-order scattering (image to slowness)
def i2s(dimag,dslow,swfld,rwfld,bslow,par):
    Flow(dslow,[dimag,swfld,rwfld,bslow],
         '''
         srmva adj=y %s nsc=0
         swf=${SOURCES[1]}
         rwf=${SOURCES[2]}
         slo=${SOURCES[3]}
         ''' % param(par))

# ------------------------------------------------------------


