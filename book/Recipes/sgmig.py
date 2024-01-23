try:    from rsf.cluster import *
except: from rsf.proj    import *

# ------------------------------------------------------------
def param(par):
    p  = ' '
    p  = p + ' --readwrite=y'
    if('verb' in par):
        p = p + ' verb='  +     par['verb']
    if('nrmax' in par):
        p = p + ' nrmax=' + str(par['nrmax'])
    if('dtmax' in par):
        p = p + ' dtmax=' + str(par['dtmax'])
    if('tmx' in par):
        p = p + ' tmx='   + str(par['tmx'])
    if('tmy' in par):
        p = p + ' tmy='   + str(par['tmy'])
    if('pmx' in par):
        p = p + ' pmx='   + str(par['pmx'])
    if('pmy' in par):
        p = p + ' pmy='   + str(par['pmy'])
    if('misc' in par):
        p = p + ' '      +     par['misc']
    p = p + ' '
    return p

# ------------------------------------------------------------
def freqs(par):
    f  = ' '
    if('nw' in par): f = f + ' nw=' + str(par['nw'])
    if('ow' in par): f = f + ' ow=' + str(par['ow'])
    if('dw' in par): f = f + ' dw=' + str(par['dw'])
    f = f + ' '
    return f

# ------------------------------------------------------------
def wflds(wfld,cmps,par):
    Flow(wfld,cmps,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         transp plane=13 memsize=500 |
         spray axis=2 n=1 o=0 d=1 |
         put label1=mx label2=my label3=hx label4=w
         ''' % par )

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
def datum(wfl1,slow,wfl0,par):
    Flow(wfl1,[wfl0,slow],
         '''
         camig mode=d inv=n %s
         slo=${SOURCES[1]}
         ''' % param(par))

def datum3(wfl1,slow,wfl0,par):
    Flow(wfl1,[wfl0,slow],
         '''
         camig3 mode=d inv=n %s verb=y
         slo=${SOURCES[1]}
         ''' % param(par))

# ------------------------------------------------------------
def model(data,slow,imag,par):
    Flow(data,[imag,slow],
         '''
         camig mode=m inv=y %s %s
         slo=${SOURCES[1]}
         ''' % (param(par),freqs(par)))

def model3(data,slow,imag,par):
    Flow(data,[imag,slow],
         '''
         camig3 mode=m inv=y %s %s verb=y
         slo=${SOURCES[1]}
         ''' % (param(par),freqs(par)))

# ------------------------------------------------------------    
def image(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         camig mode=m inv=n %s
         slo=${SOURCES[1]}
         ''' % param(par))

def image3(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         camig3 mode=m inv=n %s verb=y
         slo=${SOURCES[1]}
         ''' % param(par))

# ------------------------------------------------------------
def wfld3(wfld,slow,data,par):
    Flow(wfld,[data,slow],
         '''
         camig3 mode=w %s verb=y
         slo=${SOURCES[1]}
         ''' % param(par))

    
#def cimage(imag,slow,data,par):
#    Flow(imag,[data,slow],
#         '''
#         camig mode=m inv=n %s
#         slo=${SOURCES[1]}
#         <   ${SOURCES[0]}
#         ''' % param(par), stdin=0)

# ------------------------------------------------------------

# generate cluster execution commands
def script(EDIR,job,imag,slow,wfld,par,ngroup,nfreqs):
    mycom = 'rm ' + EDIR + '/' + job
    os.system(mycom)

    mycom = 'cp ' + bindir + '/sfcamig' + ' '+ EDIR
    os.system(mycom)
    
    allk = ['%03d' % x for x in range(ngroup)]
    for k in allk:
        _i = '_' + imag + '.' + k + '.rsf'
        _w = '_' + wfld + '.' + k + '.rsf'
        _v = '_' + slow + '.rsf'

        mycom  = 'sfcamig'
        mycom = mycom + param(par)
        mycom = mycom + ' <'    + _w
        mycom = mycom + ' slo=' + _v
        mycom = mycom + ' >'    + _i
        mycom = mycom + ' datapath=' + EDIR
        mycom = 'echo "' + mycom + '" >>' + EDIR + '/' + job
        os.system(mycom)

def execute(EDIR,JOB,ngroup,nfreqs,imag,slow,wfld,par):
    script = EDIR + '/' + JOB

    f = open(script,'w')

    allk = ['%03d' % x for x in range(ngroup)]
    for k in allk:
        _i = '_' + imag + '.' + k + '.rsf'
        _w = '_' + wfld + '.' + k + '.rsf'
        _v = '_' + slow + '.rsf'

        mycom  = 'sfcamig'
        mycom = mycom + param(par)
        mycom = mycom + ' <'    + _w
        mycom = mycom + ' slo=' + _v
        mycom = mycom + ' >'    + _i
        mycom = mycom + ' datapath=' + EDIR
        f.write(mycom+'\n')

    f.close()

# Esteban Diaz:
# I had a problem using this module because
# the "WhereIs()"function doesn't exist
# and neither command bsub

# bsub = WhereIs('bsub')

bsub = False

def run(img,wfl,slo,par,clspar):
    if bsub:
        cluster.sprd(wfl,clspar['EDIR'],4,clspar['fgroup'],clspar['felem'],clspar['fjump'])
        cluster.copy(slo,clspar['EDIR'])
        
        JOB = 'SSJOB' + '-' + slo
        
        script(clspar['EDIR'],JOB,img,slo,wfl,par,
               clspar['fgroup'],
               clspar['felem'])
        
        cluster.launch(clspar['EDIR'],JOB,
                       clspar['fgroup'],
                       clspar['time'],
                       clspar['queue'],
                       clspar['project'])
#
#        Command(        clspar['EDIR']+'/'+JOB, None,
#                execute(clspar['EDIR']  ,  JOB,
#                        clspar['fgroup'],
#                        clspar['felem'],
#                        img,slo,wfl,par),stdout=0)
#        Command(               clspar['EDIR']+'/'+JOB,None,
#                cluster.submit(clspar['EDIR']  ,  JOB,
#                               clspar['fgroup'],
#                               clspar['time'],
#                               clspar['queue'],
#                               clspar['project']))
#        
        # submit to cluster using 'bsub < $JOB'
        
        cluster.summ(img,clspar['EDIR'],clspar['fgroup'])
    else:
        image(img,slo,wfl,par)

        
