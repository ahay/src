try:    from rsf.cluster import *
except: from rsf.proj    import *

def param(par):
    p  = ' '
    p = p + ' --readwrite=y'
    if(par.has_key('verb')):
        p = p + ' verb='  +     par['verb']
    if(par.has_key('incore')):
        p = p + ' incore='+     par['incore']
    if(par.has_key('nrmax')):
        p = p + ' nrmax=' + str(par['nrmax'])
    if(par.has_key('dtmax')):
        p = p + ' dtmax=' + str(par['dtmax'])
    if(par.has_key('eps')):
        p = p + ' eps='   + str(par['eps'])
    if(par.has_key('tmx')):
        p = p + ' tmx='   + str(par['tmx'])
    if(par.has_key('tmy')):
        p = p + ' tmy='   + str(par['tmy'])
    if(par.has_key('pmx')):
        p = p + ' pmx='   + str(par['pmx'])
    if(par.has_key('pmy')):
        p = p + ' pmy='   + str(par['pmy'])
    if(par.has_key('ompnth')):
        p = p + ' ompnth='  + str(par['ompnth'])
    if(par.has_key('misc')):
        p = p + ' '       +     par['misc']
    if(par.has_key('nsc')):
        p = p + ' nsc='   + str(par['nsc'])
    p = p + ' '
    return p

def freqs(par):
    f  = ' '
    if(par.has_key('nw')): f = f + ' nw=' + str(par['nw'])
    if(par.has_key('ow')): f = f + ' ow=' + str(par['ow'])
    if(par.has_key('dw')): f = f + ' dw=' + str(par['dw'])
    if(par.has_key('jw')): f = f + ' jw=' + str(par['jw'])
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

def wflds(wfld,data,par):
    Flow(wfld,data,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         transp plane=12 memsize=500 |
         transp plane=23 memsize=500 |
         put label1=mx label2=my label3=w
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

#def slow(slow,velo,par):
#    Flow(slow,velo,
#         '''
#         math output="1/input" |
#         transp plane=12 memsize=500|
#         transp plane=23 memsize=500|
#         put label1=mx label2=my label3=z
#         ''')

# ------------------------------------------------------------

# zero-offset modeling
def model(data,slow,imag,par):
    Flow(data,[imag,slow],
         '''
         zowei mode=m inv=y %s %s
         slo=${SOURCES[1]}
         ''' % (param(par),freqs(par)))
    
# zero-offset migration
def image(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         zowei mode=m inv=n %s
         slo=${SOURCES[1]}
         ''' % param(par))

# zero-offset modeling
def model3(data,slow,imag,par):
    Flow(  data,    [imag,slow],
         '''
         zomig3 mode=m inv=y %s %s
         slo=${SOURCES[1]}
         ''' % (param(par),freqs(par)))
# zero-offset migration
def image3(imag,slow,data,par):
    Flow(  imag,    [data,slow],
         '''
         zomig3 mode=m inv=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
    
# ------------------------------------------------------------

#      causal datuming (forward in time, causal=y)
# (useful for datuming source wavefields)
def Cdtone(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=y twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Cdtone3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=d inv=n causal=y twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
    
def Cdttwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=y twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Cdttwo3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=d inv=n causal=y twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))
    
# anti-causal datuming (backward in time, causal=n)
# (useful for datuming receiver wavefields)
def Adtone(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=n twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Adtone3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=d inv=n causal=n twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
    
def Adttwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=n twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Adttwo3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=d inv=n causal=n twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))
# ------------------------------------------------------------

# causal wavefields by wavefield extrapolation
def Cwfone(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=w inv=n causal=y twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Cwftwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=w inv=n causal=y twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# anti-causal wavefields by wavefield extrapolation
def Awfone(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=w inv=n causal=n twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Awftwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=w inv=n causal=n twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# ------------------------------------------------------------

# causal wavefields by wavefield extrapolation
def Cwfone3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=w inv=n causal=y twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Cwftwo3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=w inv=n causal=y twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# anti-causal wavefields by wavefield extrapolation
def Awfone3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=w inv=n causal=n twoway=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
def Awftwo3(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig3 mode=w inv=n causal=n twoway=y %s
         slo=${SOURCES[1]}
         ''' % param(par))

# ------------------------------------------------------------

# first-order scattering (slowness to image)
def s2i(dslow,dimag,bwfld,bslow,par):
    Flow(dimag,[dslow,bwfld,bslow],
         '''
         zomva inv=n %s
         wfl=${SOURCES[1]}
         slo=${SOURCES[2]}
         ''' % param(par))

# first-order scattering (image to slowness)
def i2s(dimag,dslow,bwfld,bslow,par):
    Flow(dslow,[dimag,bwfld,bslow],
         '''
         zomva inv=y %s
         wfl=${SOURCES[1]}
         slo=${SOURCES[2]}
         ''' % param(par))

def s2i3(dslow,dimag,bwfld,bslow,par):
    Flow(dimag,[dslow,bwfld,bslow],
         '''
         zomva3 inv=n %s
         wfl=${SOURCES[1]}
         slo=${SOURCES[2]}
         ''' % param(par))

def i2s3(dimag,dslow,bwfld,bslow,par):
    Flow(dslow,[dimag,bwfld,bslow],
         '''
         zomva3 inv=y %s
         wfl=${SOURCES[1]}
         slo=${SOURCES[2]}
         ''' % param(par))

# ------------------------------------------------------------
# simulate shot-record migration
#def wem(imag,sdat,rdat,velo,custom,par):

#    sfrq = imag + sdat + '_f'
#    rfrq = imag + rdat + '_f'

#    swfl = imag + sdat + '_w'
#    rwfl = imag + rdat + '_w'

#    slow = imag + velo + '_s'
    
#    Flow(sfrq,sdat,
#         '''
#         transp |
#         fft1 inv=n opt=n |
#         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
#         transp plane=12 |
#         transp plane=23
#         ''' % par )

#    Flow(rfrq,rdat,
#         '''
#         transp |
#         fft1 inv=n opt=n |
#         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
#         transp plane=12 |
#         transp plane=23
#         ''' % par )

#    Flow(slow,velo,
#         '''
#         math output=1/input |
#         transp plane=12 |
#         transp plane=23
#         ''')

#    Cwfone(swfl,sfrq,slow,par)
#    Awfone(rwfl,rfrq,slow,par)

#    Flow(imag,[swfl,rwfl],
#         '''
#         math s=${SOURCES[0]} r=${SOURCES[1]} output="r*conj(s)" |
#         window |
#         stack axis=3 |
#         real |
#         transp
#         ''',stdin=0)

# ------------------------------------------------------------
# zero-offset migration
def zom(imag,data,velo,par):

    slow = imag + velo + '_s'
    Flow(slow,velo,
         '''
         math output=1/input |
         transp plane=12 |
         transp plane=23
         ''')

    freq = imag + data + '_f'
    Flow(freq,data,
         '''
         transp |
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         transp plane=12 |
         transp plane=23
         ''' % par )
    
    Flow(imag,[freq,slow],
         '''
         zomig3 mode=m inv=n %s
         slo=${SOURCES[1]}
         ''' % param(par))
