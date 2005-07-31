from rsfproj import *

def param(par):
    p  = ' '
    p += ' readwrite=y'
    if(par.has_key('verb')):  p += ' verb='  +     par['verb']
    if(par.has_key('nrmax')): p += ' nrmax=' + str(par['nrmax'])
    if(par.has_key('dtmax')): p += ' dtmax=' + str(par['dtmax'])
    if(par.has_key('tmx')):   p += ' tmx='   + str(par['tmx'])
    if(par.has_key('tmy')):   p += ' tmy='   + str(par['tmy'])
    if(par.has_key('pmx')):   p += ' pmx='   + str(par['pmx'])
    if(par.has_key('pmy')):   p += ' pmy='   + str(par['pmy'])
    if(par.has_key('misc')):  p += ' '       +     par['misc']
    if(par.has_key('nsc')):   p += ' nsc='   + str(par['nsc'])
    p += ' '
    return p

def freqs(par):
    f  = ' '
    if(par.has_key('nw')): f += ' nw=' + str(par['nw'])
    if(par.has_key('ow')): f += ' ow=' + str(par['ow'])
    if(par.has_key('dw')): f += ' dw=' + str(par['dw'])
    if(par.has_key('jw')): f += ' jw=' + str(par['jw'])
    f += ' '
    return f

def wflds(wfld,data,par):
    Flow(wfld,data,
         '''
         fft1 inv=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         transp plane=12 memsize=500 |
         transp plane=23 memsize=500 |
         put label1=mx label2=my label4=w
         ''' % par )

# zero-offset modeling
def model(data,slow,imag,par):
    Flow(data,[imag,slow],
        '''
        zomig mode=m inv=y %s %s
        slo=${SOURCES[1]}
        ''' % (param(par),freqs(par)))

# zero-offset migration
def image(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         zomig mode=m inv=n %s
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
def Cdttwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=y twoway=y %s
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
def Adttwo(wfld,data,slow,par):
    Flow(wfld,[data,slow],
         '''
         zomig mode=d inv=n causal=n twoway=y %s
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



