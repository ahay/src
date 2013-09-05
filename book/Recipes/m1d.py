from rsf.proj import *
from decimal import *


# --- User set --- # 

# model
model = {
    'X' :   2000,    # meter
    'dx':   10.0,
    'dt':   0.001,
    'SelT': 0.25,    # selected time for snapshot show 
    'snpintvl': 1.0, # nterval of snapshot output
    'size'    : 8,   # FD order
    'frqcut'  : 1.0,
    'pml'     : 240,
    }

# source & receiver
srp = {
     'bgn'   : 0.1,    # s, time of maximum ricker
     'frq'   : 10.0,   # source domain frequence
     'srcmms'  : 'n',  # MMS
     'inject': 'n',    # if y, inject; if n, Initiate conditon
     'slx'   : 1000.0, # source location (x), meter
     'gdep'  : 800     # receiver location (z), meter
     }

# ------------------------------------------------------------------------------
def mgraph(fin, title):
    Result(fin,
           '''
           put label1="Depth" unit1="m" |
           transp plane=23              |
           graph screenratio=0.5 title="%s"
           '''%str(title))

# ------------------------------------------------------------------------------
def setpar(mdl, srp):
    dx = mdl['dx']
    dt = mdl['dt']
    
    objpar = {
        'vel' : mdl['vel'],
        'dvel': mdl['dvel'],
        'den' : mdl['den'],
        'nx'  : mdl['X']/dx+1,
        'SelT': mdl['SelT'],
        'nt'  : int(Decimal(str(mdl['T']))/Decimal(str(dt)))+1,
        'snt' : int(Decimal(str(mdl['SelT']))/Decimal(str(dt))/ \
                Decimal(str(mdl['snpintvl']))),
        'snpi': mdl['snpintvl'],  #snap interval
        'dt'  : dt,
        'iwdt': dt*1000,          #dt for iwave
        'dx'  : dx,
        'dxhf': 0.5*dx,
        'ox'  : 0.0,
        'ot'  : 0.0,
        # source 
        'frq'   : srp['frq'],
        'wavfrq': srp['frq']/3.0,
        'bgnp'  : srp['bgn']/dt+1,
        'slx'   : srp['slx'],
        'spx'   : srp['slx']/dx+1,
        'gdep'  : srp['gdep'],
        'gp'  : int(srp['gdep']/dx+0.5),
        'srcmms': srp['srcmms'],   # MMS
        'inject': srp['inject'],   # if y, inject; if n, Initiate conditon
        # fd
        'size'  : mdl['size'],
        'fdsize': mdl['size']/2,
        'frqcut': mdl['frqcut'],
        'pml'   : mdl['pml'],
        'bd'    : mdl['pml']+int((mdl['size']+1)/2)
    }
    return objpar

def buildmodel(par, denname, velname, dvelname, denval, velval, dvelval):
    name = {
        'den' : denname,
        'vel' : velname,
        'dvel': dvelname
        }
    
    value = {
        'den' : denval,
        'vel' : velval,
        'dvel': dvelval
        }
    
    label = {
        'den': 'Density',
        'vel': 'Velocity',
        'dvel': 'Velocity'
        }
    
    unit = {
        'den': 'lg/m\^3\_',
        'vel': 'm/s',
        'dvel': 'm/s'
        }

    for m in ['den','vel','dvel']:
        Flow(name[m],None,
             '''
             spike d1=%(dx)g n1=%(nx)d  o1=0.0
                   label1=Depth unit1=m |                   
             '''%par + '''
             math output="%s" 
             '''%value[m])
        pml  = name[m]+'_pml'
        pmlt = name[m]+'_pmlt'
        pmlb = name[m]+'_pmlb'
        Flow(pmlt, name[m], 'window n1=1 f1= 0  |spray axis=1 n=%(bd)d' %par)
        Flow(pmlb, name[m], 'window n1=1 f1=-1  |spray axis=1 n=%(bd)d' %par)
        Flow(pml,[pmlt, name[m], pmlb],'cat ${SOURCES[1]} ${SOURCES[2]} axis=1')
    for m in ['den','vel','dvel']:
        Flow(name[m]+'hf',None,
             '''
             spike d1=%(dx)g n1=%(nx)d  o1=%(dxhf)g
                   label1=Depth unit1=m |                   
             '''%par + '''
             math output="%s" 
             '''%value[m])


def buildic(par, ic):
    Flow(ic,None,
         '''
         spike n1=%(nx)d d1=%(dx)g k1=%(spx)d| 
         ricker1 frequency=%(wavfrq)g | 
         scale axis=1   |
         put lable1="Depth" unit1="m" label2="Amplitude" unit2=""
         '''%par)


def buildsrcp(par, srcp):
    Flow(srcp, None,
     '''
      spike n1=%(nt)d d1=%(dt)g k1=%(bgnp)g |
      ricker1 frequency=%(frq)g | 
      scale axis=1 |math output=input*400
     '''%par)

def buildsrcd(par, srcd, prefix, subfix):
    _pf    = str(prefix) 
    sf_    = str(subfix)
    spike  = '%sspike%s'   %(_pf, sf_)
    ricker = '%sricker%s'  %(_pf, sf_)
    
    Flow(spike,None, 
         '''
         spike n1=%(nx)d n2=%(nt)d d1=%(dx)g d2=%(dt)g 
         k1=%(spx)d k2=1 
         '''%par)
    
    Flow(ricker,None,
         '''
         spike n1=%(nt)d d1=%(dt)g k1=%(bgnp)g | 
         ricker1 frequency=%(frq)g |scale axis=1 
         '''%par)

    Flow(srcd,[spike,ricker],
         '''
         convft other=${SOURCES[1]} axis=2 |
         window n2=%(nt)d | math output=input*400
         '''%par)

def buildmms(par, mms, psrc, vsrc, pint, vint, vel, dvel, den, velhf, dvelhf, denhf):
    #beta = 2*3.14159265359*par['frq']
    alpha = 2*3.1415926*par['frq']/4.0
    alpha = alpha*alpha
    Flow([mms, psrc, vsrc, pint, vint], [vel, dvel, den, velhf, dvelhf, denhf],
         '''
         sfmms1dexp nt=%d dt=%g slx=%g  alpha=%g 
                    dvel=${SOURCES[1]} den=${SOURCES[2]}
                    presrc=${TARGETS[1]} velsrc=${TARGETS[2]}
                    preinit=${TARGETS[3]} velinit=${TARGETS[4]} 
                    velhf=${SOURCES[3]} dvelhf=${SOURCES[4]} denhf=${SOURCES[5]}|
          put label1="Depth" unit1="km" label2="Time" unit2="s"
          '''%(par['nt'],par['dt'],par['slx'],alpha))
    

# ------------------------------------------------------------------------------

def lrmodel(fwf, frec, src, ic, vel, den, mmsfiles, par, prefix, suffix):
    _pf = str(prefix)
    sf_ = str(suffix)
    
    fft = '%sfft%s' %(_pf, sf_)
    rt  = '%srt%s'  %(_pf, sf_)
    lt  = '%slt%s'  %(_pf, sf_)
        
    Flow(fft, vel, 'fft1')
    Flow([rt, lt], [vel, fft], 
         '''
         isolrsg1 seed=2010 dt=%(dt)g fft=${SOURCES[1]} left=${TARGETS[1]}
         '''%par)
    if (mmsfiles == {}):
        Flow([fwf,frec], [src, lt, rt, vel, den, fft, ic],
             '''
             sfsglr1 verb=y             rec=${TARGETS[1]} 
                   left=${SOURCES[1]} right=${SOURCES[2]} 
                   vel=${SOURCES[3]}  den=${SOURCES[4]} 
                   fft=${SOURCES[5]}  ic=${SOURCES[6]}
                   gdep=%(gdep)g      slx=%(slx)g
                   inject=%(inject)s  srcmms=%(srcmms)s 
             '''%par)
    else :
        psrc = mmsfiles['presrc']
        vsrc = mmsfiles['velsrc']
        pint = mmsfiles['preinit']
        vint = mmsfiles['velinit']
        Flow([fwf,frec], [src, lt, rt, vel, den, fft, ic,
                          psrc, vsrc, pint, vint],
             '''
             sfsglr1 verb=y 
                   rec=${TARGETS[1]} 
                   left=${SOURCES[1]} right=${SOURCES[2]} 
                   vel=${SOURCES[3]}  den=${SOURCES[4]} 
                   fft=${SOURCES[5]}  ic=${SOURCES[6]}
                   presrc=${SOURCES[7]} velsrc=${SOURCES[8]}
                   preinit=${SOURCES[9]} velinit=${SOURCES[10]}
                   gdep=%(gdep)g      slx=%(slx)g
                   inject=%(inject)s  srcmms=%(srcmms)s 
             '''%par)
    
def lfdmodel(fwf, frec, src, ic, vel, den, mmsfiles, par, prefix, suffix):
    _pf = str(prefix)
    sf_ = str(suffix)
    G =  '%sG%s'  %(_pf, sf_)
    sx = '%ssx%s' %(_pf, sf_)

    Flow([G,sx],vel,
         '''
         sfsglfdc1 dt=%(dt)g eps=0.00001 npk=20 seed=2012
                   sx=${TARGETS[1]} size=%(size)d wavnumcut=%(frqcut)g
         ''' %par)
    
    if mmsfiles == {}:
        Flow([fwf, frec], [src, ic, vel, den, G, sx],    
             '''
             sfsglfd1pml rec=${TARGETS[1]} ic=${SOURCES[1]}
                         vel=${SOURCES[2]} den=${SOURCES[3]} 
                         G=${SOURCES[4]}   sx=${SOURCES[5]}
                         pmld0=20             
                         gdep=%(gdep)g slx=%(slx)g pmlsize=%(pml)d
                         inject=%(inject)s srcmms=%(srcmms)s
                         verb=y snapinter=1 
             ''' %par)
  
    else:
        psrc = mmsfiles['presrc']
        vsrc = mmsfiles['velsrc']
        pint = mmsfiles['preinit']
        vint = mmsfiles['velinit']
        Flow([fwf, frec], [src, ic, vel, den, G, sx,
                           psrc, vsrc, pint, vint],    
             '''
             sfsglfd1pml rec=${TARGETS[1]} ic=${SOURCES[1]}
                         vel=${SOURCES[2]} den=${SOURCES[3]} 
                         G=${SOURCES[4]}   sx=${SOURCES[5]}
                         presrc=${SOURCES[6]} velsrc=${SOURCES[7]}
                         preinit=${SOURCES[8]} velinit=${SOURCES[9]}
                         pmld0=20             
                         gdep=%(gdep)g slx=%(slx)g pmlsize=%(pml)d
                         inject=%(inject)s srcmms=%(srcmms)s
                         verb=y snapinter=1 
             ''' %par)

def fdmodel(fwf, frec, src, ic, vel, den, mmsfiles, par):
    if (mmsfiles == {}):
        Flow([fwf, frec], [src, ic, vel, den],    
             '''
             sfsgfd1 ic=${SOURCES[1]}
                     vel=${SOURCES[2]} den=${SOURCES[3]} rec=${TARGETS[1]}
                     pmld0=20 size=%(fdsize)d       
                     gdep=%(gdep)g slx=%(slx)g pmlsize=%(pml)d
                     inject=%(inject)s 
                     verb=y snapinter=1 
             ''' %par )
    else :
        psrc = mmsfiles['presrc']
        vsrc = mmsfiles['velsrc']
        pint = mmsfiles['preinit']
        vint = mmsfiles['velinit']
        Flow([fwf, frec], [src, ic, vel, den, psrc, vsrc, pint, vint],    
             '''
             sfsgfd1 ic=${SOURCES[1]}
                     vel=${SOURCES[2]} den=${SOURCES[3]} rec=${TARGETS[1]}
                     presrc=${SOURCES[4]} velsrc=${SOURCES[5]}
                     preinit=${SOURCES[6]} velinit=${SOURCES[7]}
                     pmld0=20 size=%(fdsize)d       
                     gdep=%(gdep)g slx=%(slx)g pmlsize=%(pml)d
                     inject=%(inject)s srcmms=%(srcmms)s
                     verb=y snapinter=1 
             ''' %par )
    
# ------------------------------------------------------------------------------

def analyticslt(fout, par, vel, prefix, subfix):
    _pf = str(prefix)
    sf_ = str(subfix)
    
    spx = par['spx']
    selt= par['SelT']
    dx  = par['dx']
    
    leftp  = spx - round(vel*selt/dx)
    rightp = spx + round(vel*selt/dx)
    
    left = '%sleft%s' %(_pf, sf_)
    right= '%sright%s'%(_pf, sf_)
    
    for fi in [left, right]:
        p = (leftp, rightp)[fi==right]
        Flow(fi,None,
             '''
             spike n1=%d d1=%g k1=%d| 
             ricker1 frequency=%g | math output="input"
             '''%(par['nx'],par['dx'],p,par['wavfrq']))

    Flow(fout,[left,right],
         '''
         math t=${SOURCES[1]} output="input+t" | 
         scale axis=2 | scale rscale=0.5 |
         put label1="Distance" unit1="km"
         ''')
    
  






