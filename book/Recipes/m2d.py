from rsf.proj import *
from decimal import *

# --- User set example --- #

# model

model = {
    'X':     2000,
    'Z':     2000,
    'T':     0.5,
    'SelT':  0.15,
    'dx':    20.0,
    'dz':    20.0,
    'snpintvl': 1.0,
    'size'    : 12,   # FD order
    'frqcut'  : 0.9,
    'pml'     : 30,
    'vel'     : '1800+0.000002*x1*x1',
    'den'     : '2000'
}
# source & receiver
srp = {
     'bgn'   : 0.15,     # s, time of maximum ricker
     'frq'   : 12.5,     # source domain frequence
     'srcmms'  : 'y',    # point source
     'inject': 'y',      # if y, inject; if n, Initiate conditon
     'slx'   : 1000,     # source location (x), meter
     'slz'   : 1000,     # source location (x), meter
     'gdep'  : 800       # receiver location (z), meter
     }

# ------------------------------------------------------------------

def setpar(mdl, srp):
    dx = mdl['dx']
    dz = mdl['dz']
    dt = mdl['dt']    
    
    objpar = {
        'vel' : mdl['vel'],
        'den' : mdl['den'], 
        'nx'  : mdl['X']/dx+1,
        'nz'  : mdl['Z']/dz+1,
        'SelT': mdl['SelT'],
        'nt'  : int(Decimal(str(mdl['T']))/Decimal(str(dt)))+1,
        'snt' : int(Decimal(str(mdl['SelT']))/Decimal(str(dt))/ \
                Decimal(str(mdl['snpintvl']))),
        'snpi': mdl['snpintvl'],  #snap interval
        'dt'  : dt,
        'iwdt': dt*1000,          #dt for iwave
        'dx'  : dx,
        'dz'  : dz,
        'ox'  : 0.0,
        'oz'  : 0.0,
        'ot'  : 0.0,
        # source 
        'frq'   : srp['frq'],
        'wavfrq': srp['frq']/3.0,
        'bgnp'  : srp['bgn']/dt+1,
        'slx'   : srp['slx'],
        'slz'   : srp['slz'],
        'iwslz' : -1*srp['slz'],
        'gdep'  : srp['gdep'],
        'iwgdep'  : -1*srp['gdep'],
        # fd
        'size'  : mdl['size'],
        'fdsize': mdl['size']/2,
        'frqcut': mdl['frqcut'],
        'pml'   : mdl['pml'],
        'bd'    : mdl['pml']+int((mdl['size']+1)/2)
        }
    return objpar

# ------------------------------------------------------------------------

def buildmodel(par, denname, velname, iwdenname, iwvelname, denval, velval):
    name = {
        'den' : denname,
        'vel' : velname,
        }
    iwname = {
        'den' : iwdenname,
        'vel' : iwvelname,
        }
    
    value = {
        'den' : denval,
        'vel' : velval,
        }
    
    label = {
        'den': 'Density',
        'vel': 'Velocity',
        }
    
    unit = {
        'den': 'kg/m\^3\_',
        'vel': 'm/s',
        }   
    for m in ['den','vel']:
        Flow(name[m],None,
             '''
             spike d1=%(dz)g n1=%(nz)d  o1=%(oz)g
                   d2=%(dx)g n2=%(nx)d  o2=%(ox)g
                   label1=Depth unit1=m label2=Distance unit2=m | 
             '''%par + '''
             math output="%s" 
             '''%value[m])
        pml  = name[m]+'_pml'
        Flow(pml,name[m],
             '''
             expand left=%(bd)d right=%(bd)d 
                    top=%(bd)d  bottom=%(bd)d
             '''%par) 
        Flow(iwname[m],name[m],'math output="input/1000" ')

def buildsrcp(par, srcp) :
    Flow(srcp, None,
         '''
         spike n1=%(nt)d d1=%(dt)g k1=%(bgnp)g |
         ricker1 frequency=%(frq)g | 
         scale axis=1 
         '''%par)



# ----------------------------------------------------------------------------------
def sglr2(Fwf, Frec, Fsrc, Fvel, Fden, par, prefix, suffix):
    _pf = str(prefix)
    sf_ = str(suffix)

    fft = '%sfft%s' %(_pf, sf_)
    rt  = '%srt%s'  %(_pf, sf_)
    lt  = '%slt%s'  %(_pf, sf_)
    
    Flow(fft, Fvel, 'fft1 |fft3 axis=2 pad=1')

    Flow([rt, lt], [Fvel, fft], 
         '''
         isolrsg2 seed=2010 dt=%(dt)g fft=${SOURCES[1]} left=${TARGETS[1]}
         '''%par)

    Flow([Fwf,Frec], [Fsrc, lt, rt, Fvel, Fden, fft],
             '''
             sfsglr2 verb=y             rec=${TARGETS[1]} 
                     left=${SOURCES[1]} right=${SOURCES[2]} 
                     vel=${SOURCES[3]}  den=${SOURCES[4]} 
                     fft=${SOURCES[5]}  slz=%(slz)g 
                     gdep=%(gdep)d      slx=%(slx)g
                     snapinter=1 
             '''%par)

def sglfd2(Fwf, Frec, Fsrc, Fvel, Fden, par, prefix, suffix):
    _pf = str(prefix)
    sf_ = str(suffix)

    Gx =  '%sGx%s'  %(_pf, sf_)
    sxx = '%ssxx%s' %(_pf, sf_)
    sxz = '%ssxz%s' %(_pf, sf_)
    Gz =  '%sGz%s'  %(_pf, sf_)
    szx = '%sszx%s' %(_pf, sf_)
    szz = '%sszz%s' %(_pf, sf_)
    
    Flow([Gx,sxx,sxz],Fvel,
         '''
         sfsglfdcx2_7 dt=%(dt)g eps=0.0001 npk=150 
                      size=%(size)d sx=${TARGETS[1]} sz=${TARGETS[2]}
                      wavnumcut=%(frqcut)g
         '''%par)

    Flow([Gz,szx,szz],Fvel,
         '''
         sfsglfdcz2_7 dt=%(dt)g eps=0.0001 npk=150 
                      size=%(size)d sx=${TARGETS[1]} sz=${TARGETS[2]}
                      wavnumcut=%(frqcut)g
         '''%par)

    Flow([Fwf,Frec], [Fsrc,Fvel,Fden,Gx,sxx,sxz,Gz,szx,szz],
     '''
     sfsglfd2pml rec=${TARGETS[1]} 
                 vel=${SOURCES[1]} den=${SOURCES[2]}
                 Gx=${SOURCES[3]} sxx=${SOURCES[4]} sxz=${SOURCES[5]}
                 Gz=${SOURCES[6]} szx=${SOURCES[7]} szz=${SOURCES[8]}
                 freesurface=n  verb=y
                 slx=%(slx)g slz=%(slz)g pmlsize=%(pml)d snapinter=1 
                 srcdecay=y  gdep=%(gdep)d srctrunc=0.2
     ''' %par)


def iw(Fwf, Frec, Fsrc, Fvel, Fden, par, prefix, suffix):
    _pf = str(prefix)
    sf_ = str(suffix)

    # iwave source
    srctr = '%ssrctr%s' %(_pf, sf_)
    susrc = Fsrc+'%s.su'%sf_
    Flow(srctr,Fsrc,'segyheader ')
    Flow(susrc,[Fsrc, srctr],'suwrite tfile=${SOURCES[1]} endian=0')

    # iwave par
    rcdpar = '%srcdpar%s' %(_pf, sf_)
    slc    = '%sslice%s'  %(_pf, sf_)
    Flow(rcdpar,None,'spike n1=%(nt)d n2=%(nx)d d1=%(dt)g d2=1 o2=0'%par)
    Flow(slc,rcdpar,'window n1=1')
    hkeys = dict(sx=par['slx'],
                 gx='x1*%g'%par['dx'],
                 delrt=0,
                 selev=par['iwslz'],
                 gelev=par['iwgdep'])
    
    hkeysname = {'sx'   : '%siwsx%s'    %(_pf, sf_),
                 'gx'   : '%siwgx%s'    %(_pf, sf_),
                 'delrt': '%siwdelrt%s' %(_pf, sf_),
                 'selev': '%siwselev%s' %(_pf, sf_),
                 'gelev': '%siwgelev%s' %(_pf, sf_)
                 }

    hks = hkeys.keys()
    hksn =[]      
    i = 0
    hstr = ''
    for hk in hks:
        hksn += [hkeysname[hk]]
        Flow(hkeysname[hk],slc,'math output=%s | dd type=int' % str(hkeys[hk]))
        i += 1
        hstr += '%s=${SOURCES[%d]} ' % (hk,i)

    trecord = '%strecord%s' %(_pf, sf_)
    hdr = '%shdr%s.su'      %(_pf, sf_)
    
    Flow(trecord,[rcdpar]+hksn,'segyheader ' + hstr)
    Flow(hdr,[rcdpar, trecord],'suwrite tfile=${SOURCES[1]} endian=0')

    iwpar = '%siwave%s.par' %(_pf, sf_)
    suiwrec = Frec +'.su'
    Flow(iwpar,[Fden,Fvel,hdr,'demo.par',susrc],
         '''
         /bin/cat ${SOURCES[3]} > $TARGET &&
         echo
         nt=%(nt)d dt=%(iwdt)g'''%par +
        '''
        moviestep=%g
        hdrfile  = ${SOURCES[2]}
        datafile  = %s
        velocity = ${SOURCES[1]}
        density = ${SOURCES[0]}
        source = ${SOURCES[4]} >> $TARGET
        ''' %(par['iwdt']*par['snpi'],suiwrec),stdin=0,stdout=-1)

    # --- IWAVE modeling --- #
    Flow([Frec,suiwrec],iwpar,
         '''
         asg par=$SOURCE &&
         suread < ${TARGETS[1]} read=data endian=0
         ''',stdin=0)
    



