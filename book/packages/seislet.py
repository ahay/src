from rsfproj import *
import random, string, math

random.seed(2005)
nr = 0

def rnd(x):
    global nr
    r = str(random.randint(1,nr))
    return r

def seislet(data,              # data name
            n1,n2,             # data dimensions
            o2=0,d2=1,         # lateral scale
            rect1=10,rect2=10, # smoothing for dip estimation
            p0=0, pmin=-100,   # initial and minimum dips
            clip=3,            # clip percentile
            eps=0.1,           # regularization
            nsp=200            # number of spikes
            ):
    'Seislet transform fun'

    global nr
    
    Result(data,'grey  title=Input')

    dip = data+'dip'
    Flow(dip,data,
         'dip rect1=%d rect2=%d p0=%g pmin=%g' % (rect1,rect2,p0,pmin))
    Result(dip,'grey  color=j title=Slope scalebar=y')

    seis = data+'seis'
    Flow(seis,[data,dip],'seislet dip=${SOURCES[1]} eps=%g adj=y inv=y' % eps)
    Result(seis,'put o2=0 d2=1 | grey  title="Seislet Transform" label2=Scale unit2=')

#    sseis = data+'sseis'
#    Flow(sseis,[data,dip],
#         'seislet dip=${SOURCES[1]} eps=%g adj=y niter=100' % eps)
#    Result(sseis,'grey  title="Sparse Seislet Transform" ')

    sinv = data+'sinv'
#    ssinv = data+'ssinv'

    Flow(sinv,[seis,dip],'seislet dip=${SOURCES[1]} eps=%g' % eps)
    Result(sinv,'grey  title="Inverse Seislet Transform" ')

#    Flow(ssinv,[sseis,dip],'seislet dip=${SOURCES[1]} eps=%g' % eps)
#    Result(ssinv,'grey  title="Inverse Seislet Transform" ')

    wvlt = data+'wvlt'

    Flow(wvlt,data,'transp | dwt | transp')
    Result(wvlt,'put o2=0 d2=1 | grey  title="Wavelet Transform" label2=Scale unit2=')

    for c in (1,clip,25):
        rec = '%ssrec%d' % (data,c)
        Flow(rec,[seis,dip],
             '''
             threshold pclip=%d |
             seislet dip=${SOURCES[1]} eps=%g
             ''' % (c,eps))
        Result(rec,'grey  title="Inverse Seislet Transform (%d%%)" ' % c)
        wrec = '%swrec%d' % (data,c)
        Flow(wrec,wvlt,
             'threshold pclip=%d | transp | dwt adj=y inv=y | transp' % c)
        Result(wrec,'grey  title="Inverse Wavelet Transform (%d%%)" ' % c)

    max=int(math.log(n2)/math.log(2))
    for m in xrange(max):
        scale = int(math.pow(2,m))
        slet = '%sslet%d' % (data,scale)
        Flow(slet,[seis,dip],
             'cut f2=%d | seislet dip=${SOURCES[1]} eps=%g' % (scale,eps))
        Result(slet,'grey  title="Scale=%d" ' % scale)
        diff = '%sdiff%d' % (data,scale)
        Flow(diff,[seis,dip],
             'cut n2=%d | seislet dip=${SOURCES[1]} eps=%g' % (scale,eps))
        Result(diff,'grey  title="Scale=%d" ' % scale)

    nr = n1
    k1 = string.join(map(rnd,range(nsp)),',')
    nr = n2
    k2 = string.join(map(rnd,range(nsp)),',')

    imps = data+'imps'
    Flow(imps,dip,
     '''
     spike nsp=%d k1=%s k2=%s n1=%d n2=%d o2=%g d2=%g |
     seislet dip=$SOURCE eps=%g
     ''' % (nsp,k1,k2,n1,n2,o2,d2,eps),stdin=0)
    Result(imps,'grey  title=Seislets')

    
    impw = data+'impw'
    Flow(impw,dip,
     '''
     spike nsp=%d k1=%s k2=%s n1=%d n2=%d o2=%g d2=%g |
     transp | dwt eps=%g adj=y inv=y | transp
     ''' % (nsp,k1,k2,n1,n2,o2,d2,eps),stdin=0)
    Result(impw,'grey  title=Wavelets')

def diplet(data,              # data name
           n1,n2,             # data dimensions
           o2=0,d2=1,         # lateral scale
           pmin=-2.0,         # minimum dip
           pmax=2.0,          # maximum dip
           np=101,            # number of dips
           clip=3,            # clip percentile
           eps=0.1,           # regularization
           nsp=200            # number of spikes
           ):
    'Seislet frame fun'
    
    global nr
    
    dips = data+'dips'
    Flow(dips,data,'spray axis=3 n=%d o=%g d=%g | math output=x3' % (np,pmin,(pmax-pmin)/(np-1)))
    
    dipl = data+'dipl'
    Flow(dipl,[data,dips],'diplet dips=${SOURCES[1]} eps=%g adj=y inv=y' % eps)

    for c in (1,clip,25):
        rec = '%sdrec%d' % (data,c)
        Flow(rec,[dipl,dips],
             '''
             threshold pclip=%d |
             diplet dips=${SOURCES[1]} eps=%g
             ''' % (c,eps))
        Result(rec,'grey  title="Inverse Seislet Frame (%d%%)" ' % c)

    nr = n1
    k1 = string.join(map(rnd,range(nsp)),',')
    nr = n2
    k2 = string.join(map(rnd,range(nsp)),',')
    nr = np
    k3 = string.join(map(rnd,range(nsp)),',')

    imps = data+'dimps'
    Flow(imps,dips,
         '''
         spike nsp=%d k1=%s k2=%s k3=%s n1=%d n2=%d n3=%d o2=%g d2=%g |
         diplet dips=$SOURCE eps=%g
         ''' % (nsp,k1,k2,k3,n1,n2,np,o2,d2,eps),stdin=0)
    Result(imps,'grey  title="Seislet Frame Members" ')

   
