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
            rect1=10,rect2=10, # smoothing for dip estimation
            p0=0, pmin=-100,   # initial and minimum dips
            clip=3,            # clip percentile
            eps=0.1,           # regularization
            nsp=200            # number of spikes
            ):
    'Seislet transform fun'

    global nr
    
    Result(data,'grey title=Input')

    dip = data+'dip'
    Flow(dip,data,
         'dip rect1=%d rect2=%d p0=%g pmin=%g' % (rect1,rect2,p0,pmin))
    Result(dip,'grey color=j title=Slope scalebar=y')

    seis = data+'seis'
    Flow(seis,[data,dip],'seislet dip=${SOURCES[1]} eps=%g adj=y inv=y' % eps)
    Result(seis,'grey title="Seislet Transform" ')

#    sseis = data+'sseis'
#    Flow(sseis,[data,dip],
#         'seislet dip=${SOURCES[1]} eps=%g adj=y niter=100' % eps)
#    Result(sseis,'grey title="Sparse Seislet Transform" ')

    sinv = data+'sinv'
#    ssinv = data+'ssinv'

    Flow(sinv,[seis,dip],'seislet dip=${SOURCES[1]} eps=%g' % eps)
    Result(sinv,'grey title="Inverse Seislet Transform" ')

#    Flow(ssinv,[sseis,dip],'seislet dip=${SOURCES[1]} eps=%g' % eps)
#    Result(ssinv,'grey title="Inverse Seislet Transform" ')

    wvlt = data+'wvlt'

    Flow(wvlt,data,'transp | dwt | transp')
    Result(wvlt,'grey title="Wavelet Transform" ')

    for c in (1,clip,25):
        rec = '%ssrec%d' % (data,c)
        Flow(rec,[seis,dip],
             '''
             threshold pclip=%d |
             seislet dip=${SOURCES[1]} eps=%g
             ''' % (c,eps))
        Result(rec,'grey title="Inverse Seislet Transform (%d%%)" ' % c)
        wrec = '%swrec%d' % (data,c)
        Flow(wrec,wvlt,
             'threshold pclip=%d | transp | dwt adj=y inv=y | transp' % c)
        Result(wrec,'grey title="Inverse Wavelet Transform (%d%%)" ' % c)

    max=int(math.log(n2)/math.log(2))
    for m in xrange(max):
        scale = int(math.pow(2,m))
        slet = '%sslet%d' % (data,scale)
        Flow(slet,[seis,dip],
             'cut f2=%d | seislet dip=${SOURCES[1]} eps=%g' % (scale,eps))
        Result(slet,'grey title="Scale=%d" ' % scale)
        diff = '%sdiff%d' % (data,scale)
        Flow(diff,[seis,dip],
             'cut n2=%d | seislet dip=${SOURCES[1]} eps=%g' % (scale,eps))
        Result(diff,'grey title="Scale=%d" ' % scale)

    nr = n1
    k1 = string.join(map(rnd,range(nsp)),',')
    nr = n2
    k2 = string.join(map(rnd,range(nsp)),',')

    imps = data+'imps'
    Flow(imps,dip,
     '''
     spike nsp=%d k1=%s k2=%s n1=%d n2=%d |
     seislet dip=$SOURCE eps=%g
     ''' % (nsp,k1,k2,n1,n2,eps),stdin=0)
    Result(imps,'grey title=Seislets')
