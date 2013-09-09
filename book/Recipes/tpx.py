from rsf.proj import *

def FPX(fpx,data,
        np,               # number of slopes
        nw,               # number of frequencies
        p0=-1,            # first slope
        dp=None,          # slope increment
        v0=0,             # velocity continuation
        ):

    if not dp:
        dp=-2.0*p0/(np-1)

    # TX -> FX
    fx = 'fx-'+data
    if (v0 > 0):
        Flow(fx,data,
             '''
             fft1 | window n1=%d | fft3 axis=2 | 
             vczo2 v0=0 nv=1 dv=%g | 
             window | fft3 axis=2 inv=y
             ''' % (nw,v0))
    else:
        Flow(fx,data,'fft1 | window n1=%d' % nw)

    # FX -> XPF
    xpf = 'xpf-'+data
    basis = 'basis-'+data
    Flow([xpf,basis],fx,
         '''
         transp |
         cltft basis=${TARGETS[1]} dip=y 
         p0=%g dp=%g np=%d 
         rect=3 niter=1000 verb=n
         ''' % (p0,dp,np),split=[1,nw],
         reduce='cat axis=3')
    
    Flow(fpx,[xpf,basis],
         'mul ${SOURCES[1]} | transp plane=13',
         split=[2,np])

def TPX(tpx,data,
        nt,               # number of time samples
        np,               # number of slopes
        nw=0,             # number of frequencies
        p0=-1,            # first slope
        dp=None,          # slope increment
        ):

    fpx = 'fpx-'+data

    nt2=nt
    if nt2%2:
        nt2 += 1
    nw0=nt2/2+1
    if not nw:
        nw = nw0
        
    FPX(fpx,data,np,nw,p0,dp)

    Flow(tpx,fpx,
         '''
         pad n1=%d | fft1 inv=y
         ''' % nw0,split=[3,'omp'])
    
