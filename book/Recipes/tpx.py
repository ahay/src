from rsf.proj import *

def TPX(tpx,data,
        nt,               # number of time samples
        np,               # number of slopes
        nw=0,             # number of frequencies
        p0=-1,            # first slope
        ):

    dp=-2.0*p0/(np-1)

    nt2=nt
    if nt2%2:
        nt2 += 1
    nw0=nt2/2+1
    if not nw:
        nw = nw0

    # TX -> FX
    fx = 'fx-'+data
    Flow(fx,data,'fft1 | window n1=%d' % nw) 

    # FX -> XPF
    xpf = 'xpf-'+data
    basis = 'basis-'+data
    Flow([xpf,basis],fx,
         '''
         transp |
         cltft basis=${TARGETS[1]} dip=y p0=%g dp=%g np=%d rect=3 niter=1000 verb=n
         ''' % (p0,dp,np),split=[1,nw],reduce='cat axis=3')
    
    Flow(tpx,[xpf,basis],
         '''
         mul ${SOURCES[1]} |
         transp plane=13   |
         pad n1=%d |
         fft1 inv=y
         ''' % nw0,split=[2,np])
