try:    from rsf.cluster import *
except: from rsf.proj    import *
from math import *

# ------------------------------------------------------------
def shot2gridw(swfl,rwfl,wave,shot,par):
    _ssss = '_' + swfl
    _rrrr = '_' + rwfl

    # _wave(nw)
    Flow(wave+'-frq',wave,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         put label1=w
         ''' % par )    
    Flow(shot+'-frq',shot,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         put label1=w 
         ''' % par )

    shot2grid(swfl,rwfl,wave+'-frq',shot+'-frq',par)

# ------------------------------------------------------------
def shot2grid(swfl,rwfl,wave,shot,par):
    _ssss = '_' + swfl
    _rrrr = '_' + rwfl
    
    Flow(shot+'-tmp',shot,
         '''
         spray axis=3 n=1 o=0 d=1 |
         spray axis=5 n=1 o=0 d=1 |
         put label1=t label2=rx label3=ry label4=sx label5=sy
         ''' % par )

    # _rrrr(nw,nx,ny,ne)
    # _ssss(nw,nx,nt,ne)
    Flow([_rrrr,_ssss],[shot+'-tmp',wave],
         '''
         shot2grid verb=y
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
def time2freq(dtime,dfreq,par):
    # input  is t-x-y-e
    # output is x-y-w-e

    Flow(dfreq,dtime,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         transp plane=12 | transp plane=23 | 
         put label1=x label2=y label3=w label4=e
         o2=%(oy)g d2=%(dy)g unit1=%(ux)s label1=%(lx)s unit2=%(uy)s label2=%(ly)s
         ''' % par)

def freq2time(dfreq,dtime,par):
    # input  is x-y-w-e
    # output is t-x-y-e

    Flow(dtime,dfreq,
         '''
         transp plane=23 | transp plane=12 |
         pad beg1=%(fw)d n1out=%(nt)d |
         fft1 inv=y opt=n |
         put label1=t label2=x label3=y label4=e
         ''' % par)


def t2f(dfreq,dtime,par):
    # input  is t-x
    # output is f-x

    Flow(dfreq,dtime,
         '''
         fft1 inv=n opt=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         put label1=f unit1=Hz
         ''' % par)

def f2t(dtime,dfreq,par):
    # input  is f-x
    # output is t-x

    Flow(dtime,dfreq,
         '''
         pad beg1=%(fw)d n1out=%(nt)d |
         fft1 inv=y opt=n |
         window j1=2 | pad n1out=%(nt)d |
         put label1=t unit1=s
         ''' % par)

         
# ------------------------------------------------------------
def delay(dou,din,delay,par):

    Flow(dou,[din,delay],
         '''
         encode del=${SOURCES[1]} verb=y
         ''' % par)

# ------------------------------------------------------------
# construct random delays
def random(delay,ns,os,ds,ne,dtmax,dtlag,par):

    Flow(delay,None,
         '''
         spike nsp=1 mag=0 n1=%d o1=%g d1=%g n2=%d |
         noise type=n |
         scale axis=123 |
         scale rscale=%g |
         add add=%g
         ''' % (ns,os,ds,ne,dtmax,dtlag) )

# ------------------------------------------------------------
# construct linear delays
def linear(delay,ns,os,ds,ne,dtmax,dtlag,par):

    dt = 2*dtmax/(ne-1)
    dd = 0.5*ns*ds
    mm = os+dd

    op = - dtmax / dd
    dp = dt/dd

    Flow(delay,None,
         '''
         math
         n1=%d o1=%g d1=%g
         n2=%d o2=0  d2=1
         output="%g+(%g+x2*%g)*(x1-%g)"
         ''' % (ns,os,ds,ne,
                dtlag,
                op,dp,
                mm
                ))
    
# ------------------------------------------------------------
def linrand(delay,ns,os,ds,ne,dtmax,dtlag,perc,par):

    dt = 2*dtmax*perc/(ne-1)
    dd = 0.5*ns*ds
    mm = os+dd

    op = - dtmax*perc / dd
    dp = dt/dd

    Flow(delay+'-lin',None,
         '''
         math
         n1=%d o1=%g d1=%g
         n2=%d o2=0  d2=1
         output="%g+(%g+x2*%g)*(x1-%g)"
         ''' % (ns,os,ds,ne,
                dtlag,
                op,dp,
                mm
                ))

    Flow(delay+'-rnd',None,
         '''
         spike nsp=1 mag=0 n1=%d o1=%g d1=%g n2=%d |
         noise type=n |
         scale axis=123 |
         scale rscale=%g |
         add add=%g
         ''' % (ns,os,ds,ne,(1-perc)*dtmax,0) )

    Flow(delay,[delay+'-lin',delay+'-rnd'],'add ${SOURCES[1]}')

# ------------------------------------------------------------
# construct harmonic delays
def harmonic(delay,ns,os,ds,ne,dtmax,dtlag,par):

    kk = 2.0*pi/(ns*ds)
    dd = 2.0*pi/ne

    Flow(delay,None,
         '''
         math
         n1=%d o1=%g d1=%g
         n2=%d o2=0  d2=1
         output="%g+%g*sin((%g)*x1-x2*%g)"         
         ''' % (ns,0,ds,ne,
                dtlag,dtmax,
                kk,dd
         ))

# ------------------------------------------------------------
def amplitude(delay,ns,os,ds,ne,dtmax,dtlag,par):

    kk = 2.0*pi/(ns*ds)
    dd = 2.0*pi/ne

    Flow(delay,None,
         '''
         math
         n1=%d o1=%g d1=%g
         n2=%d o2=0  d2=1
         output="%g+(%g*x2)*sin((%g)*x1-x2*%g)"         
         ''' % (ns,0,ds,ne,
                dtlag,dtmax/(ne-1),
                kk,dd
         ))
