from rsf.proj import *

def igrey(custom):
    return '''
    grey %s 
    ''' % (custom)

def wflds(swfl,rwfl,wave,shot,par):
    _wave = '_' + wave
    _shot = '_' + shot

    _ssss = '_' + swfl
    _rrrr = '_' + rwfl
    
    # _wave(nw)
    Flow(_wave,wave,
         '''
         fft1 inv=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g |
         put label1=w
         ''' % par )
    # _shot(nw,nrx,nry,nsx,nsy)
    Flow(_shot,'shot',
         '''
         window n3=%(ns)d f3=%(fs)d j3=%(js)d |
         fft1 inv=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g |
         spray axis=3 n=1 o=0 d=1 |
         spray axis=5 n=1 o=0 d=1 |
         put label1=w label2=rx label3=ry label4=sx label5=sy
         ''' % par )

    # _rrrr(nw,nx,ny,ne)
    # _ssss(nw,nx,nt,ne)
    Flow([_rrrr,_ssss],[_shot,_wave],
         '''
         srsyn verb=y
         nx=%(nx)d ox=%(ox)g dx=%(dx)g
         wav=${SOURCES[1]}
         swf=${TARGETS[1]}
         ''' % par )

    # d0s(nx,ny,nw,ne)
    # d0r(nx,ny,nw,ne)
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

def datum(swf1,rwf1,slow,swf0,rwf0,par):
    _slow = '_' + slow

    # half slowness because zomig is two-way
    Flow(_slow,slow,'math "output=input/2"') 

    Flow(swf1,[swf0,_slow],
         '''
         %(ZOM)s mode=d inv=y
         slo=${SOURCES[1]}
         ''' % par )
    Flow(rwf1,[rwf0,_slow],
         '''
         %(ZOM)s mode=d inv=n
         slo=${SOURCES[1]}
         ''' % par )


def image(imag,slow,swlf,rwfl,par):
    Flow(imag,[swlf,slow,rwfl],
         '''
         %(SPM)s
         slo=${SOURCES[1]}
         rwf=${SOURCES[2]}
         ''' % par )

    Result(imag,imag,'window | transp |'+ igrey('pclip=99'))

def cimage(imag,slow,swlf,rwfl,par):
    Flow(imag,[swlf,slow,rwfl],
         '''
         %(SPM)s
         slo=${SOURCES[1]}
         <   ${SOURCES[0]}
         rwf=${SOURCES[2]}
         ''' % par, stdin=0)
