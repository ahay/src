from rsf.proj import *
import os


# ----------------------------------------------------------------------
# FD DATA SHOT GATHERS
# ----------------------------------------------------------------------


def makeshots(shots,ddata,tdata,par):

    # Shot positions
    Flow('tss',tdata,'dd type=float | headermath output="%(os)g + (fldr-1)*(%(ds)g)" | window' % par,local=1)
    Flow('tsi','tss','math output=input/%(ds)g' % par,local=1)

    # Offset positions
    Flow('too',tdata,'dd type=float | headermath output="(tracf-1)*(%(dr)g) - (fldr-1)*(%(ds)g)" | window' % par,local=1)
    Flow('toi','too','math output=input/%(dr)g' % par,local=1)

    # Create sraw by binning
    Flow('tos','toi tsi','cat axis=2 space=n ${SOURCES[1]} | transp | dd type=int',local=1)
    Flow('sraw',[ddata, 'tos'],'intbin head=${SOURCES[1]} xkey=0 ykey=1',local=1)
    # t,o,s

    # Prepare shots by muting and windowing
    Flow(shots,'sraw',
         '''
         put label1=Time unit1=s
         unit2=km d2=%(dh)g o2=%(oh)g label2=Offset
         unit3=km d3=%(ds)g o3=%(os)g label3=Shot |
         mutter half=false t0=%(mutt)g v0=%(mutv)g |
         window min2=%(minh)g max2=%(maxh)g
         ''' % par,local=1)
    # t,h,s


# SHOTS VISUALIZATION


def shotplot(shp):
    return '''
    window n3=1 min3=%g |
    grey pclip=100 color=I gainpanel=a wantframenum=y unit1=s label1=Time
    label2=Offset unit2=km labelsz=8
    ''' % shp

def zsplot(zp):
    return ''' 
    window n2=1 f2=%d |
    grey pclip=100 color=I gainpanel=a label1=Time label2=Position
    title="Zero offset" unit1=s unit2=km labelsz=5
    ''' % zp

def viewshots(shots,par):

    Plot('shot1','shots',shotplot(par['s1']),local=1)
    Plot('shot2','shots',shotplot(par['s2']),local=1)
    Plot('shot3','shots',shotplot(par['s3']),local=1)

    # Shot gathers
    Result('shotimg','shot1 shot2 shot3','SideBySideAniso',local=1)

    # Zero offset
    Result('zero','shots',zsplot(par['zero']),local=1)


# ----------------------------------------------------------------------
# WAVEFIELDS AND SLOWNESS FOR WAVE EQUATION MIGRATION (PSPI)
# ----------------------------------------------------------------------


def fwaveslow(shot,velocity,par):

    # Receiver shot gather wavefield FFT
    Flow('rfft',shot,
         '''
         fft1 | window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d |
         spray axis=3 n=1 o=0 d=1 label=hy | spray axis=5 n=1 o=0 d=1 label=sy
         ''' % par, local=1)
    # w, hx, hy, sx, sy

    # Source wavefield FFT
    Flow('sfft',None,
         '''
         spike k1=1 n1=%(nt)d d1=%(dt)g |
         fft1 | window squeeze=n n1=%(nw)d min1=%(ow)g j1=%(jw)d | put label1=w
         ''' % par,local=1)
    # w

    # Interpolate wavefields on surface grid
    Flow('rwav swav','rfft sfft','srsyn nx=%(nmx)d dx=%(dmx)g ox=%(omx)g wav=${SOURCES[1]} swf=${TARGETS[1]}' % par,local=1)
    # w, x, y, s

    # Transpose and setting coordinates for 3-D migration
    Flow('rtra','rwav','transp plane=12 | transp plane=23',local=1)
    Flow('stra','swav','transp plane=12 | transp plane=23',local=1)
    # x, y, w, s

    # Prepare slowness on 3-D grid
    Flow('slow',velocity,
         '''
         transp | window j1=%(jsx)d | math output=1/input |
         spray axis=2 n=1 d=1 o=0 |
         window j3=%(jsz)d squeeze=n
         ''' % par,local=1)
    # x, y, z


# WAVEFIELD FREQUENCY CONTENT VISUALIZATION


def viewfreqwaves(rwav,par):

    Plot('rwav1',rwav,
         '''
         window n4=1 f4=%(fp1)d |  math output="abs(input)" | real |
         grey title=Wavefield label1=Frequency unit1=Hz label2=Distance unit2=km
         ''' % par, local=1)

    Plot('rwav2',rwav,
         '''
         window n4=1 f4=%(fp2)d |  math output="abs(input)" | real |
         grey title=Wavefield label1=Frequency unit1=Hz label2=Distance unit2=km
         ''' % par, local=1)

    Plot('rwav3',rwav,
         '''
         window n4=1 f4=%(fp3)d |  math output="abs(input)" | real |
         grey title=Wavefield label1=Frequency unit1=Hz label2=Distance unit2=km
         ''' % par, local=1)

    Result('rwavimg','rwav1 rwav2 rwav3','SideBySideAniso',local=1)
    
    
# ----------------------------------------------------------------------
# COMMON ANGLE GATHERS
# ----------------------------------------------------------------------


# COMMON OFFSET VISUALIZATION


def viewodcig(hlist,odcimg):

    img0 = odcimg+'_0'
    img1 = odcimg+'_1'
    img2 = odcimg+'_2'
    
    Plot(img0,hlist[0],'grey title=ODCIG pclip=100 wantframenum=y label1=Depth unit1=km label2=Offset unit2=km',local=1)
    Plot(img1,hlist[1],'grey title=ODCIG pclip=100 wantframenum=y label1=Depth unit1=km label2=Offset unit2=km',local=1)
    Plot(img2,hlist[2],'grey title=ODCIG pclip=100 wantframenum=y label1=Depth unit1=km label2=Offset unit2=km',local=1)
    
    Result(odcimg,[img0, img1, img2],'SideBySideAniso',local=1)
    
    
# OFFSET VECTOR VISUALIZATION


def viewvodcig(hxzlist,vodcimg):
    
    vimg0 = vodcimg+'_0'
    vimg1 = vodcimg+'_1'
    vimg2 = vodcimg+'_2'

    bias=0.0
    
    Plot(vimg0,hxzlist[0],'grey title="Local offset vector" pclip=100 bias=%g color=j label1="hx" unit1=km label2="hz" unit2=km' % bias ,local=1)
    Plot(vimg1,hxzlist[1],'grey title="Local offset vector" pclip=100 bias=%g color=j label1="hx" unit1=km label2="hz" unit2=km' % bias ,local=1)
    Plot(vimg2,hxzlist[2],'grey title="Local offset vector" pclip=100 bias=%g color=j label1="hx" unit1=km label2="hz" unit2=km' % bias ,local=1)
    
    Result(vodcimg,[vimg0, vimg1, vimg2],'SideBySideAniso',local=1)
    
    
# SLANT STACK VISUALIZATION


def sstiplot(custom):
    return '''
    grey unit1=km unit2=  wantframenum=y label1=z title= label2="tan \F10 q\F3" %s
    ''' % (custom)

def sstdplot(custom):
    return '''
    grey unit1=km unit2=  wantframenum=y label1=hx label2="tan \F10 g\F3"  title= %s
    ''' % (custom)

def viewsst(sstlist,sstimg):

    simg0 = sstimg+'_0'
    simg1 = sstimg+'_1'
    simg2 = sstimg+'_2'
    
    Plot(simg0,sstlist[0],sstiplot(' '),local=1)
    Plot(simg1,sstlist[1],sstiplot(' '),local=1)
    Plot(simg2,sstlist[2],sstiplot(' '),local=1)
    
    Result(sstimg,[simg0, simg1, simg2],'SideBySideAniso',local=1)

def viewdsst(dsstlist,dsstimg):

    dsimg0 = dsstimg+'_0'
    dsimg1 = dsstimg+'_1'
    dsimg2 = dsstimg+'_2'
    
    Plot(dsimg0,dsstlist[0],sstdplot(' '),local=1)
    Plot(dsimg1,dsstlist[1],sstdplot(' '),local=1)
    Plot(dsimg2,dsstlist[2],sstdplot(' '),local=1)
    
    Result(dsstimg,[dsimg0, dsimg1, dsimg2],'SideBySideAniso',local=1)
    
    
# COMMON ANGLE GATHERS


def adciginc(hcig,flist,par):
    
    # Slant stack
    Flow(flist[0],hcig,'slant adj=y p0=%(op)g np=%(np)d dp=%(dp)g' % par,local=1)
    # cz, tan
    # hx, tan
    
    Flow(flist[1],flist[0],'tan2ang a0=%(oa)g na=%(na)d da=%(da)g' % par,local=1)
    # cz, incidence
    # hx, dip

    
# CAG VISUALIZATION


def cigiplot(custom):
    return '''
    grey pclip=100 min2=-40 max2=40 gainpanel=a label1=z unit1=km title=CIAG label2="\F10 q\F3" unit2="\^o\_" %s
    ''' % (custom)

def cigdplot(custom):
    return '''
    grey pclip=100 label1=hx unit1=km title=CIAG label2="\F10 g\F3" unit2="\^o\_" %s
    ''' % (custom)
    
def viewcag(aciglist,par,cigimg):

    cimg0 = cigimg+'_0'
    cimg1 = cigimg+'_1'
    cimg2 = cigimg+'_2'

    Plot(cimg0,aciglist[0],cigiplot(' '),local=1)
    Plot(cimg1,aciglist[1],cigiplot(' '),local=1)
    Plot(cimg2,aciglist[2],cigiplot(' '),local=1)
    
    Result(cigimg,[cimg0, cimg1, cimg2],'SideBySideAniso',local=1)

def viewdcag(daciglist,par,dcigimg):
    
    dcimg0 = dcigimg+'_0'
    dcimg1 = dcigimg+'_1'
    dcimg2 = dcigimg+'_2'
    
    Plot(dcimg0,daciglist[0],cigdplot(' '),local=1)
    Plot(dcimg1,daciglist[1],cigdplot(' '),local=1)
    Plot(dcimg2,daciglist[2],cigdplot(' '),local=1)
    
    Result(dcigimg,[dcimg0, dcimg1, dcimg2],'SideBySideAniso',local=1)
    
    
# ----------------------------------------------------------------------
# TAU-P GATHERS
# ----------------------------------------------------------------------


# INVERSE TAU-P NMO


def invtaupnmo(isst,loc,velocity,par):
    
    Flow('vdp1',velocity,'window n2=1 min2=%g squeeze=y' % loc,local=1)
    Flow('vsst','vdp1','put d2=%(dp)g o2=%(op)g | spray axis=2 n=%(np)d' % par,local=1)
    
    Flow('pgath',[isst, 'vsst'],
         '''
         tan2ang na=%(ntp)d da=%(dtp)g a0=%(otp)g top=y velocity=${SOURCES[1]} |
         depth2time velocity=${SOURCES[1]} nt=%(nt)d dt=%(dt)g
         ''' % par,local=1)
    
    Flow('vint','vdp1','depth2time velocity=$SOURCE nt=%(nt)d dt=%(dt)g t0=%(ot)g' % par,local=1)
    Flow('itaup2','pgath vint','itaupmo velocity=${SOURCES[1]}',local=1)

    
# SEG-Y CONVERSION


def segytaup(tpfile,par,outsegy):

    # Number of traces
    ndata = par['nxtp']*par['ntp']
    
    # Data = tau, x, p

    # tracl = trace number (1,2)
    Flow('tracl',tpfile,'window n1=1 | put n2=1 n1=%d o1=1 d1=1 | math output=x1 | dd type=int' % ndata,local=1)

    # cdp = x position (6)
    Flow('cdp',tpfile,'window n1=1 squeeze=y | put o1=1 d1=1 | math output=x1 | put n2=1 n1=%d o1=1 d1=1 | dd type=int' % ndata,local=1)

    # cdpt = p position (7)
    Flow('cdpt',tpfile,'window n1=1 squeeze=y | put o2=1 d2=1 | math output=x2 | put n2=1 n1=%d o1=1 d1=1 | dd type=int' % ndata,local=1)

    # offsp = p-offset (12)
    Flow('offsp',tpfile,'window n1=1 squeeze=y | math output=x2 | put n2=1 n1=%d o1=1 d1=1 | dd type=int' % ndata,local=1)
    
    # Header traces
    Flow('theader',[tpfile,'tracl','cdp','cdpt','offsp'],
         '''
         segyheader tracl=${SOURCES[1]} tracr=${SOURCES[1]} cdp=${SOURCES[2]} cdpt=${SOURCES[3]} offset=${SOURCES[4]} 
         ''',local=1)
    
    # SEG-Y = trace header + data + (ascii header + binary header)
    Flow(outsegy,[tpfile, 'theader'],'segywrite su=n tape=$TARGET verb=y tfile=${SOURCES[1]}', local=1)
    
