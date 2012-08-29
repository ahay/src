import os
Import('MADAGASCAR', 'TMPDATAPATH')

if MADAGASCAR:
    print 'building SEGY input data, Madagascar mode'
    from rsf.proj import *
    
    Flow('wavelet',None,
         '''
         spike n1=401 o1=-0.4 d1=0.002 k1=201 |
         trapez frequency=0.,2.5,15.0,20.0 |
         sfscale dscale=1.e6 
         ''')
    Flow('thead','wavelet','segyheader')
    Flow(TMPDATAPATH + '/wavelet.su','wavelet thead','suwrite tfile=${SOURCES[1]} endian=0')
         
    Flow('data',None,'spike n1=1501 n2=301 d1=0.002 d2=1 o2=0')
    Flow('slice','data','window n1=1')
    
    hkeys = dict(sx=3300,gx='100+20*x1',delrt=0,selev=-40,gelev=-20)
    hks = hkeys.keys()
    
    i=0
    hstr = ''
    for hk in hks:
        Flow(hk,'slice','math output=%s | dd type=int' % str(hkeys[hk]))
        i += 1
        hstr += '%s=${SOURCES[%d]} ' % (hk,i)
    
    Flow('tdata',['data']+hks,'segyheader ' + hstr) 
    Flow(TMPDATAPATH + '/hdr.su','data tdata','suwrite tfile=${SOURCES[1]} endian=0')
    
else:
    print 'building SEGY input data, standalone mode (using SU)'
    os.system('suspike ntr=1 nt=501  offset=0 nspk=1 ix1=1 it1=201 dt=0.002 | sushw key=delrt  a=-400 |sufilter f=0.,2.5,15.0,20.0 amps=0.,1.,1.,0. | sugain scale=1.e6 > ' + TMPDATAPATH + '/wavelet.su')
    os.system('sunull nt=1501 ntr=301 dt=0.002 | sushw key=sx a=3300 c=0 j=301| sushw key=gx a=100 b=20 j=301 | sushw key=delrt a=0| sushw key=selev a=-40 | sushw key=gelev a=-20 > ' + TMPDATAPATH + '/hdr.su')


