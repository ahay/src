import os 
import shutil
Import('MADAGASCAR TMPDATAPATH ROOT')

if MADAGASCAR:
    from rsf.proj import *
    os.system(ROOT + 'sfspike n1=401 o1=-0.4 d1=0.002 k1=201 | ' + ROOT + 'sftrapez frequency=0.,2.5,15.0,20.0 | ' + ROOT + 'sfscale dscale=1.e6 > wavelet.rsf')
    os.system('< wavelet.rsf ' + ROOT + 'sfsegyheader > thead.rsf')
    os.system('< wavelet.rsf ' + ROOT + 'sfsuwrite tfile=thead.rsf endian=0 > ' + TMPDATAPATH + '/wavelet.su')
    os.system(ROOT + 'sfspike n1=1501 n2=301 d1=0.002 d2=1 o2=0 > data.rsf')
    os.system('< data.rsf ' + ROOT + 'sfwindow n1=1 > slice.rsf')
    hkeys = dict(sx=3300,gx='100+20*x1',delrt=0,selev=-40,gelev=-20)
    hks = hkeys.keys()
    i=0
    hstr = ''
    for hk in hks:
        cmda = 'sfmath output=' + str(hkeys[hk])
        cmd = ROOT + cmda +' | ' + ROOT + 'sfdd type=int'
        os.system('< slice.rsf ' + cmd + ' > ' + hk + '.rsf')
        i += 1
        hstr += hk + '=' + hk + '.rsf '
#'%s=${SOURCES[%d]} ' % (hk,i)
    os.system('< data.rsf ' + ROOT + 'sfsegyheader ' + hstr + ' > tdata.rsf')    
#    Flow('tdata',['data']+hks,'segyheader ' + hstr) 
    os.system('< data.rsf ' + ROOT + 'sfsuwrite tfile=tdata.rsf endian=0 > ' + TMPDATAPATH + '/hdr.su')
#    Flow('hdr.su','data tdata','suwrite tfile=${SOURCES[1]} endian=0')
    os.system(ROOT + 'sfrm wavelet.rsf; ' + ROOT + 'sfrm thead.rsf; ' + ROOT + 'sfrm data.rsf; ' + ROOT + 'sfrm slice.rsf; ' + ROOT + 'sfrm tdata.rsf; ')
    for hk in hks:
        os.system(ROOT + 'sfrm ' + hk + '.rsf')
else:
    print 'building SEGY input data, standalone mode (using SU)'
    os.system('suspike ntr=1 nt=501  offset=0 nspk=1 ix1=1 it1=201 dt=0.002 | sushw key=delrt  a=-400 |sufilter f=0.,2.5,15.0,20.0 amps=0.,1.,1.,0. | sugain scale=1.e6 > ' + TMPDATAPATH + '/wavelet.su')
    os.system('sunull nt=1501 ntr=301 dt=0.002 | sushw key=sx a=3300 c=0 j=301| sushw key=gx a=100 b=20 j=301 | sushw key=delrt a=0| sushw key=selev a=-40 | sushw key=gelev a=-20 > ' + TMPDATAPATH + '/hdr.su')



