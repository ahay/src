$CWPROOT/bin/suspike ntr=1 nt=551  offset=0 nspk=1 ix1=1 it1=201 dt=0.002 | $CWPROOT/bin/sushw key=delrt  a=-400 |$CWPROOT/bin/sufilter f=0.,2.5,15.0,20.0 amps=0.,1.,1.,0. | $CWPROOT/bin/sugain scale=1.e6 > $DATAPATH/wavelet.su

$CWPROOT/bin/sunull nt=1501 ntr=1204 dt=0.002 | $CWPROOT/bin/sushw key=sx a=3000 c=250 j=301 | $CWPROOT/bin/sushw key=sy a=3300 c=0 j=301 | $CWPROOT/bin/sushw key=gx a=100 b=20 j=301 | $CWPROOT/bin/sushw key=gy a=3100 c=0 j=301 | $CWPROOT/bin/sushw key=delrt a=0| $CWPROOT/bin/sushw key=selev a=-40 | $CWPROOT/bin/sushw key=gelev a=-20 > $DATAPATH/hdr3d.su




