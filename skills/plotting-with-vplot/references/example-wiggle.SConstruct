from rsf.proj import *

Flow('traces', None,
     '''
     spike n1=500 n2=15 k1=100,300 nsp=2 mag=1,0.6 |
     ricker1 frequency=25
     ''')

Result('traces',
       '''
       wiggle clip=0.8 title="Wiggle display"
              label1=Time unit1=s label2=Trace unit2=""
              poly=y zplot=1.5
       ''')

End()
