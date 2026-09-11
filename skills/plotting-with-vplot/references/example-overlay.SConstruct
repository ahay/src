from rsf.proj import *

# Data: a 2D gather.
Flow('cmp', None,
     '''
     spike n1=400 n2=60 k1=80,180,300 nsp=3 mag=1,0.7,0.5 |
     ricker1 frequency=30
     ''')
Plot('cmp', 'grey title="CMP" label1=Time unit1=s label2=Offset unit2=m')

# A velocity curve drawn on top.
Flow('vcurve', None,
     '''
     math n1=60 output="1500+x1*10" d1=1 o1=0
     ''')
Plot('vcurve',
     '''
     graph transp=y yreverse=y min2=0 max2=2000
           wantaxis=n wanttitle=n plotcol=2 plotfat=3 pad=n
     ''')

Result('cmp-with-vcurve', 'cmp vcurve', 'Overlay')

End()
