from rsf.proj import *

Flow('synth', None,
     '''
     spike n1=200 n2=100 k1=50,100,150 nsp=3 mag=1,0.8,0.6 |
     ricker1 frequency=40 |
     noise seed=2025 var=0.005
     ''')

Result('synth',
       '''
       grey title="Synthetic gather"
            label1="Time" unit1="s" label2="Trace" unit2=""
            color=j scalebar=y
       ''')

End()
