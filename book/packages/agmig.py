# Testing angle-gather time migration

from rsfproj import *

def agmig(name,dx=0.008):
    # Plot model
    Plot(name,'grey label1=Time label2=Space unit1=s unit2=km title=Model')
    # Modeling
    data = 'data-'+name
    Flow(data,name,
         '''
         halfint inv=1 |
         preconstkirch zero=y inv=y h0=0 dh=%g nh=48 vel=1.5 |
         window
         ''' % dx)
    zero = 'zero-'+name
    Plot(zero,data,
         '''
         window n3=1 f3=0 |
         grey label1=Time label2=Midpoint unit1=s unit2=km
         title="Zero Offset"
         ''')
    far = 'far-'+name
    Plot(far,data,
         '''
         window n3=1 f3=47 |
         grey label1=Time label2=Midpoint unit1=s unit2=km
         title="Far Offset"
         ''')

    Result(data,[name,zero,far],'SideBySideIso')
    # Common-offset migration
