from rsf.proj import *

# plot model fields
Result('vp2d_5m', './vp2d_5m/vp2d_5m.rsf','put label1=Depth unit1=m label2 = Distance unit2=m label=Velocity unit=m/ms | grey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y')

Result('dn2d_5m', './dn2d_5m/dn2d_5m.rsf','put label1=Depth unit1=m label2 = Distance unit2=m label=Density unit=g/cm^3 | grey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y')

# 2-4 results
traces = []
gathers =['demo20m', 'demo10m', 'demo5m', 'demo2.5m']
for r in range(4):
    data = './' + gathers[r] + '.rsf'
    Flow(data,'./' + gathers[r] + '/data.su','suread < $SOURCE read=data endian=0',stdin=0)
    trace = 'trace' + str(r)
    Flow(trace,data,'window n2=1 f2=100')
    traces.append(trace)
Result('data5m','./demo5m.rsf','grey title="Dome Data, 2-4 scheme, dx = dz = 5 m"')
Result('trace',traces,
   '''
   cat axis=2 ${SOURCES[1:4]} | 
   graph plotcol=4,2,6,7 wanttitle=n label2=Pressure unit2=MPa
   ''')
Result('wtrace',traces,
   '''
   cat axis=2 ${SOURCES[1:4]} | window min1=1.8 max1=2.5 |
   graph plotcol=4,2,6,7 wanttitle=n label2=Pressure unit2=MPa
   ''')

# 2-8 results
traces8k = []
gathers8k =['demo20m8k', 'demo10m8k', 'demo5m8k', 'demo2.5m8k']
for r in range(4):
    data = './' + gathers8k[r] + '.rsf'
    Flow(data,'./' + gathers8k[r] + '/data.su','suread < $SOURCE read=data endian=0',stdin=0)
    trace = 'trace8k' + str(r)
    print trace
    Flow(trace,data,'window n2=1 f2=100')
    traces8k.append(trace)
Result('data5m8k','./demo5m8k.rsf','grey title="Dome Data, 2-8 scheme, dx = dz = 5 m"')
Result('trace8k',traces8k,
   '''
   cat axis=2 ${SOURCES[1:4]} | 
   graph plotcol=4,2,6,7 wanttitle=n label2=Pressure unit2=MPa
   ''')
Result('wtrace8k',traces8k,
   '''
   cat axis=2 ${SOURCES[1:4]} | window min1=1.8 max1=2.5 |
   graph plotcol=4,2,6,7 wanttitle=n label2=Pressure unit2=MPa
   ''')

# enable rsf.proj functionality
End()
