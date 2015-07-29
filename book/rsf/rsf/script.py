>>> import m8r
>>> spike = m8r.spike(n1=1000,n2=100)[0]
>>> spike
<m8r.File object at 0x4038b10>
>>> m8r.clip(clip=0.5)
<m8r.Filter object at 0x9976690>
>>> cliped = m8r.clip(clip=0.5)[spike]
>>> cliped2 = m8r.spike(n1=1000,n2=100).clip(clip=0.5)[0]
>>> import numpy
>>> cliped = numpy.clip(spike,-0.5,0.5)
