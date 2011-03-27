from rsf.proj import Flow

Flow('spike',None,'spike n1=1000 n2=100 | bandpass fhi=10')
Flow('cliped','spike','clip clip=0.5')
