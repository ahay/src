import m8r

# top of reservoir
top = m8r.math(n1=301,d1=1,o1=0,output=0.06,label1='Trace')[0]
# bottom of reservoir
bot = m8r.math(output='input*(2-x1/300)')[top]
# wedge model
wedge=m8r.cat(axis=2).unif2(n1=181,o1=0,d1=0.001,label1='Time',unit1='s',v00=[10,15,20])[top,bot]
# seismic convolution modeling
seis = m8r.ai2refl.ricker1(frequency=25)[wedge]
# plotting
seis.grey(title='Wedge Model',color='G').show()
m8r.window(j2=10)[seis].wiggle(title="Wedge Model",poly=True,yreverse=True,transp=True).show()
