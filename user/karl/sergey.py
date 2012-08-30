import m8r
import matplotlib.pyplot as plt

model = m8r.math(n1=201,o1=-6.283185,d1=.06283185,
                 n2=201,o2=-6.283185,d2=.06283185,
                 n3=101,o3=-3.141592,d3=.06283185,
                 output='sin(x1*x3) * cos(x2)*sin(x3) + sin(0.5*x3)')[0]
model.grey(color='j',gainpanel='a',title='RSF').show()

plt.imshow(model[0,:,:])
plt.show()
