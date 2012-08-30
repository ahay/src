#!/usr/bin/env python
import m8r
import matplotlib.pyplot as plt

model = m8r.Input('model.rsf')

model.grey(color='j',gainpanel='a',title='RSF').show()

plt.imshow(model[0,:,:])
plt.show()
