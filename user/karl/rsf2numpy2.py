import rsf.api as rsf
import numpy as np
import scipy as sp
import  matplotlib.pyplot as plt

fin=rsf.Input('/home/karl/m8r/madagascar-1.3/user/karl/model.rsf')

print "fin.shape()=",fin.shape()
print "list(fin.shape())=",list(fin.shape())
print "(list(fin.shape())).reverse()=",(list(fin.shape())).reverse()
mylist=list(fin.shape())
mylist.reverse()

data = np.zeros(mylist,'f')
print "data array shape=",data.shape
print "\n"
fin.read(data)
print "data array shape=",data.shape
print "data[0,:,:].shape=",data[0,:,:].shape

# in python and c you index data[i3,i2,i1].   Shape need to be (n3,n2,n1), 
# reverse from that it in fin.shape()


fin1=rsf.Input('/home/karl/m8r/madagascar-1.3/user/karl/model.rsf')
data1 = np.zeros(fin1.shape(),'f')
fin1.read(data1)

plt.subplot(121)
plt.imshow(data[0,:,:].T) # good left plot has shape reversed

plt.subplot(122)
plt.imshow(data1[0,:,:].T) # bad right plot has shape directly from fin.Input 

plt.show()
