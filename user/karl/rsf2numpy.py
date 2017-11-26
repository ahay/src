import rsf.api as rsf
import numpy as np
import scipy as sp
import  matplotlib.pyplot as plt

fin=rsf.Input('/home/karl/m8r/madagascar-1.3/user/karl/model.rsf')

class Header:
    def __init__(self, RSFFile):
        self.dims = {}
        self.shape = []
        for i in range(1,4):
            nkey = 'n%d'%i
            dkey = 'd%d'%i
            okey = 'o%d'%i
            self.dims[nkey] = RSFFile.int(nkey)
            self.dims[dkey] = RSFFile.float(dkey)
            self.dims[okey] = RSFFile.float(okey)
            self.shape.append(RSFFile.int(nkey))
        self.shape.reverse()
#header = fin.shape() # Header(fin)
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

p=plt.imshow(data[0,:,:])

plt.show()
