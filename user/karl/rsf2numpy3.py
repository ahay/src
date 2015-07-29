#!/usr/bin/env python
import rsf.api as rsf
import numpy as np
import scipy as sp
import  matplotlib.pyplot as plt
import sys
import os
import c_m8r as c_rsf

# next challenge is to get input file names from command line
# use Msfin.c as example

print "program name",sys.argv[0]
print "type sys.argv=",type(sys.argv)
if len(sys.argv)>1 :
   print "parameters:",sys.argv[1:]

filenames=[]
if not(os.isatty(file.fileno(sys.stdin))):
    filenames.append("in")

for parameter in sys.argv[1:]:
    print "processing parameter",parameter
    if parameter.find("=")==-1 :
        print "no = in parameter"
        filenames.append(parameter)

if len(filenames)<1:
    print "just to help me test, if there are no files in the list, I will"
    print "append the file:"
    print "/home/karl/m8r/madagascar-1.3/user/karl/model.rsf"
    filenames.append('/home/karl/m8r/madagascar-1.3/user/karl/model.rsf')

print "list of file names:",filenames


fin=rsf.Input(filenames[0])


#data = np.zeros(fin.shape(),'f')
#print "data array shape=",data.shape
#print "\n"
#fin.read(data)
#c_rsf.sf_floatread(np.reshape(data,(data.size)),c_rsf.sf_input(fin.tag))
#to get the name of the data file in rsf: fin.string("in")

if not(fin.type == "float"):
    print "oh no!! data is not float"
else:
    print "all right!  data is float.  I can handle that."

data=np.memmap(fin.string("in"), dtype=np.float32, mode="r",shape=fin.shape())
fin1=rsf.Input('/home/karl/m8r/madagascar-1.3/user/karl/model.rsf')
data1 = np.zeros(fin1.shape(),'f')
fin1.read(data1)

plt.subplot(1,2,1)
plt.imshow(data[0,:,:].T) # good left plot has shape reversed

plt.subplot(1,2,2)
plt.imshow(data1[0,:,:].T) # bad right plot has shape directly from fin.Input 

plt.show()
