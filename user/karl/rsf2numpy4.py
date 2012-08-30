#!/usr/bin/env python
import rsf.api as rsf
import numpy as np
import scipy as sp
import  matplotlib.pyplot as plt
import sys
import os
import c_m8r as c_rsf
import m8r

# next challenge is to get input file names from command line
# use Msfin.c as example

print "program name",sys.argv[0]
print "type sys.argv=",type(sys.argv)
if len(sys.argv)>1 :
   print "parameters:",sys.argv[1:]

filenames=[]
# i do not think I want to have option to input data on stdin
#if not(os.isatty(file.fileno(sys.stdin))):
#    filenames.append("in")

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

fin={}
for filename in filenames:
   fin[filename]=m8r.Input(filename)

nsubplot=4

(n3,n2,n1)=fin[filename].shape()
for subplot in range(nsubplot):
   plt.subplot(1,nsubplot+1,subplot+1)
   plt.imshow(fin[filenames[0]][int(n3*float(subplot)/nsubplot),:,:].T)


plt.show()
