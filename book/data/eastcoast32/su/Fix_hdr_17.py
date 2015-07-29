#!/usr/bin/env python
import sys
for line in sys.stdin.readlines():  
    tokens=line.split()
#    print "line=" + line  + "tokens=" + str(tokens)
    tracl =long(tokens[0])
    fldr  =long(tokens[1])
    tracf =long(tokens[2]) 

#    print "old =" + str(tracl) + " " + str(fldr) + " " + str(tracf)
    if fldr<700 : 
        fldr += 200
#    print "new =" + str(tracl) + " " + str(fldr) + " " + str(tracf)
        
    print str(tracl)  \
        + " " + str(fldr)      \
        + " " + str(tracf) + " "

 
