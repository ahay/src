#!/usr/bin/env python
'''quick view of binary data (1D/2D/3D), usage: sfviewbin data.bin n1 n2 n3
'''

##   Copyright (C) 2019 
##   coded by Hongyu Zhou @ China University of Petroleum 
##   (Email:carrot_cup@163.com)
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
import os, sys
#import m8r
print ("\nquick view of binary data (1D/2D/3D), Usage: \n\n\n 1D: sfviewbin data.bin n1 \n 2D: sfviewbin data.bin n1 n2 \n 3D: sfviewbin data.bin n1 n2 n3\n") 
if len(sys.argv)>=2: 
	filename=sys.argv[1];
	filename=filename+".rsf"

	f = open(filename,'w' )

	f.write('in="%s"\n' %sys.argv[1])
	f.write('data_format="native_float"\n')
	f.write('esize=4\n')

	if len(sys.argv)<3:            # 
   		print ("more parameters required\n") 
    	  
	if len(sys.argv) == 3:
		f.write('n1=%s\n' %sys.argv[2])
	if len(sys.argv) == 4:
		f.write('n1=%s\n' %sys.argv[2])
		f.write('n2=%s\n' %sys.argv[3])
	if len(sys.argv) == 5:
		f.write('n1=%s\n' %sys.argv[2])
		f.write('n2=%s\n' %sys.argv[3])
		f.write('n3=%s\n' %sys.argv[4])
	f.close()




if len(sys.argv) == 3:
	os.system("sfgraph<%s|sfpen  bgcolor=w" %filename)

if len(sys.argv) == 4:
	os.system("sfgrey<%s color=j|sfpen  bgcolor=w" %filename)

if len(sys.argv)==5:
	n1_half=int(sys.argv[2])/2
	n2_half=int(sys.argv[3])/2
	n3_half=int(sys.argv[4])/2
	print("n1_half=%d" %n1_half) 
	print("n2_half=%d" %n2_half) 
	print("n3_half=%d" %n3_half) 
	os.system("sfbyte < %s gainpanel=a bar=bar.rsf | sfgrey3  frame1=%d frame2=%d frame3=%d color=j flat=n scalebar=y  | sfpen  bgcolor=w" %(filename, n1_half, n2_half,n3_half))
