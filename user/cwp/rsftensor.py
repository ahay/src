try:    
  from rsf.cluster import *
  from rsf.proj import WhereIs
except: from rsf.proj    import *
import os

JHOME=WhereIs('java')


# ------------------------------------------------------------
# ------------------------------------------------------------
#
# Function Definitions
#
# ------------------------------------------------------------
# ------------------------------------------------------------

#################################################################
## Computes the structure tensors
def metricTensors3D(ten3D,img3D,hw1,hw2,sig1,sig2=0.0,sig3=0.0,jmem='2g'):
	Flow(ten3D,img3D,
		'''
		%s -server -ea -Xmx%s cip.RSFTensor
		hw1=%d hw2=%d sig1=%f sig2=%f sig3=%f 
		opt=metric dim=3
		''' %(JHOME,jmem,hw1,hw2,sig1,sig2,sig3))

def metricTensors2D(ten2D,img2D,hw1,hw2,sig1,sig2=0.0,jmem='2g'):
	Flow(ten2D,img2D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 hw1=%d hw2=%d sig1=%f sig2=%f opt=metric dim=2''' 
			 %(JHOME,jmem,hw1,hw2,sig1,sig2))

def structureTensors3D(ten3D,img3D,sig1,sig2=0.0,sig3=0.0,jmem='2g'):
	Flow(ten3D,img3D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f sig3=%f opt=struc dim=3''' 
			 %(JHOME,jmem,sig1,sig2,sig3))

def structureTensors2D(ten2D,img2D,sig1,sig2=0.0,jmem='2g'):
	Flow(ten2D,img2D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f opt=struc dim=2''' 
			 %(JHOME,jmem,sig1,sig2))

def orientedTensors3D(ten3D,img3D,sig1,sig2=0.0,sig3=0.0,jmem='2g'):
	Flow(ten3D,img3D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f sig3=%f opt=orien dim=3''' 
			 %(JHOME,jmem,sig1,sig2,sig3))

def orientedTensors2D(ten2D,img2D,sig1,sig2=0.0,jmem='2g'):
	Flow(ten2D,img2D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f opt=orien dim=2''' 
			 %(JHOME,jmem,sig1,sig2))

def exclusionTensors3D(ten3D,img3D,sig1,sig2=0.0,sig3=0.0,ecc=1.0,jmem='2g'):
	Flow(ten3D,img3D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f sig3=%f opt=excld ecc=%f dim=3''' 
			 %(JHOME,jmem,sig1,sig2,sig3,ecc))

def exclusionTensors2D(ten2D,img2D,sig1,sig2=0.0,ecc=1.0,jmem='2g'):
	Flow(ten2D,img2D,
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=%f sig2=%f opt=excld ecc=%f dim=2''' 
			 %(JHOME,jmem,sig1,sig2,ecc))

def getTensorsAtCoord2D(tenC2D,tenAll2D,coord2D,par,jmem='2g'):
	Flow(tenC2D,[tenAll2D, coord2D],
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=1 sig2=1 opt=coord dim=2 coordf=${SOURCES[1]}
			 o1=%f o2=%f d1=%f d2=%f''' 
			 %(JHOME,jmem,par['o1'],par['o2'],par['d1'],par['d2']))

def getTensorsAtCoord3D(tenC3D,tenAll3D,coord3D,par,jmem='2g'):
	Flow(tenC3D,[tenAll3D, coord3D],
			 '''%s -server -ea -Xmx%s cip.RSFTensor
			 sig1=1 sig2=1 sig3=1 opt=coord dim=3 coordf=${SOURCES[1]}
			 o1=%f o2=%f o3=%f d1=%f d2=%f d3=%f''' 
			 %(JHOME,jmem,
			   par['o1'],par['o2'],par['o3'],
				 par['d1'],par['d2'],par['d3']))
