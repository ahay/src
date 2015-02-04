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
## computes the linearity image using structure tensors, and
## their corresponding eigenvalues. See MinesJTK for more detail
def linearity2D(lin2D,img2D,sig1,sig2=0.0,zclip=0.0,jmem='2g'):
	Flow(lin2D,img2D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f zclip=%f opt=lin dim=2
	     '''%(JHOME,jmem,sig1,sig2,zclip))

def linearity3D(lin3D,img3D,sig1,sig2=0.0,sig3=0.0,zclip=0.0,jmem='2g'):
	Flow(lin3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f sig3=%f zclip=%f opt=lin dim=3
	     ''' %(JHOME,jmem,sig1,sig2,sig3,zclip))


#################################################################
## computes the planarity image using structure tensors, and
## their corresponding eigenvalues. See MinesJTK for more detail
def planarity3D(pln3D,img3D,sig1,sig2=0.0,sig3=0.0,zclip=0.0,jmem='2g'):
	Flow(pln3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f sig3=%f zclip=%f opt=pln dim=3
	     ''' %(JHOME,jmem,sig1,sig2,sig3,zclip))


#################################################################
## Computes the eigen tensors
def eigenTensors3D(eig3D,img3D,sig1,sig2=0.0,sig3=0.0,jmem='2g'):
	Flow(eig3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f sig3=%f opt=tensor dim=3
	     ''' %(JHOME,jmem,sig1,sig2,sig3))

def eigenTensors2D(eig2D,img2D,sig1,sig2=0.0,jmem='2g'):
	Flow(eig2D,img2D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f opt=tensor dim=2
	     ''' %(JHOME,jmem,sig1,sig2))
	

#################################################################
## Get the normal vectors for the entire image
## 2D: Normals are stored x1,x2,u[1,2] where x1 is fastest in memory
## 3D: Normals are stored x1,x2,x3,u[1-3] where x1 is fastest in memory
def getAllNormal2D(norm2D,img2D,sig1,sig2=0.0,jmem='2g'):
	Flow(norm2D,img2D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f opt=nvec dim=2
	     ''' %(JHOME,jmem,sig1,sig2))


def getAllNormal3D(norm3D,img3D,sig1,sig2=0.0,sig3=0.0,jmem='2g'):
	Flow(norm3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=%f sig2=%f sig3=%f opt=nvec dim=3
	     ''' %(JHOME,jmem,sig1,sig2,sig3))


#################################################################
## Get the normal vectors at coordinates
def getNormalAtCoord2D(normC2D,normAll2D,coord2D,par,jmem='2g'):
	Flow(normC2D,[normAll2D, coord2D],
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=1 sig2=1 opt=natc dim=2 coordf=${SOURCES[1]}
	     o1=%f o2=%f d1=%f d2=%f
	     '''%(JHOME,jmem,par['o1'],par['o2'],par['d1'],par['d2']))


def getNormalAtCoord3D(normC3D,normAll3D,coord3D,par,jmem='2g'):
	Flow(normC3D,[normAll3D, coord3D],
	     '''
	     %s -server -ea -Xmx%s cip.RSFLOF
	     sig1=1 sig2=1 sig3=1 opt=natc dim=3 coordf=${SOURCES[1]}
	     o1=%f o2=%f o3=%f d1=%f d2=%f d3=%f
	     '''%(JHOME,jmem,
		  par['o1'],par['o2'],par['o3'],
		  par['d1'],par['d2'],par['d3']))
