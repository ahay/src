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
## computes the linear semblance using structure oriented 
## semblance, and their corresponding eigenvalues. 
## See MinesJTK for more detail
def linearSemblance2D(lsem2D,img2D,hw1,hw2,sig1,sig2=0.0,zclip=0.0,jmem='2g'):
	Flow(lsem2D,img2D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLSF
	     hw1=%d hw2=%d sig1=%f sig2=%f zclip=%f opt=lin dim=2
	     '''%(JHOME,jmem,hw1,hw2,sig1,sig2,zclip))


#################################################################
## computes the linear semblance using structure oriented 
## semblance, and their corresponding eigenvalues. 
## See MinesJTK for more detail
def linearSemblance3D(lsem3D,img3D,hw1,hw2,sig1,sig2=0.0,sig3=0.0,zclip=0.0,jmem='2g'):
	Flow(lsem3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLSF
	     hw1=%d hw2=%d sig1=%f sig2=%f sig3=%f zclip=%f opt=lin dim=3
	     '''%(JHOME,jmem,hw1,hw2,sig1,sig2,sig3,zclip))


#################################################################
## computes the planar semblance using structure oriented 
## semblance, and their corresponding eigenvalues. 
## See MinesJTK for more detail
def planarSemblance3D(psem3D,img3D,hw1,hw2,sig1,sig2=0.0,sig3=0.0,zclip=0.0,jmem='2g'):
	Flow(psem3D,img3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFLSF
	     hw1=%d hw2=%d sig1=%f sig2=%f sig3=%f zclip=%f opt=pln dim=3
	     '''%(JHOME,jmem,hw1,hw2,sig1,sig2,sig3,zclip))
