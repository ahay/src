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
## This clips sets the edges of a image to values of zero
def applyPadding2D(odata2D,idata2D,n1,n2,pad1,pad2):
	Flow(idata2D+'padding',idata2D,
	     'spike n1=%d n2=%d k1=%d l1=%d k2=%d l2=%d'
	     %(n1,n2,pad1+1,n1-pad1,pad2+1,n2-pad2))
	Flow(odata2D,[idata2D,idata2D+'padding'],
	     'math i=${SOURCES[0]} j=${SOURCES[1]} output="i*j"')
	
def applyPadding3D(odata3D,idata3D,n1,n2,n3,pad1,pad2,pad3):
	Flow(idata3D+'padding',idata3D,
	     'spike n1=%d n2=%d n3=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d'
	     %(n1,n2,n3,pad1+1,n1-pad1,pad2+1,n2-pad2,pad3+1,n3-pad3))
	Flow(odata3D,[idata3D,idata3D+'padding'],
	     'math i=${SOURCES[0]} j=${SOURCES[1]} output="i*j"')


#################################################################
## computes the the cippicks using an exlusion radius and a 
## compound probability image
def pickCoord2D(cippicks2D,cprob2D,r1,r2=0.0,jmem='2g'):
	Flow(cippicks2D,cprob2D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFGreedyCIPPicker 
	     r1=%f r2=%f opt=coord dim=2
	     '''%(JHOME,jmem,r1,r2))

def pickCoord3D(cippicks3D,cprob3D,r1,r2=0.0,r3=0.0,jmem='2g'):
	Flow(cippicks3D,cprob3D,
	     '''
	     %s -server -ea -Xmx%s cip.RSFGreedyCIPPicker 
	     r1=%f r2=%f r3=%f opt=coord dim=3
	     '''%(JHOME,jmem,r1,r2,r3))

def pickTensorCoord2D(cippicks2D,cprob2D,tensor2D,rmax,jmem='2g'):
	Flow(cippicks2D,[cprob2D,tensor2D],
	     '''
	     %s -server -ea -Xmx%s cip.RSFGreedyCIPPicker 
	     r1=%f opt=coord tfile=${SOURCES[1]} etype=tensor 
	     dim=2
	     '''%(JHOME,jmem,rmax))
	
def pickTensorCoord3D(cippicks3D,cprob3D,tensor3D,rmax,jmem='2g'):
	Flow(cippicks3D,[cprob3D,tensor3D],
	     '''
	     %s -server -ea -Xmx%s cip.RSFGreedyCIPPicker 
	     r1=%f opt=coord tfile=${SOURCES[1]} etype=tensor 
	     dim=3
	     '''%(JHOME,jmem,rmax))

#################################################################
## transposes the picks by dims (picks are stored in 2D array)
##
def transpPicksDims2D(opicks,ipicks):
	opn1 = ipicks+'-n1'
	opn2 = ipicks+'-n2'
	Flow(opn1,ipicks,'window f2=0 n2=1')
	Flow(opn2,ipicks,'window f2=1 n2=1')
	if plane=='12':
		Flow(opicks,[opn2,opn1],
		     'cat axis=2 space=n ${SOURCES[1]}')
		
def transpPicksDims3D(opicks,ipicks,plane='12'):
	opn1 = ipicks+'-n1'
	opn2 = ipicks+'-n2'
	opn3 = ipicks+'-n3'
	Flow(opn1,ipicks,'window f2=0 n2=1')
	Flow(opn2,ipicks,'window f2=1 n2=1')
	Flow(opn3,ipicks,'window f2=2 n2=1')
	if plane=='12':
		Flow(opicks,[opn2,opn1,opn3],
		     'cat axis=2 space=n ${SOURCES[1:3]}')
	elif plane=='13':
		Flow(opicks,[opn3,opn2,opn1],
		     'cat axis=2 space=n ${SOURCES[1:3]}')
	elif plane=='23':
		Flow(opicks,[opn1,opn3,opn2],
		     'cat axis=2 space=n ${SOURCES[1:3]}')
		

#################################################################
## Show ExclusionZones
##
def showTensorZones3D(oimage,iimage,coord,priority,tensor,rmax=0.0,f=0,j=1,jmem='2g'):
	Flow(oimage,[iimage,coord,priority,tensor],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     rmax=%f opt=tensor cfile=${SOURCES[1]} pfile=${SOURCES[2]} tfile=${SOURCES[3]} 
	     first=%d jump=%d dim=3
	     '''%(JHOME,jmem,rmax,f,j))
	
def showTensorZones2D(oimage,iimage,coord,priority,tensor,rmax=0.0,f=0,j=1,jmem='2g'):
	Flow(oimage,[iimage,coord,priority,tensor],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     rmax=%f opt=tensor cfile=${SOURCES[1]} pfile=${SOURCES[2]} tfile=${SOURCES[3]} 
	     first=%d jump=%d dim=2
	     '''%(JHOME,jmem,rmax,f,j))

def showEllipsoidZones3D(oimage,iimage,coord,priority,r1=0.0,r2=0.0,r3=0.0,f=0,j=1,jmem='2g'):
	Flow(oimage,[iimage,coord,priority],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     r1=%f r2=%f r3=%f opt=ellp cfile=${SOURCES[1]} pfile=${SOURCES[2]}
	     first=%d jump=%d dim=3
	     '''%(JHOME,jmem,r1,r2,r3,f,j))
	
def showEllipseZones2D(oimage,iimage,coord,priority,r1=0.0,r2=0.0,f=0,j=1,jmem='2g'):
	Flow(oimage,[iimage,coord,priority],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     r1=%f r2=%f opt=ellp cfile=${SOURCES[1]} pfile=${SOURCES[2]}
	     first=%d jump=%d dim=2
	     '''%(JHOME,jmem,r1,r2,f,j))
	
def showSingleTensorZones3D(oimage,iimage,priority,tensor,c1,c2,c3,rmax=0.0,jmem='2g'):
	Flow(oimage,[iimage,priority,tensor],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     c1=%f c2=%f c3=%f rmax=%f opt=tsingle pfile=${SOURCES[1]} tfile=${SOURCES[2]} 
	     dim=3
	     '''%(JHOME,jmem,c1,c2,c3,rmax))
	
def showSingleTensorZones2D(oimage,iimage,priority,tensor,c1,c2,rmax=0.0,jmem='2g'):
	Flow(oimage,[iimage,priority,tensor],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     c1=%f c2=%f rmax=%f opt=tsingle pfile=${SOURCES[1]} tfile=${SOURCES[2]} 
	     dim=2
	     '''%(JHOME,jmem,c1,c2,rmax))
	
def showSingleEllipsoidZones3D(oimage,iimage,priority,c1,c2,c3,r1=0.0,r2=0.0,r3=0.0,jmem='2g'):
	Flow(oimage,[iimage,priority],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     c1=%f c2=%f c3=%f r1=%f r2=%f r3=%f opt=esingle pfile=${SOURCES[1]}
	     dim=3
	     '''%(JHOME,jmem,c1,c2,c3,r1,r2,r3))
	
def showSingleEllipseZones2D(oimage,iimage,priority,c1,c2,r1=0.0,r2=0.0,jmem='2g'):
	Flow(oimage,[iimage,priority],
	     '''
	     %s -server -ea -Xmx%s cip.RSFExclusionZones 
	     c1=%f c2=%f r1=%f r2=%f opt=esingle pfile=${SOURCES[1]}
	     dim=2
	     '''%(JHOME,jmem,c1,c2,r1,r2))
	

#################################################################
## Sort picks by exclusion-zone semblance
##
def sortTensorEZSemblance3D(opicks,ipicks,priority,tensor,rmax=0.0,jmem='2g'):
	Flow(opicks,[priority,ipicks,tensor],
	  '''
	  %s -server -ea -Xmx%s cip.RSFEZSemblanceSort 
	  rmax=%f opt=tensor cfile=${SOURCES[1]} tfile=${SOURCES[2]} 
	  dim=3
	  '''%(JHOME,jmem,rmax))

def sortTensorEZSemblance2D(opicks,ipicks,priority,tensor,rmax=0.0,jmem='2g'):
	Flow(opicks,[priority,ipicks,tensor],
	  '''
	  %s -server -ea -Xmx%s cip.RSFEZSemblanceSort 
	  rmax=%f opt=tensor cfile=${SOURCES[1]} tfile=${SOURCES[2]} 
	  dim=2
	  '''%(JHOME,jmem,rmax))

def sortEZSemblance3D(opicks,ipicks,priority,r1=0.0,r2=0.0,r3=0.0,jmem='2g'):
	Flow(opicks,[priority,ipicks],
	  '''
	  %s -server -ea -Xmx%s cip.RSFEZSemblanceSort 
	  r1=%f r2=%f r3=%f opt=ellp cfile=${SOURCES[1]}
	  dim=3
	  '''%(JHOME,jmem,r1,r2,r3))

def sortEZSemblance2D(opicks,ipicks,priority,r1=0.0,r2=0.0,jmem='2g'):
	Flow(opicks,[priority,ipicks],
	  '''
	  %s -server -ea -Xmx%s cip.RSFEZSemblanceSort 
	  r1=%f r2=%f opt=ellp cfile=${SOURCES[1]}
	  dim=2
	  '''%(JHOME,jmem,r1,r2))
