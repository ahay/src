from rsf.proj import *
import math

#--------------------------#
#Basic Utils - Raanan Dafni#
#--------------------------#

#*************************************************************************************************#
#******************************************** VELOCITY *******************************************#
#*************************************************************************************************#

#-------------------------------------------------------------------------------
# Velocity Model Building
def MakeVel (fileOut='velModel', **params) :

	Flow (fileOut, None,
		'''
		math n1=%(nx)d o1=0 d1=%(dx)f output="%(Z0)f+%(dzdx)f*x1" |
		unif2 n1=%(nz)d d1=%(dz)f n2=%(nx)d d2=%(dx)f v00=%(V1)f,%(V2)f |
		put n3=%(ny)d d3=%(dy)f o3=%(origy)f
		''' % params )

#-------------------------------------------------------------------------------
# Square Velocity
def SquareVel (fileOut='csqModel', fileIn='velModel') :

	Flow (fileOut, fileIn,
		'''
		mul ${SOURCES[0]} | put data_type=csq
		''')

#-------------------------------------------------------------------------------
# Smooth Velocity
def SmoothVel (fileOut='csqModelSm', fileIn='csqModel', **params) :

	Flow (fileOut, fileIn,
		'''
		smooth rect1=%(rect1)d rect2=%(rect2)d repeat=%(repeat)d 
		''' % params )

#-------------------------------------------------------------------------------
# Pertub Velocity
def PerturbVel (fileOut='csqPerturb', fileIn1='csqModel', fileIn2='csqModelSm') :

	Flow (fileOut, [fileIn1, fileIn2], 
		'''
		add scale=1,-1 ${SOURCES[1]}
		''')

#-------------------------------------------------------------------------------
# Inverse Velocity
def InvVel (fileOut='SlModel', fileIn='velModel') :

	Flow(fileOut, fileIn,
	     	'''
	     	math output=1/input
	     	''')

#-------------------------------------------------------------------------------
# Transpose Velocity
def TranspVel (fileOut='transpModel', fileIn='Model') :

	Flow(fileOut, fileIn,
	     	'''
	     	transp plane12 |
		put id1=1 id2=0
	     	''')

#-------------------------------------------------------------------------------
# Transpose Extended Velocity
def TranspExtVel (fileOut='transpExtModel', fileIn='Model') :

	Flow(fileOut, fileIn,
	     	'''
	     	transp plane12 |
		put id1=1 id2=0 id3=100 dim=2 gdim=3
	     	''')

#-------------------------------------------------------------------------------
# Extend Velocity
def ExtVel (fileOut1='ExtModel', fileOut2='zeroModel', fileIn='Model', **params) :

	Flow (fileOut2, None,
	    	'''
		math n1=%(nx)d o1=0 d1=%(dx)f n2=%(nh)d d2=%(dh)f o2=0 output="%(Z0)f+%(dzdx)f*x1" |
		unif3 n1=%(nz)d d1=%(dz)f n2=%(nx)d d2=%(dx)f n3=%(nh)d d3=%(dh)f o3=%(oh)f
		v00=%(V1)f,%(V2)f
		''' % params)
	Flow (fileOut1, [fileIn, fileOut2],
		'''
		cat axis=3 ${SOURCES[1]} ${SOURCES[1]}
		''')


#*************************************************************************************************#
#******************************************** WAVELET ********************************************#
#*************************************************************************************************#

#-------------------------------------------------------------------------------
# Create Base Wavelet
def CreateWavelet (fileOut1='wavelet', fileOut2='wavelethead', fileOut3='wavelet_base.su', **params) :

	#Gaussian derivative wavelet
	Flow (fileOut1, None,
		'''
		math n1=%(nt)d d1=%(dt)f output="%(f_s1)g*(x1-%(t0)g)*exp(-%(f_s2)g*(x1-%(t0)g)^2)" |
		put o1=%(orig)f
		''' % params )

	#create header file
	Flow (fileOut2, fileOut1,
		'''
		segyheader
		''')

	#convert to SU format
	Flow (fileOut3, [fileOut1, fileOut2],
		'''
		suwrite endian=0 tfile=${SOURCES[1]}
		''')


#*************************************************************************************************#
#************************************** ACQUSITION GEOMETRY **************************************#
#*************************************************************************************************#

#-------------------------------------------------------------------------------
# Design Shot-Geophone Geometry
def DesignShotGeoGeom (fileOut1='geomGrid', fileOut2='SGMap', **params) :

	Flow (fileOut1, None, 
		'''
		spike n1=%(nt)d d1=%(dt)f n2=%(ng)d o2=%(origg)d d2=%(dg)d
		      n3=%(ns)d o3=%(origs)d d3=%(ds)d mag=0
		''' % params)

	Flow (fileOut2, fileOut1,'window n1=1')

#-------------------------------------------------------------------------------
# Create Shot-Gathers Headers
def CreateSGHeaders (fileOut1='SGhead', fileOut2='SGhead.su', fileIn1='geomGrid', fileIn2='SGMap', **params) :

	#create the header fields
	if params['ns'] == 1:
		Flow ('tracl', fileIn2, 'math output=x1                           | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output=%(sht0)f                     | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)" | dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')
	else:
		Flow ('tracl', fileIn2, 'math output="x1+%(ng)d*(x2-1)"                             | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output="%(sht0)f+%(shtint)f*(x2-1)"                   | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)+%(shtint)f*(x2-1)" | dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')

	#create header file
	Flow (fileOut1, [fileIn1, 'tracl', 'selev', 'gelev', 'sx', 'gx', 'ofst'],
		'''
		segyheader tracl=${SOURCES[1]} selev=${SOURCES[2]}  gelev=${SOURCES[3]}
			      sx=${SOURCES[4]}    gx=${SOURCES[5]} offset=${SOURCES[6]}
		''')

	#convert to SU format
	Flow (fileOut2, [fileIn1, fileOut1],
		'''
		suwrite tfile=${SOURCES[1]} endian=0
		''')

#-------------------------------------------------------------------------------
# Create Towed Streamer Source Trces
def CreatTowedArray (fileOut='waveletSG.su', fileIn1='wavelet_base.su', fileIn2='SGhead.su') :

	Flow (fileOut, [fileIn1, fileIn2],
		'''
		towed_array src=${SOURCES[0]} data=${SOURCES[1]} towed=${TARGETS[0]}
		''', stdin=0, stdout=0)


#*************************************************************************************************#
#******************************************* MODELING ********************************************#
#*************************************************************************************************#

# Non-Linearized Modeling
def Modeling_acd (fileOut='bornSG.su', fileIn1='waveletSG.su', fileIn2='csqModel', fileIn3='SGhead.su', **params) :

	Flow (fileOut,[fileIn1, fileIn2, fileIn3],
		'''
		/bin/cp ${SOURCES[2]} $TARGET &&
		acd deriv=%(deriv)d adjoint=%(adj)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f
		source=${SOURCES[0]} csq=${SOURCES[1]} data=$TARGET
		''' % params, 
		stdin=0,stdout=-1,workdir='bornSG.work')

#-------------------------------------------------------------------------------
# Linearized Modeling (Born)
def ModelingBorn_acd (fileOut='bornSG.su', fileIn1='waveletSG.su', fileIn2='csqModelSm', fileIn3='csqPerturb', fileIn4='SGhead.su', **params) :

	Flow (fileOut,[fileIn1, fileIn2, fileIn3, fileIn4],
		'''
	  	/bin/cp ${SOURCES[3]} $TARGET &&
		acd deriv=%(deriv)d adjoint=%(adj)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f
		source=${SOURCES[0]} csq=${SOURCES[1]} csq_d1=${SOURCES[2]} data=$TARGET
		dump_lda=%(dump_lda)d dump_term=%(dump_term)d''' % params, 
		stdin=0,stdout=-1,workdir='bornSG.work')


#*************************************************************************************************#
#******************************************* MIGRATION *******************************************#
#*************************************************************************************************#

# RTM
def RTM_acd (fileOut='Mig', fileIn1='bornSG.su', fileIn2='csqModelSm', fileIn3='waveletSG.su' , fileIn4='csqPerturbExtTransp', **params) :

	Flow(fileOut,[fileIn1, fileIn2, fileIn3, fileIn4],
     		'''
     		scale < ${SOURCES[3]} > $TARGET dscale=0.0 &&
     		acd deriv=%(deriv)d adjoint=%(adj)d nsnaps=%(nsnaps)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f 
     		data=${SOURCES[0]} csq=${SOURCES[1]} source=${SOURCES[2]} csq_b1=$TARGET sampord=1
     		''' % params,
		stdin=0,stdout=-1,workdir='migSG.work')


