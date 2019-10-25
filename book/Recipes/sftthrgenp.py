from rsf.proj import *

import math

# only part need to change
# while changing from one case study
# to another

### tapering after each iteration 
### to get rid of edge effects
### not sure if it is efficient
tap1=100
tap2=100

# catenating files and making a movie over iterations 
def plotmovie(it,title,display,extra=''):
	return  '''
		rcat axis=3 ${SOURCES[1:%d]} |
		put o3=1 d3=1 n3=%g label3="iteration" | window min1=%g max1=%g min2=%g max2=%g |
		grey gainpanel=each wanttitle=n scalebar=n title="%s" %s
		'''%(it,it,display[0],display[1],display[2],display[3],title,extra)


# Scaled plot for data misfit
def plotgrl2d(it,title,extra=''):
	return  '''
		rcat axis=1 ${SOURCES[1:%d]} | scale axis=1 |
		put o1=1 d1=1 n1=%g label1="iteration" |
		graph label1=iteration label2="l2 data residual" unit1= unit2= wanttitle=y min2=0.0 max2=1.0
		title="%s" %s
		'''%(it,it,title,extra)

# Not scaled plot for data misfit
def plotgrl2dns(it,title,extra=''):
	return  '''
		rcat axis=1 ${SOURCES[1:%d]} | math output="input*100" |
		put o1=1 d1=1 n1=%g label1="iteration" |
		graph label1=iteration label2="l2 data residual/l2d*100" unit1= unit2= wanttitle=y
		title="%s" %s
		'''%(it,it,title,extra)

# Plotting model update
def plotgrl2m(it,title,extra=''):
	return 	'''
		rcat axis=1 ${SOURCES[1:%d]} | scale axis=1 | 
		put o1=1 d1=1 n1=%g label1="iteration" |
		graph label1=iteration label2="l2 model residual" unit1= unit2= wanttitle=n title="%s" min2=0.0 max2=1.0
                %s
		'''%(it-1,it-1,title,extra)

# Plotting difference btw current model and specified true model
def plotgrl2mt(it,title,extra=''):
	return 	'''
		rcat axis=1 ${SOURCES[1:%d]} | scale axis=1 | 
		put o1=1 d1=1 n1=%g label1="iteration" |
		graph label1=iteration label2="l2 model residual" unit1= unit2= wanttitle=y title="%s" min2=0.0 max2=1.0
                %s
		'''%(it,it,title,extra)

# Plot final result
def plotfinalx(title,display):
	return '''
               grey title="%s" min1=%g max1=%g min2=%g max2=%g
	       '''%(title,display[0],display[1],display[2],display[3])

def sftthr(op1,op1cg,outer,inner,name,data,model,thr0,dthr,plottitle,soft,clip,display,norm,op2='',truemod=None,normbydata=False,lcurve=False):
	
	# op1      - operator for forward modelling. It has the following structure:
	#            op1 = (velocity_for_migration.rsf, dip_for_PWD.rsf, 'forward modeling string', 'dependencies format 1', 'dependencies format 2',' | precond') 
	# op1cg    - operator for conjugate gradients. It has the following structure:
	#            op1cg = (velocity_for_migration.rsf, dip_for_PWD.rsf, 'operator for CG', 'dependencies format 1', 'dependencies format 2',' | precond')
	#            precond - additional shaping to avoid instabilities: e.g. cosine tapering and bandpass; make sure it has ' | '
	# outer    - number of outer iterations
	# inner    - number of inner iterations
	# name     - preceding suffix of internal files
	# data     - data to fit
	# model    - starting model, corresponds to reflectivity distribution, usually zero
	# thr0     - first threshold
	# dthr     - increment
	# plottitl - graphs' and plots' names
	# soft     - if use sfthreshold with pclip parameter, if no use sfthr with mode=hard
	# clip     - clip value for plotting difference between real and synthetic data
	# display  - windowing parameters organized in a list display = [min1,max1,min2,max2]
	# norm     - windowing for norm estimation (to get rid of edge effects) norm = [min1,max1,min2,max2]
	# op2      - second forward modelling operator (when we run a chain we use it for kirchhoff)
	# truemod  - true model if available
	# normbyda - if normalize l2 data misfit by l2 data norm. Not sure if it is reasonable (makes sense if misfit is almost 0% of data) 
	# lcurve   - disables plotting and other functionality not necessary for lcurve analysis

	movie_x = [] # movie in the model space
	movie_diffr = [] # movie of Kirchhoff modelling of the model space
	movie_data_pi = [] # movie of dataspace, which is path integral
	movie_l2 = [] # residual estimation
	movie_pdd = [] # difference between pi data and synthetic pi data

	# soft thresholding iterations
	# applied after conjgrad optimization
	# x_k+1 = STHR(conjgrad(x_k))

	# Ax = d
	# where d is data
	# A is op1

	# A(x0 + y) = data
	# Ay = data - Ax0 = d_0
	# y = result

	# x = sftthr(x0 + y)

	# x0 is a starting model

	shrp = thr0

	if (not lcurve):

		# Normalize convergence curves by l2-data norm

		if normbydata :

			Flow(name+'-l2dnorm',data,
				'''
				math output="input*input" |
				stack axis=2 norm=n |
				stack axis=1 norm=n |
				math output="sqrt(input)"
				''')

		else : 

			Flow(name+'-l2dnorm',data,
				'''
				window n2=1 n1=1 | math output="1.0"
				''')

	for i in range(outer):

		if i==0:
			# read starting model and apply threshold to it
			#?# or thr should be applied after the first iteration?
			#?# does not matter if the input model is 0
			
			# if perform percentage clip
			#if soft == 1 :
			#	Flow(name+'-'+model+'-0',model,'threshold pclip=%g'%shrp)

			# if perform value clip
			#if soft !=1 :
			#	Flow(name+'-'+model+'-0',model,'thr thr=%g mode=hard'%shrp)

			#prev = name+'-'+model+'-0'

			prev = model

		if i>0:
			prev = name+'-x-it-%d'%(i-1)


		x = name+'-x-it-%d'%i

		result = name+'-result-it-%d'%i

		d_0 = name+'-d0-it-%d'%i

		# data - Ax0 = d_0		
		
		Flow(d_0,[prev,data,op1[0],op1[1]],op1[2] + op1[3] + '| add scale=-1,1 ${SOURCES[1]}')
	
		# forward pathint of the model - d_pi (pi of diffractions)
		# Ay = data - Ax0 = d_0
		# y = result
	
		Flow(result,[d_0,op1cg[0],op1cg[1]],
			'''
			conjgrad 
			%s %s
			niter=%d mod=$SOURCE'''%(op1cg[2],op1cg[4],inner))
		#?# mod = model dimensions

		# thr(y + x0) = x

		# if perform soft thresholding based on
		# % of coefficients surviving
		# might be non-linear procedure
		if soft==1 :
			Flow(x,[result,prev],
			     '''
                             add scale=1,1 ${SOURCES[1]} |
                             threshold pclip=%g %s
                             '''%(shrp,op1cg[5]))

		# if perform value clip		
		else :
			Flow(x,[result,prev],
                             '''
                             add scale=1,1 ${SOURCES[1]} |
                             thr thr=%g mode=hard %s
                             '''%(shrp,op1cg[5]))

		shrp = shrp + dthr
