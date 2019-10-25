from rsf.proj import *

import math

# This script is compatible with
# sfchainde and thresholds
# only diffraction part of the model

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

def plotmovie2(it,n2,o2,title,display,extra=''):
	return  '''
		rcat axis=3 ${SOURCES[1:%d]} |
		window n2=%d f2=%d | put o2=%g |
		put o3=1 d3=1 n3=%g label3="iteration" | window min1=%g max1=%g min2=%g max2=%g |
		grey gainpanel=each wanttitle=n scalebar=y title="%s" %s
		'''%(it,n2,n2,o2,it,display[0],display[1],display[2],display[3],title,extra)

def plotmovienw(it,title,extra=''):
	return  '''
		rcat axis=3 ${SOURCES[1:%d]} |
		put o3=1 d3=1 n3=%g label3="iteration" |
		grey gainpanel=each wanttitle=n scalebar=y title="%s" %s
		'''%(it,it,title,extra)

# Scaled plot for data misfit
def plotgrl2d(it,title,extra=''):
	return  '''
		rcat axis=1 ${SOURCES[1:%d]} | scale axis=1 |
		put o1=1 d1=1 n1=%g label1="iteration" |
		graph label1=iteration label2="l2 data residual" unit1= unit2=
		plotfat=8 plotcol=7wanttitle=y min2=0.0 max2=1.0
		title="%s" wanttitle=n %s
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
		graph label1=iteration label2="l2 model residual" unit1= unit2=
		plotfat=8 plotcol=7 wanttitle=n title="%s" min2=0.0 max2=1.0
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

def sftthr(op1,op1cg,outer,inner,name,data,model,thr0,dthr,dipim,radius,plottitle,soft,clip,nx,ox,display,norm,op2='',movie=False,truediff=None,orthogon=(False,-999.999,-999.999,-100,False),mask=None,truemod=None,normbydata=False,lcurve=False):

	# op1      - operator for forward modelling. It has the following structure:
	#            op1 = (velocity_for_migration.rsf, dip_for_PWD.rsf, 'forward modeling string', 'dependencies format 1', 'dependencies format 2',' | precond')
	# op1cg    - operator for conjugate gradients. It has the following structure:
	#            op1cg = (velocity_for_migration.rsf, dip_for_PWD.rsf, 'operator for CG', 'dependencies format 1', 'dependencies format 2',' | precond')
	#            precond - additional shaping to avoid instabilities: e.g. cosine tapering and bandpass; make sure it has ' | '
	# outer    - number of outer iterations
	# outer    - number of outer iterations
	# inner    - number of inner iterations
	# name     - preceding suffix of internal files
	# data     - data to fit
	# model    - starting model, corresponds to [ reflections, diffractions ] (usually zero)
	# thr0     - first threshold
	# dthr     - increment for threshold
	# dipim    - dip estimated in an image domain
	# radius   - regularization on reflections (spraying radius)
	# plottitl - graphs' and plots' names
	# soft     - if use sfthreshold with pclip parameter, if no use sfthr with mode=hard
	# clip     - clip value for plotting difference between real and synthetic data
	# nx, ox   - number of traces in reflection/diffraction model and starting value - needed for windowing for display
	# display  - windowing parameters organized in a list display = [min1,max1,min2,max2]
	# norm     - windowing for norm estimation (to get rid of edge effects) norm = [min1,max1,min2,max2]
	# op2      - second forward modelling operator (when we run a chain we use it for Kirchhoff)
	# movie    - if create movies
	# truediff - diffraction hyperbolas for instance with noise to compute the difference with de-noised ones
	# orthogon - tuple containing signal and noise orthogonalization parameters (used to remove reflections from a diffraction model)
	#            (perform orthogon,rect1,rect2,niter,stopper)
	# mask     - known locations of diffractions (usually determined beforehand)
	# truemod  - true model if available
	# normbyda - if normalize l2 data misfit by l2 data norm. Not sure if it is reasonable (makes sense if misfit is almost 0% of data)
	# lcurve   - disables plotting and other functionality not necessary for lcurve analysis

	movie_x = [] # movie in the model space
	movie_diffr = [] # movie of Kirchhoff modelling of the model space
	movie_noise = [] # difference btw diffr and true diffractions
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

	# diffractions - soft thresholding
	# reflections - spraying

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

		Flow(result,[d_0,model,op1[0],op1[1]],
			'''
			conjgrad
			%s
			niter=%d mod=${SOURCES[1]}'''%(op1cg[2] + op1cg[3],inner))

		# thr(y + x0) = x

		# if perform soft thresholding based on
		# % of coefficients surviving
		# might be non-linear procedure

		sumx = x+'-sum' # both models
		xthr = x+'-thr' # regularization on diffractions
		xspr = x+'-spr' # regularization on reflections

		# if apply orthogonalization only once
		if ( orthogon[4] ):

			stop = True

		# in any case orthogonalization will be applied
		# on the first iteration
		if ( i == 0 ):

			stop = False

		# diffractions and reflections orthogonalization
		# reflections leaked to the diffraction domain
		# will be moved from diffractions to reflections
		if ( orthogon[0] and ( stop == False ) ):

			resultRefl     = result+'Refl'  # Reflection part of an update
			resultReflOrth = result+'ReflO' # Orthogonalized reflection update
			resultDiffOrth = result+'DiffO' # Orthogonalized diffraction update

			# First, update is assigned to reflections only
			# recall that update to both parts of the model are the same
			# so orthogonalization will put all the energy to relfections' model
			if i == 0:

				Flow(resultRefl,result,'window n2=%d'%(nx))

			# Here we are measuring similarity between an update to diffractions
			# and previous reflection model: everything what looks like reflections
			# in the update to diffractions will be deleted (moved to resultReflOrth,
			# which is not used anywhere)
			else :

				Flow(resultRefl,prev,'window n2=%d'%(nx))

			Flow([resultDiffOrth,resultReflOrth],[result,resultRefl],
					'''
					window n2=%d f2=%d |
					ortho rect1=%d rect2=%d niter=%d
					sig=${SOURCES[1]} sig2=${TARGETS[1]}
					'''%(nx,nx,orthogon[1],orthogon[2],orthogon[3]))

			Flow(sumx,[result,resultDiffOrth,prev],
				'''
				window n2=%d |
				cat axis=2 ${SOURCES[1]} |
				add scale=1,1 ${SOURCES[2]}
				'''%(nx))

		else :

			Flow(sumx,[result,prev],
                        	 '''
                        	 add scale=1,1 ${SOURCES[1]}
                        	 ''')

		if soft==1 :

			# if we want to apply
			# diffraction location mask
			if ( mask != None ):

				Flow(xthr,[sumx,mask],
					'''
					window n2=%d f2=%d |
					math K=${SOURCES[1]} math output="input*K" |
					thr thr=%g mode=soft
					'''%(nx,nx,shrp))

			else:

				Flow(xthr,sumx,
					'''
					window n2=%d f2=%d |
					thr thr=%g mode=soft
					'''%(nx,nx,shrp))

		# if perform value clip
		else :

			# if we want to apply
			# diffraction location mask
			if ( mask != None ):

				Flow(xthr,[sumx,mask],
					'''
					window n2=%d f2=%d |
					math K=${SOURCES[1]} math output="input*K" |
					thr thr=%g mode=hard
					'''%(nx,nx,shrp))

			else:

				Flow(xthr,sumx,
					'''
					window n2=%d f2=%d |
					thr thr=%g mode=hard
					%s
					'''%(nx,nx,shrp,op1cg[5]))

		Flow(xspr,[sumx,dipim],
			'''
			window n2=%d |
			pwspray dip=${SOURCES[1]} ns=%d reduce=stack %s
			'''%(nx,radius,op1cg[5]))

		Flow(x,[xspr,xthr],
			'''
			cat axis=2 ${SOURCES[1]}
			''')

		# cooling
		shrp = shrp + dthr

		""" currently disabled - wavenumbers need to be restored be4 looking at the differences
		# diffraction hyperbolas modelling
		# should specify operator op2
		if op2 != '':

			# synthetic diffractions = op2 (diffractivity model)
			diffr = name+'-diffr-it-%d'%i

			# noise = true diffractions - synthetic diffractions
			noise = name+'-noise-it-%d'%i

			Flow(diffr,x,op2)

			# truediff = true diffractions
			if truediff!=None:
				truediffs = truediff+'s'
				Flow(truediffs,truediff,'scale axis=2')
				Flow(noise,[diffr,truediffs],'scale axis=2 | add scale=1,-1 ${SOURCES[1]}')
		"""

		# synthetic P PWD L (reflection + diffraction models)
		data_pi = name+'-data-pi-it-%d'%i

		Flow(data_pi,[x,op1[0],op1[1]],op1[2] + op1[4] + ' | put o2=%g'%ox)

		# Computing difference between path-integral data and synthetic path-integral data (pdd)

		pdd = name+'-pdd-it-%d'%i

		Flow(pdd,[data,data_pi],
			'''
			add scale=1,-1 ${SOURCES[1]}
			''')

		# creating movies
		if ((not lcurve) and movie):

			movie_x.append(x)
			"""
			if op2 != '':

				movie_diffr.append(diffr)
				movie_noise.append(noise)
			"""
			movie_data_pi.append(data_pi)

			movie_pdd.append(pdd)

			# computing the residual in the data domain
			l2 = name+'-l2-it-%d'%i

			# accounts for extended model and data dimensions
			# scaled by l2 data norm if normbydata=True
			Flow(l2,[data,data_pi,name+'-l2dnorm'],
				'''
				add scale=1,-1 ${SOURCES[1]} |
				window min1=%g max1=%g min2=%g max2=%g |
				math output="input*input" |
				stack axis=2 norm=n |
				stack axis=1 norm=n |
				math K=${SOURCES[2]} output="sqrt(input)/K"
				'''%(norm[0],norm[1],norm[2],norm[3]))

			movie_l2.append(l2)

		#||| for end

	if ((not lcurve) and movie):

		print('Following movies will be generated')
		print('model movie %s'%('a-'+name+'-movie-x'))
		#print 'diffr movie %s'%('a-'+name+'-movie-diffr')
		print('data pi movie %s'%('a-'+name+'-movie-dpi'))

		# diffractivity only
		#Result('a-'+name+'-movie-x',movie_x,plotmovie2(outer,nx,ox,plottitle,display))

		# diffractivity and reflectivity
		Plot(name+'-movie-x',movie_x,plotmovienw(outer,plottitle))

		"""
		if op2 != '':
			# diffractions
			Result('a-'+name+'-movie-diffr',movie_diffr,plotmovie(outer,plottitle,display))
			if truediff!=None:
				# difference with true diffractions
				Result('a-'+name+'-movie-noise',movie_noise,plotmovie(outer,plottitle,display))
		"""
		# data domain
		Plot('a-'+name+'-movie-dpi',movie_data_pi,plotmovie(outer,plottitle,display))

		# data domain difference
		Plot('a-'+name+'-movie-pdd',movie_pdd,plotmovie(outer,plottitle,display,
					'''
					min1=%g max1=%g min2=%g max2=%g gainpanel=all clip=%g
					'''%(norm[0],norm[1],norm[2],norm[3],clip)))

		Plot('a-'+name+'-movie-l2',movie_l2,plotgrl2d(outer,plottitle))

		#Plot('a-'+name+'-movie-l2',movie_l2,plotgrl2d(outer,plottitle))

		# Not scaled plot of l2d misfit

		#Result('a-'+name+'-movie-l2ns',movie_l2,plotgrl2dns(outer,plottitle))

		# sequential iterations model update estimation
		# done for both model spaces (reflections and diffractions) separately
		movie_l2m = []
		movie_l2mr = []

		# we do not need to subtract the final model
		# from the final model
		for j in range (outer - 1):

			l2m = name+'-l2m-it-%d'%j
			l2mr = name+'-l2mr-it-%d'%j

			# accounts for extended model and data dimensions
			Flow(l2m,[name+'-x-it-%d'%j,name+'-x-it-%d'%(j+1)],
				'''
				add scale=1,-1 ${SOURCES[1]} | window n2=%d f2=%d | put o2=%g |
				window min1=%g max1=%g min2=%g max2=%g |
				math output="input*input" |
				stack axis=2 norm=n |
				stack axis=1 norm=n |
				math output="sqrt(input)"
				'''%(nx,nx,ox,norm[0],norm[1],norm[2],norm[3]))

			Flow(l2mr,[name+'-x-it-%d'%j,name+'-x-it-%d'%(j+1)],
				'''
				add scale=1,-1 ${SOURCES[1]} | window n2=%d | put o2=%g |
				window min1=%g max1=%g min2=%g max2=%g |
				math output="input*input" |
				stack axis=2 norm=n |
				stack axis=1 norm=n |
				math output="sqrt(input)"
				'''%(nx,ox,norm[0],norm[1],norm[2],norm[3]))

			movie_l2m.append(l2m)

			movie_l2mr.append(l2mr)

		Plot('a-'+name+'-movie-l2m',movie_l2m,plotgrl2m(outer,plottitle))

		Plot('a-'+name+'-movie-l2mr',movie_l2mr,plotgrl2m(outer,plottitle))

		#Plot('a-'+name+'-movie-l2m',movie_l2m,plotgrl2m(outer,plottitle,'wanttitle=n'))

		#Plot('a-'+name+'-finalx',name+'-x-it-%d'%(outer-1),plotfinalx(plottitle,display))

		"""
		# comparison with "true" model if we have it
		if truemod != None:

			movie_l2mt = []

			for j in range (outer):

				l2mt = name+'-l2mt-it-%d'%j

				Flow(l2mt,[name+'-x-it-%d'%j,truemod],
				'''
				add scale=1,-1 ${SOURCES[1]} |
				window min1=%g max1=%g min2=%g max2=%g |
				math output="input*input" |
				stack axis=2 norm=n |
				stack axis=1 norm=n |
				math output="sqrt(input)"
				'''%(l2min1,l2max1,l2min2,l2max2))

				movie_l2mt.append(l2mt)

			Result('a_'+name+'_movie_l2mt',movie_l2mt,plotgrl2mt(outer,plottitle))

			Plot('a_'+name+'_movie_l2mt',movie_l2mt,plotgrl2mt(outer,plottitle))
		"""
