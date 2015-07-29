## This is a script for the test of deblending using sparsity promotion with shaping regularization (Using the SNR definition from Araz).

## Copyright (C) 2013
## Yangkang Chen, The Univeristy of Texas at Austin

## The are several assumptions for implementing the following scripts:
## 1. data1 and data2 have the same dimensions 
## 2. data1 and data2 are sparse in fourier or seislet domains
## 3. the unblended data is known (in other words it's used for testing). 
##    Otherwise, the unblended1(2) is substituted by blended1(2), in which case
##    The error diagram has no physical meanings.



from rsf.proj import*
def Grey(data,other): 
	Result(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" wherexlabel=b wheretitle=t color=b %s'%other)

def Greyplot(data,other): 
	Plot(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" wherexlabel=b wheretitle=t color=b %s'%other)

def Graph(data,other):
	Result(data,'graph label1="" label2="" unit1="" unit2=""  title="" wherexlabel=b wheretitle=t %s' %other)

def Graphplot(data,other):
	Plot(data,'graph label1="" label2="" unit1="" unit2=""  title="" wherexlabel=b wheretitle=t %s' %other)

def shapeslet(sig1,sig2,dip1,dip2,mute1,mute2,padno,n2,mode,thr,i): #mute1(2) is a string, when there is no mute, mute1(2) is "cp"
    Flow(sig1+"p%g"%(i+1),[sig1,dip1],
	 '''
	 %s |pad n2=%g  |
         seislet dip=${SOURCES[1]} eps=0.1 
		  adj=y inv=y unit=y type=b |
         threshold1 type=%s ifperc=1 thr=%f |
         seislet dip=${SOURCES[1]} eps=0.1  
			inv=y unit=y type=b | 
         window n2=%g |%s
         '''%(mute1,padno,mode,thr,n2,mute1))
    Flow(sig2+"p%g"%(i+1),[sig2,dip2],
	 ''' 
	 %s|pad n2=%g  |
         seislet dip=${SOURCES[1]} eps=0.1 
		  adj=y inv=y unit=y type=b |
         threshold1 type=%s ifperc=1 thr=%f |
         seislet dip=${SOURCES[1]} eps=0.1  
			inv=y unit=y type=b | 
         window n2=%g | %s
         '''%(mute2,padno,mode,thr,n2,mute2))
    return (sig1+"p%g"%(i+1),sig2+"p%g"%(i+1))

def shapefft(sig1,sig2,mute1,mute2,mode,thr,i):#mute1(2) is a string, when there is no mute, mute1(2) is "cp"
    Flow(sig1+"p%g"%(i+1),sig1,
	 '''
	 %s|fft1 | fft3 axis=2 pad=1|
         threshold1 type=%s ifperc=1 thr=%f |
         fft3 inv=y axis=2 | fft1 inv=y |%s
         '''%(mute1,mode,thr,mute1))
    Flow(sig2+"p%g"%(i+1),sig2,
	 '''
	 %s|fft1 | fft3 axis=2 pad=1|
         threshold1 type=%s ifperc=1 thr=%f |
         fft3 inv=y axis=2 | fft1 inv=y | %s
         '''%(mute2,mode,thr,mute2))
    return (sig1+"p%g"%(i+1),sig2+"p%g"%(i+1))

def shapefxdecon(sig1,sig2,mute1,mute2,n2,i): #mute1(2) is a string, when there is no mute, mute1(2) is "cp"
    Flow(sig1+"p%g"%(i+1),[sig1],
	 '''
	 %s | fxdecon n2w=%g 
         |%s
         '''%(mute1,n2,mute1))
    Flow(sig2+"p%g"%(i+1),[sig2],
	 ''' 
 	 %s | fxdecon n2w=%g 
         |%s
         '''%(mute2,n2,mute2))
    return (sig1+"p%g"%(i+1),sig2+"p%g"%(i+1))

def step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,n2,fraction):
#    Flow(init1+'sigtemp'+'%g'%(i+1),[sig1,sig2],'cat axis=2 ${SOURCES[1]}')
#    Flow(init1+'shottimetemp'+'%g'%(i+1),[shottime1,shottime2],'cat axis=2 ${SOURCES[1]}')
#    Flow(init1+'blendedtemp'+'%g'%(i+1),[init1+'sigtemp'+'%g'%(i+1),init1+'shottimetemp'+'%g'%(i+1)],
#    	'''
#    	blend shot_time_in=${SOURCES[1]} shot_time_out=${SOURCES[1]} 
#    	''')
#    	
#    Flow(init1+"%g"%(i+1),[init1+'blendedtemp'+'%g'%(i+1),blended1,sig1], 
#         '''
#         window n2=%g | add scale=-1,1 ${SOURCES[1]} | add scale=%g,%g ${SOURCES[2]}
#         '''%(n2,fraction,1))
#    Flow(init2+"%g"%(i+1),[init1+'blendedtemp'+'%g'%(i+1),blended2,sig2], 
#         '''
#         window n2=%g f2=%g | add scale=-1,1 ${SOURCES[1]} | add scale=%g,%g ${SOURCES[2]}
#         '''%(n2,n2-1,fraction,1))
    
    Flow(init1+"%g"%(i+1),[sig2,shottime1,shottime2,blended1,sig1], 
         '''
          blend shot_time_in=${SOURCES[2]} shot_time_out=${SOURCES[1]}
	    | add scale=%f,%f,%f ${SOURCES[3]} ${SOURCES[4]}
         '''%(-fraction,fraction,1-fraction))

    Flow(init2+"%g"%(i+1),[sig1,shottime2,shottime1,blended2,sig2],
         '''
          blend shot_time_in=${SOURCES[2]} shot_time_out=${SOURCES[1]}
	    | add scale=%f,%f,%f ${SOURCES[3]} ${SOURCES[4]}
         '''% (-fraction,fraction,1-fraction))
    return (init1+"%g"%(i+1),init2+"%g"%(i+1))    
    

def deblendfft(unblended1,
	  unblended2,
	  blended1,
          blended2,
	  init1,
	  init2,
	  deblended1,
	  deblended2,
	  shottime1,
	  shottime2,
	  mute1,
	  mute2,
	  n1,
	  n2,
	  niter,
	  mode,
	  thr,
	  clip,
	  fraction):

    sig1, sig2 = (init1,init2)
    sigs1=[]
    sigs2=[]
    snrsa=[]
    snrsb=[]
    Flow('zero',unblended1,'math output=0')
    Flow('energy1',[unblended1, 'zero'],'diff match=${SOURCES[1]}')
    Flow('energy2',[unblended2, 'zero'],'diff match=${SOURCES[1]}')
    for i in range(niter):
 	old1,old2 = (sig1,sig2)
	snra=init1+'-snra%g'%i
	snrb=init2+'-snrb%g'%i

	Flow(snra,[unblended1,sig1,'energy1'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	Flow(snrb,[unblended2,sig2,'energy2'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	
	(nsig1,nsig2)=step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,n2,fraction)
	(sig1,sig2)=shapefft(nsig1,nsig2,mute1,mute2,mode,thr,i)

	
    	Greyplot(sig1,'title="Esti R 1 (iter=%g)" clip=%g'% (i+1,clip))
    	Greyplot(sig2,'title="Esti R 2 (iter=%g)" clip=%g' % (i+1,clip))	

	sigs1.append(sig1)
    	sigs2.append(sig2)
        snrsa.append(snra)
        snrsb.append(snrb)

    Flow(init1+'-snrsa',snrsa,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp | put n1=%g o1=0 d1=1 '%(len(snrsa),len(snrsa)))
    Flow(init2+'-snrsb',snrsb,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp |  put n1=%g o1=0 d1=1'%(len(snrsb),len(snrsb)))

    Graph(init1+'-snrsa','title="Data error 1" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)	
    Graph(init2+'-snrsb','title="Data error 2" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)		

#Making movie
    Plot(init1+'-sigs1',sigs1,'Movie')
    Plot(init2+'-sigs2',sigs2,'Movie')
    Flow(deblended1,sig1,'cp')
    Flow(deblended2,sig2,'cp')

def deblendslet(unblended1,
	  unblended2,
	  blended1,
          blended2,
	  init1,
	  init2,
	  deblended1,
	  deblended2,
	  shottime1,
	  shottime2,
	  mute1,
	  mute2,
	  n1,
	  n2,
	  r1,
	  r2,
	  padno,
	  niter,
	  ddip,
	  mode,
	  thr,
	  clip,
	  fhi,
	  fraction):
    Flow('mask',blended1,'math output=1 | pad n2=%g'%(padno))
    sig1, sig2 = (init1,init2)
    sigs1=[]
    sigs2=[]
    sigs=[]
    snrsa=[]
    snrsb=[]
    Flow('zero',unblended1,'math output=0')
    Flow('energy1',[unblended1, 'zero'],'diff match=${SOURCES[1]}')
    Flow('energy2',[unblended2, 'zero'],'diff match=${SOURCES[1]}')
    for i in range(niter):
#######################################################################################
# update dip map every "ddip" iterations
#######################################################################################
    	if i % ddip ==0 :
		dip1=init1+'dip%g'%int(i/ddip)
		dip2=init2+'udip%g'%int(i/ddip)
		Flow(dip1,[sig1,'mask'],
     			'''
			bandpass fhi=%g | pad n2=%g | 
			dip mask=${SOURCES[1]} rect1=%g rect2=%g
			'''%(fhi,padno,r1,r2))
		Flow(dip2,[sig2,'mask'],
     			'''
			bandpass fhi=%g | pad n2=%g | 
			dip mask=${SOURCES[1]} rect1=%g rect2=%g
			'''%(fhi,padno,r1,r2))
#######################################################################################
 	old1,old2 = (sig1,sig2)
	snra=init1+'-snra%g'%i
	snrb=init2+'-snrb%g'%i

	Flow(snra,[unblended1,sig1,'energy1'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	Flow(snrb,[unblended2,sig2,'energy2'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	
	(nsig1,nsig2)=step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,n2,fraction)
	(sig1,sig2)=shapeslet(nsig1,nsig2,dip1,dip2,mute1,mute2,padno,n2,mode,thr,i)

    	Greyplot(sig1,'title="Esti R 1 (iter=%g)" clip=%g'% (i+1,clip))
    	Greyplot(sig2,'title="Esti R 2 (iter=%g)" clip=%g' % (i+1,clip))	

	sigs1.append(sig1)
    	sigs2.append(sig2)
	
	Plot(sig1+'t',[sig1,sig2],'SideBySideAniso')
	sigs.append(sig1+'t')
	
        snrsa.append(snra)
        snrsb.append(snrb)

    Flow(init1+'-snrsa',snrsa,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp | put n1=%g o1=0 d1=1 '%(len(snrsa),len(snrsb)))
    Flow(init2+'-snrsb',snrsb,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp |  put n1=%g o1=0 d1=1'%(len(snrsb),len(snrsb)))

    Graph(init1+'-snrsa','title="Data error 1" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)	
    Graph(init2+'-snrsb','title="Data error 2" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)	

#Making movie
    Plot(init1+'-sigs1',sigs1,'Movie')
    Plot(init2+'-sigs2',sigs2,'Movie')

    Plot(init1+'-sigs',sigs,'Movie')

    Flow(deblended1,sig1,'cp')
    Flow(deblended2,sig2,'cp')
   

def deblendfxdecon(unblended1,
	  unblended2,
	  blended1,
          blended2,
	  init1,
	  init2,
	  deblended1,
	  deblended2,
	  shottime1,
	  shottime2,
	  mute1,
	  mute2,
	  n1,
	  n2,
	  niter,
	  clip,
	  fraction):
    sig1, sig2 = (init1,init2)
    sigs1=[]
    sigs2=[]
    snrsa=[]
    snrsb=[]
    Flow('zero',unblended1,'math output=0')
    Flow('energy1',[unblended1, 'zero'],'diff match=${SOURCES[1]}')
    Flow('energy2',[unblended2, 'zero'],'diff match=${SOURCES[1]}')
    for i in range(niter):
#######################################################################################
 	old1,old2 = (sig1,sig2)
	snra=init1+'-snra%g'%i
	snrb=init2+'-snrb%g'%i

	Flow(snra,[unblended1,sig1,'energy1'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	Flow(snrb,[unblended2,sig2,'energy2'],'diff match=${SOURCES[1]} | math a=${SOURCES[2]} output="a/input"' )
	
	(nsig1,nsig2)=step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,n2,fraction)
	(sig1,sig2)=shapefxdecon(nsig1,nsig2,mute1,mute2,n2,i)


    	Greyplot(sig1,'title="Esti R 1 (iter=%g)" clip=%g'% (i+1,clip))
    	Greyplot(sig2,'title="Esti R 2 (iter=%g)" clip=%g' % (i+1,clip))	

	sigs1.append(sig1)
    	sigs2.append(sig2)
        snrsa.append(snra)
        snrsb.append(snrb)

    Flow(init1+'-snrsa',snrsa,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp | put n1=%g o1=0 d1=1 '%(len(snrsa),len(snrsb)))
    Flow(init2+'-snrsb',snrsb,'cat axis=2 ${SOURCES[1:%g]}|math output="10*log(input)/log(10)" | transp |  put n1=%g o1=0 d1=1'%(len(snrsb),len(snrsb)))

    Graph(init1+'-snrsa','title="Data error 1" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)	
    Graph(init2+'-snrsb','title="Data error 2" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%g d1=1'%niter)	

#Making movie
    Plot(init1+'-sigs1',sigs1,'Movie')
    Plot(init2+'-sigs2',sigs2,'Movie')

    Flow(deblended1,sig1,'cp')
    Flow(deblended2,sig2,'cp')


