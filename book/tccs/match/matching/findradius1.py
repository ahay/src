from rsf.proj import *

def findradius1(D1, D2,              # original data (high) and filtered data (low)
            niter=5,                 # number of corrections
            rec1=15, rec2=15,        # smooth division parameters
            recmin=0,                # smooth division minimum radius (for decreasing radius)
            recmax=1000,             # smooth division radius 1st iteration (choose large value to estimate the stationary radius)
            rec1f=15, rec2f=15,      # smooth local frequency parameters
            maxvalf=80,              # maximum value for plotting difference in local frequency
            biasf=40,                # bias for plotting local frequency difference
            theor=False,             # use smoothed version of theoretical smoothing radius as starting model 
            rectheor=1000,           # smoothing radius for theoretical radius 
            scale=9,                 # scale for theoretical smoothing radius
            d1=0.001,                # interval spacing in 1st dimension (time)
            initial=10,              # initial value for constant smoothing radius 
            titlehigh="D1",          # original data (high) title
            titlelow="D2",           # filtered data (low) title 
            localfreq=False,         # estimate based on local frequency attribute 
            it=0):                   # correction number 

    '''
    Estimating the vertical non-stationary triangle smoothing radius
    with least-squares iterative inversion 
    '''

    # local frequency calculation 
    locfreq = '''iphase order=10 rect1=%d rect2=%d hertz=y complex=y | 
                 put label="Frequency" unit=Hz'''%(rec1f,rec2f)

    # abs difference in local frequency calculation 
    freqdif = 'add scale=-1,1 ${SOURCES[1]} | math output="abs(input)" | put label=Frequency'
    
    # Plotting   
    def locfreqplot(name):
        return 'grey minval=10 maxval=80 bias=45 color=j scalebar=y title="%s" '%(name)

    def specplot(name):
        return '''cat axis=2 ${SOURCES[1]} | 
                  scale axis=1 | window max1=120 | 
                  graph title="%s" label2="Amplitude" unit2="" '''%(name)

    def freqdifplot(name):
        return 'grey color=j scalebar=y maxval=%d bias=%d title="%s"' %(maxvalf,biasf,name)

    def seisplot(name):
        return 'grey color=seismic title="%s" '%name

    def rectplot(name):
        return '''
        grey color=viridis mean=y title="%s" scalebar=y 
        barlabel=Radius barunit=samples 
        '''%name

    # get spectral content of D1 and D2 
    Flow('%s%i-spec'%(D1,it),D1,'spectra all=y')
    Flow('%s%i-spec'%(D2,it),D2,'spectra all=y')
    Result('nspectra%i'%it,'%s%i-spec %s%i-spec'%(D1,it,D2,it),specplot("Normalized Spectra"))
    
    # local frequency of D1 
    Flow('%s%i-locfreq'%(D1,it),D1,locfreq)
    Result('%s%i-locfreq'%(D1,it),locfreqplot('%s local frequency'%titlehigh))

    # local frequency of D2 
    Flow('%s%i-locfreq'%(D2,it),D2,locfreq)
    Result('%s%i-locfreq'%(D2,it),locfreqplot('%s local frequency'%titlelow))
    
    # initial difference in local frequencies 
    Flow('locfreqdif%i'%(it),'%s%i-locfreq %s%i-locfreq'%(D1,it,D2,it),freqdif)
    Result('locfreqdif%i'%(it),freqdifplot("initial difference in local frequency"))

    # initial smoothing radius
    if (theor):
        from math import pi
        Flow('rect0%i'%it,'%s%i-locfreq %s%i-locfreq'%(D2,it,D1,it),
            '''
            math f1=${SOURCES[1]} output="sqrt(%g*(1/(input*input)-1/(f1*f1)))/%g" |
            smooth rect1=%g rect2=%g 
              '''%(scale,2*pi*d1,rectheor,rectheor))
    else:
        Flow('rect0%i'%it,D1,'math output="input*0+%g" '%(initial))

    Result('rect%d%i'%(0,it),rectplot("smoothing radius %d"%(0)))
    Plot('rect%d%i'%(0,it),rectplot("smoothing radius %d"%(0)))
    
    # initial smoothed D1
    smoothed = '%s-smooth0%i'%(D1,it)
    Flow(smoothed,[D1,'rect0%i'%(it)],'nsmooth rect1=${SOURCES[1]} ')
    Result(smoothed,seisplot("%s smoothed %d"%(D1,0))) 
    
    # initial smoothing derivative applied to D1
    der = '%s-smoothder0%i'%(D1,it)
    Flow(der,[D1,'rect0%i'%(it)],'ndsmooth rect1=${SOURCES[1]}')
    
    # local frequency 
    Flow('%s-locfreq'%(smoothed),smoothed,locfreq)
    
    # difference in local frequency 
    Flow('locfreqdif%d%i'%(0,it),'%s-locfreq %s%i-locfreq'%(smoothed,D2,it),freqdif) 
    Plot('locfreqdif%d%i'%(0,it),freqdifplot("difference in local frequency %d"%(0))) 

    # set initial smooth division parameters 
    rect1=recmax
    rect2=recmax

    for i in range(1, niter+1):
        j = i-1

        # calculate new radius
        Flow('rect%d%i'%(i,it),[D2,smoothed,der,'rect%d%i'%(j,it)],
            '''
            add ${SOURCES[1]} scale=1,-1 |
            divn den=${SOURCES[2]} rect1=%g rect2=%g|
            add ${SOURCES[3]} scale=1,1
            '''%(rect1,rect2))
        if i==5:
        	Result('rectn%d%i'%(i,it),'rect%d%i'%(i,it),rectplot("smoothing radius %d"%i))
        else:
        	Result('rect%d%i'%(i,it),rectplot("smoothing radius %d"%i))
        
        # smoothing D1
        smoothed = '%s-smooth%d%i'%(D1,i,it)
        Flow(smoothed,[D1,'rect%d%i'%(i,it)],'nsmooth rect1=${SOURCES[1]}')
        Result(smoothed,seisplot("%s smoothed %d"%(D1,i)))
        
        # smoothing derivative applied to D1 
        der = '%s-smoothder%d%i'%(D1,i,it)
        Flow(der,[D1,'rect%d%i'%(i,it)],'ndsmooth rect1=${SOURCES[1]}')
        
        # normalized spectra after smoothing 
        Flow('%s-spec'%(smoothed),smoothed,'spectra all=y')
        Result('%s-spec'%(smoothed),'%s-spec %s%i-spec'%(smoothed,D2,it),specplot("Normalized Spectra %d"%(i)))
        
        # local frequency of smoothed D1 
        Flow('%s-locfreq'%(smoothed),smoothed,locfreq)
        Result('%s-locfreq'%(smoothed),locfreqplot("%s smoothed local frequency %d" %(titlehigh,i)))
        
        # abs difference in local frequency after smoothing  
        Flow('locfreqdif%d%i'%(i,it),'%s-locfreq %s%i-locfreq'%(smoothed,D2,it),freqdif) 
        Result('locfreqdif%d%i'%(i,it),freqdifplot("difference in local frequency %d"%(i))) 

        # set new smooth division parameters 
        rect1=rec1/(2**j)
        rect2=rec1/(2**j)

        if rect1 < recmin:
            rect1=recmin
        if rect2 < recmin:
            rect2=recmin






