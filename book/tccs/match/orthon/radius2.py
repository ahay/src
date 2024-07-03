from rsf.proj import *

def radius2(high, low,               # initial high-resolution and low-resolution images
            niter,                   # number of corrections
            c,                       # step length for radius corrections. Can be type int or float for constant c 
                                     #      or type array for changing c.
            bias=-15, clip=30,       # bias and clip for display
            rect1=40, rect2=80,      # radius for local frequency calculation
            maxrad=1000,             # maximum allowed radius
            theor=True,              # use theoretical smoothing radius
            scale=9,                 # scale for theoretical smoothing radius
            initial=10,              # initial value for contant smoothing radius 
            minval=0, maxval=25,     # minval and maxval for freqdif-filt display
            titlehigh="High",        # high title
            titlelow="Low",          # low title
            it=0):                   # correction number (from newarp.py)

    if type(c) is float or type(c) is int:
        c = [c]*niter

    def seisplot(name):
        return 'grey color=g title="%s" screenratio=1.2'%name
    
    locfreq = '''iphase order=10 rect1=%d rect2=%d hertz=y complex=y | 
                 put label="Frequency" unit=Hz'''%(rect1,rect2)

    def locfreqplot(name):
        return 'grey mean=y color=j scalebar=y title="%s" screenratio=1.2'%name

    freqdif = 'add scale=-1,1 ${SOURCES[1]} | put label=Frequency'

    def freqdifplot(num):
        return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" 
                  screenratio=1.2''' %(num) 

    specplot = '''cat axis=2 ${SOURCES[1]} | 
                  scale axis=1 | window max1=125 |
                  graph title="Normalized Spectra" label2="Amplitude" unit2=""'''

    def rectplot(name):
        return '''
        grey color=viridis title="%s" scalebar=y 
        barlabel=Radius barunit=samples screenratio=1.2
        bias=4.3 minval=2 maxval=7 
        '''%name
    
    def rectplot2(name):
        return '''
        grey color=viridis title="%s" scalebar=y 
        barlabel=Radius barunit=samples screenratio=1.2
        mean=y minval=2 maxval=7 
        '''%name
    
    smooth = 'nsmooth rect1=${SOURCES[1]}'

    # local frequency of high-resolution image
    Flow('high-freq%i'%it,high,locfreq)
    #Result('high-freq%i'%it,locfreqplot('%s Local Frequency'%titlehigh))

    # local frequency of low-resolution image
    Flow('low-freq%i'%it,low,locfreq)
    #Result('low-freq%i'%it,locfreqplot('%s Local Frequency'%titlelow))

    # initial difference in local freqeuncy 
    Flow('freqdif%i'%it,'low-freq%i high-freq%i'%(it,it),freqdif)
    #Result('freqdif%i'%it,freqdifplot(-1))

    # initial guess for the smoothing radius - Line-Search method
    if (theor):
        from math import pi
        Flow('rect0%i'%it,'low-freq%i high-freq%i'%(it,it),'''math f1=${SOURCES[1]} 
              output="sqrt(%g*(1/(input*input)-1/(f1*f1)))/%g" '''%(scale,2*pi*0.001))
    else:
        Flow('rect0%i'%it,'low-freq%i'%it,'math output=%f'%initial)
    Result('rect0%i'%it,rectplot("Smoothing Radius 0"))

    # initial smoothing of high-frequency data 
    Flow('high-smooth0%i'%it,'%s rect0%i'% (high,it),smooth)
    Result('high-smooth0%i'%it, seisplot("%s Smooth 0"%titlehigh))

    # local frequency of initial smoothing of high-frequency data 
    Flow('high-smooth-freq0%i'%it,'high-smooth0%i'%it,locfreq)
    #Result('high-smooth-freq0%i'%it,locfreqplot("%s Local Frequency Smoothed %d" %(titlehigh,0)))

    # difference in local frequency 
    Flow('freqdif-filt0%i'%it,'low-freq%i high-smooth-freq0%i'%(it,it),freqdif) 
    #Result('freqdif-filt0%i'%it,freqdifplot(0)) 

    # initial guess for the smoothing radius - Gauss-Newton method
    Flow('new_rect00','rect10','math output="input*0+4" ')
    Result('new_rect00',' grey label2="Trace" label1="Time" unit1=s unit2="" labelsz=10 title="" wherexlabel=b wheretitle=t screenratio=1.2 labelfat=4 font=2  mean=y scalebar=y barlabel=Radius barunit=samples color=j scalebar=y')
    
    # smooth divison parameters 
    rec1=10
    rec2=10

    prog=Program('radius2.c')
    for i in range(1, niter+1): 
        j = i-1

        ################################
        # Line-Search Radius Estimation 
        ################################
        Flow('rect%d%i rect-low%d%i rect-high%d%i'%(i,it,i,it,i,it),'rect%d%i freqdif-filt%d%i %s'%(j,it,j,it,prog[0]),
              './${SOURCES[2]} freq=${SOURCES[1]} low=${TARGETS[1]} high=${TARGETS[2]} c=%f'%(c[j]))
        Result('rect%d%i'%(i,it),rectplot(""))
        
        Flow('high-smooth%d%i'%(i,it),'%s rect-high%d%i'%(high,i,it),smooth)
        #Result('high-smooth%d%i'%(i,it), seisplot('%s Smooth %d'%(titlehigh,i)))

        Flow('low-smooth%d%i'%(i,it),'%s rect-low%d%i'%(low,i,it),smooth)
        #Result('low-smooth%d%i'%(i,it), seisplot('%s Smooth %d'%(titlelow,i)))

        Flow('high-smooth-freq%d%i'%(i,it),'high-smooth%d%i'%(i,it),locfreq)
        #Result('high-smooth-freq%d%i'%(i,it),locfreqplot('%s Local Frequency Smoothed %d'%(titlehigh,i)))

        Flow('low-smooth-freq%d%i'%(i,it),'low-smooth%d%i'%(i,it),locfreq)
        #Result('low-smooth-freq%d%i'%(i,it),locfreqplot('%s Local Frequency Smoothed %d'%(titlelow,i)))

        Flow('freqdif-filt%d%i'%(i,it),'low-smooth-freq%d%i high-smooth-freq%d%i'%(i,it,i,it),freqdif)
        #Result('freqdif-filt%d%i'%(i,it),freqdifplot(i))
       

        ################################
        # Gauss-Newton Radius Estimation  
        ################################
        smoothed = 'smooth_%d'%(j) 
        der = "der_%d"%(j) 
        Flow(smoothed,[high, 'new_rect%d%i'%(j,it)],smooth)
        Flow(der,[high, 'new_rect%d%i'%(j,it)],'ndsmooth rect1=${SOURCES[1]} ider=1')
        Flow('new_rect%d%i'%(i,it),[low,smoothed,der,'new_rect%d%i'%(j,it)],
            '''
            add ${SOURCES[1]} scale=1,-1 |
            divn den=${SOURCES[2]} rect1=%g rect2=%g|
            add ${SOURCES[3]} scale=1,1
            '''%(rec1,rec2))
        Result('new_rect%d%i'%(i,it),rectplot2(""))


