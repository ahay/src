from rsf.proj import *

def find_radius(high, low,           # initial high-resolution and legacy images
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
            titlehigh="Hires",       # high title
            titlelow="Legacy",       # low title
            it=0):                   # correction number (from newarp.py)

    if type(c) is float or type(c) is int:
        c = [c]*niter

    def seisplot(name):
        return 'grey color=rwb scalebar=y title="%s" '%name
    
    locfreq = '''iphase order=10 rect1=%d rect2=%d hertz=y complex=y | 
                 put label="Frequency" unit=Hz'''%(rect1,rect2)

    def locfreqplot(name):
        return 'grey mean=y color=j scalebar=y title="%s" '%name

    freqdif = 'add scale=-1,1 ${SOURCES[1]} | put label=Frequency'

    def freqdifplot(num):
        #return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" ''' %num
        #return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" 
        #          clip=%d bias=%d minval=%d maxval=%d''' %(num,clip,bias,minval,maxval) 
        return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" 
                  clip=%d bias=%d minval=%d maxval=%d screenratio=1.2''' %(num,22,-15,-15,8) 

    specplot = '''cat axis=2 ${SOURCES[1]} | 
                  scale axis=1 | window max1=125 |
                  graph title="Normalized Spectra" label2="Amplitude" unit2=""'''

    def rectplot(name):
        return 'grey color=viridis mean=y title="%s" scalebar=y barlabel=Radius barunit=samples '%name

    smooth = 'nsmooth1 rect=${SOURCES[1]}'

    Flow('high-freq%i'%it,high,locfreq)
    Result('high-freq%i'%it,locfreqplot('%s Local Frequency'%titlehigh))

    Flow('low-freq%i'%it,low,locfreq)
    Result('low-freq%i'%it,locfreqplot('%s Local Frequency'%titlelow))

    Flow('freqdif%i'%it,'low-freq%i high-freq%i'%(it,it),freqdif)
    Result('freqdif%i'%it,freqdifplot(-1))


    # initial guess - gauss-newton
    #rec1=10000
    #rec2=10000
    rec1 = rec2 = 30 # 30 

    Flow('new_rect01','rect10','math output="input*0+5" ')
    Plot('new_rect01',' grey label2="Trace" label1="Time" unit1=s unit2="" labelsz=10 title="" wherexlabel=b wheretitle=t screenratio=1.2 labelfat=4 font=2  mean=y scalebar=y barlabel=Radius barunit=samples color=j scalebar=y')
    
    for i in range(1, niter+1): 
        j = i-1
        # Gauss-Newton Method 
        # D1 = high
        # D2 = low
        # Ri+1 = Ri + (D2 - S(Ri)*D1)/(S'(Ri)*D1)

        smoothed = 'smooth_%d%i'%(j,it) 
        der = "der_%d%i"%(j,it) 
        Flow(smoothed,'%s new_rect%d%i'%(high,j,it),smooth)
        Flow(der,'%s new_rect%d%i'%(high,j,it),'ndsmooth rect1=${SOURCES[1]} ider=1')

        Flow('new_rect%d%i'%(i,it),[low,smoothed,der,'new_rect%d%i'%(j,it)],
            '''
            add ${SOURCES[1]} scale=1,-1 |
            divn den=${SOURCES[2]} rect1=%g rect2=%g|
            add ${SOURCES[3]} scale=1,1
            '''%(rec1,rec2))
        Plot('new_rect%d%i'%(i,it),rectplot("Smoothing Radius %d, GN")%i)

        Flow('new_high-smooth%d%i'%(i,it),'%s new_rect%d%i'%(high,i,it),smooth)
        Plot('new_high-smooth%d%i'%(i,it), seisplot('%s Smooth %d'%(titlehigh,i)))

        #Flow('new_high-smooth-spec%d%i'%(i,it),'new_high-smooth%d%i'%(i,it),'spectra all=y')
        #Plot('new_high-smooth-spec%d%i'%(i,it),'new_high-smooth-spec%d%i low-spec%i'%(i,it,it),specplot)
        
        #Flow('new_high-smooth-freq%d%i'%(i,it),'new_high-smooth%d%i'%(i,it),locfreq)
        #Result('new_high-smooth-freq%d%i'%(i,it),locfreqplot('%s Local Frequency Smoothed %d'%(titlehigh,i)))

        #Flow('new_freqdif-filt%d%i'%(i,it),'low-freq%i new_high-smooth-freq%d%i'%(it,i,it),freqdif)
        #Plot('new_freqdif-filt%d%i'%(i,it),freqdifplot(i))


        rec1=30/(2**(j)) # 30
        rec2 = 30/(2**(j)) # 30 

        rec_min1 = rec_min2 = 5

        if rec1 < rec_min1:
          rec1 = rec_min1 
        if rec2 < rec_min2:
          rec2 = rec_min2 








