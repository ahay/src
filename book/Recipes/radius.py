from rsf.proj import *

def radius(high, low,               # initial high-resolution and legacy images
           niter,                   # number of corrections
           c,                       # 'step length' for radius corrections. Can be type int or float for constant c 
                                    # or type array for changing c.
           bias=-15, clip=30,       # bias and clip for display
           rect1=40, rect2=80, rect3=1,     # radius for local frequency calculation
           maxrad=1000,             # maximum allowed radius
           theor=True,              # use theoretical smoothing radius
           scale=9,                 # scale for theoretical smoothing radius
           initial=10,              # initial value for contant smoothing radius 
           minval=0,
           maxval=25,
           titlehigh="Hires",
           titlelow="Legacy",
           ind = ''):

    if type(c) is float or type(c) is int:
        c = [c]*niter

    def seisplot(name):
        return 'grey title="%s"'%name
    
    locfreq = '''iphase order=10 rect1=%d rect2=%d rect3=%d hertz=y complex=y | 
                 put label="Frequency" unit=Hz'''%(rect1,rect2,rect3)

    def locfreqplot(name):
        return 'grey mean=y color=j scalebar=y title="%s" '%name

    freqdif = 'add scale=-1,1 ${SOURCES[1]} | put label=Frequency'
    def freqdifplot3(num):
        return '''byte gainpanel=all bar=bar.rsf allpos=y | 
                  grey3 frame1=400 frame2=40 frame3=30 allpos=y color=j scalebar=y mean=y 
                  title="Difference in Local Frequenies %s"''' %(num)

    def freqdifplot(num):
        return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %s" 
                  clip=%d bias=%d minval=%d maxval=%d''' %(num,clip,bias,minval,maxval) 

    specplot = '''cat axis=2 ${SOURCES[1]} | 
                  scale axis=1 | window max1=180 |
                  graph title="Normalized Spectra" label2="Amplitude" unit2=""'''

    def rectplot(name):
        return 'grey color=j mean=y title="%s" scalebar=y barlabel=Radius barunit=samples'%name

    smooth = 'nsmooth1 rect=${SOURCES[1]}'

    #Result(high, seisplot(titlehigh))
    #Result(low, seisplot(titlelow))
    def n(name):
            return name+str(ind)
    high_freq = n('high-freq')
    low_freq = n('low-freq')

    Flow(high_freq,high,locfreq)
    Result(high_freq,locfreqplot('%s Local Frequency'%titlehigh))

    Flow(low_freq,low,locfreq)
    Result(low_freq,locfreqplot('%s Local Frequency'%titlelow))

    freqdif_init = n('freqdif-init')

    Flow(freqdif_init,[low_freq, high_freq],freqdif)
    Plot(freqdif_init,freqdifplot3(''))

    # initial smoothing radius
    rect = n('rect0')
    if (theor):
        from math import pi
        Flow(rect,[low_freq, high_freq],'''math f1=${SOURCES[1]} 
              output="sqrt(%g*(1/(input*input)-1/(f1*f1)))/%g" '''%(scale,2*pi*0.001))
    else:
        Flow(rect,low_freq,'math output=%f'%initial)

    Plot(rect,rectplot("Smoothing Radius 0"))

    high_smooth = n('high-smooth0')
    Flow(high_smooth,[high,rect],smooth)
    Result(high_smooth, seisplot("%s Smooth 0"%titlehigh))

    high_spec = n('high-spec')
    low_spec = n('low-spec')
    Flow(high_spec,high,'spectra all=y')
    Flow(low_spec,low,'spectra all=y')

    high_ss = n('high-smooth-spec0')
    nspec = n('nspectra-diff')
    Flow(high_ss,high_smooth,'spectra all=y')
    Result(nspec,[high_spec, low_spec],specplot)
    Result(high_ss,[high_ss, low_spec],specplot)
    
    high_sf = n('high-smooth-freq0')
    Flow(high_sf,high_smooth,locfreq)
    Result(high_sf,locfreqplot("%s Local Frequency Smoothed %d" %(titlehigh,0)))
           

    freqdif_filt = n('freqdif-filt0')
    Flow(freqdif_filt,[low_freq, high_sf],freqdif) 
    #Result('freqdif-filt0',freqdifplot('0')) 
    Result(freqdif_filt,freqdifplot3('0')) 

    for i in range(1, niter+1): 

        j = i-1
        rect = n('rect%i'%i)
        rectj = n('rect%i'%j)
        fdfj = n('freqdif-filt%i'%j)

        Flow(rect,[rectj,fdfj],'radius freq=${SOURCES[1]} c=%f'%c[j])
        Plot(rect,rectplot("Smoothing Radius %d")%i)
        
        hs = n('high-smooth%i'%i)
        hss = n('high-smooth-spec%i'%i)
        hsf = n('high-smooth-freq%i'%i)
        fdf = n('freqdif-filt%i'%i)

        Flow(hs,[high,rect],smooth)
        Plot(hs, seisplot('%s Smooth %d'%(titlehigh,i)))

        Flow(hss,hs,'spectra all=y')
        Plot(hss,[hss,low_spec],specplot)
        
        Flow(hsf,hs,locfreq)
        Plot(hsf,locfreqplot('%s Local Frequency Smoothed %d'%(titlehigh,i)))

        Flow(fdf,[low_freq, hsf],freqdif)
        Plot(fdf,freqdifplot(str(i)))
        
