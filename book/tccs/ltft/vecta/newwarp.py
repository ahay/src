from rsf.proj import *
import math, string, sys
import rsf.recipes.version as version

warp0 = '''
warp1 other=${SOURCES[1]} warpin=${SOURCES[2]}
verb=1 nliter=0
'''

getmask = 'add scale=1,-1 ${SOURCES[1]} | mask min=0 | dd type=float'

def psrect(rect):
    return '''
    math min=${SOURCES[2]} max=${SOURCES[1]}
    output="sqrt(1+%d*(1/min^2-1/max^2)*input)" | dd type=int
    ''' % rect

def pprect(rect):
    return '''
    math min=${SOURCES[2]} max=${SOURCES[1]}
    output="sqrt(1+%d*(1/max^2-1/min^2)*(1-input))" | dd type=int
    ''' % rect

balance = '''
nsmooth1 rect=${SOURCES[1]} |
abalance rect1=100 order=100 other=${SOURCES[2]}
'''

def simil(rect1=1,rect2=1):
    return '''
    similarity other=${SOURCES[1]} rect1=%d rect2=%d
    ''' % (rect1,rect2)

def warping(niter,rect1=1,rect2=1,rect3=1):
    return '''
    warp1 other=${SOURCES[1]} warpin=${SOURCES[2]}
    warpout=www.rsf
    verb=1 nliter=%d noamp=1 rect1=%d rect2=%d rect3=%d > ${TARGETS[1]} &&
    warpadd < ${SOURCES[3]} add=www.rsf > $TARGET &&
    rm www.rsf
    ''' % (niter,rect1,rect2,rect3)

def pick(min2,max2,rect1=1,rect2=1,rect3=1,an=0.5):
    return '''
    window min2=%g max2=%g | 
    pick rect1=%d rect2=%d rect3=%d an=%g |
    window''' % (min2,max2,rect1,rect2,rect3,an)

def warpscan(ng,g0,gmax,rect1=1,rect2=1,rect3=1,rect4=1):
    dg = (gmax-g0)/(ng-1)
    return '''
    warpscan other=${SOURCES[1]} niter=100
    ng=%d dg=%g g0=%g rect1=%d rect2=%d rect3=%d rect4=%d |
    math output='(1+input)^4' |
    window''' % (ng,dg,g0,rect1,rect2,rect3,rect4)

def warp2gamma(ss):
    return '''
    math output="input+x1" |
    smoothder %s
    ''' % ('| math output="2*input-1" ','')[ss]

def warp2egamma(ss):
    return '''
    math output="(input+x1)/x1" %s 
    ''' % ('| math output="2*input-1" ','')[ss]

def nwarp2(name,       # name prefix
           pp,ps,      # PP and PS images
           warp,       # initial warp
           nx,         # number of traces
           tmax,       # maximum time for display
           tmin=0,     # minimum time for display
           j2=1,       # trace sumsampling
           trace=None, # seleted trace
           o2=0,       # trace start
           gmin=1,     # minimum gamma
           gmax=4,     # maximum gamma
           niter=20,   # warping iterations
           dt=0.004,   # time sampling
           fmin=0,     # minimum frequency
           fmax=40,    # maximum frequency
           frect=12,   # frequency smoothing
           frame1=10,  # time frame          
           ng=101,     # number of gammas
           g0=0.96,    # first gamma
           pmin=0,     # minimum gamma for picking
           pmax=2,     # maximum gamma for picking
           an=0.5,     # anisotropy for picking
           rect1=50,   # vertical smoothing
           rect2=50,   # lateral smoothing
           iter=2,     # number of iterations
           ss=0,       # PP-PS (0) or PP-SS (1)
           inter=1,    # interleaving
           clip=6      # display clip
           ):
    

    if version.old_version():
        return # think how to do it better

    interg = '''
    pad n2=%d | put n2=%d n3=%d | stack
    ''' % ((nx/inter+1)*inter,inter,nx/inter+1)
    inter = 2*inter    
    interw = '''
    pad n2=%d | put n2=%d n3=%d | stack
    ''' % ((nx/inter+1)*inter,inter,nx/inter+1)
      
    if trace:
        for case in (pp,ps,warp):
            Flow(case+'1',case,'window n2=1 min2=%d' % trace)
        warp1(name+'t',pp+'1',ps+'1',warp+'1',tmax,tmin,
              gmin,gmax,niter,
              dt,fmin,fmax,frect,ng,g0,pmin,pmax,an,rect1,iter,ss)
    else:
        trace=10

    def plot(title):
        return '''
        window min1=%g max1=%g |
        grey title="%s" label1=Time unit1=s clip=%g
        labelfat=3 font=2 titlefat=3 label2=Trace unit2=
        ''' % (tmin,tmax,title,clip)

    vplot = plot('Vp/Vs') + '''
    clip=%g scalebar=y color=j bias=%g minval=%g maxval=%g
    ''' % (0.5*(gmax-gmin),0.5*(gmin+gmax),gmin,gmax)

    balance = '''
    nsmooth1 rect=${SOURCES[1]} |
    abalance rect1=%d rect2=%d order=100 other=${SOURCES[2]}
    ''' % (rect1,rect2)

    ifreq = 'iphase rect1=%d rect2=%d order=100' % (2*rect1,2*rect2)

    def freqplot(title):
        return '''
        scale scale dscale=%g |
        %s clip=%g bias=%g color=j scalebar=y barlabel="Frequency (Hz)"
        ''' % (0.5/(math.pi*dt),plot(title),(fmax-fmin)*0.25,(fmax+fmin)*0.5)

    def specplot(title):
        return '''
        cat axis=2 ${SOURCES[1]} |
        graph title="%s" max1=%g label1="Frequency (Hz)"
        dash=0,1 plotfat=7 label2= 
        ''' % (title,4*fmax)

    def giplot(title):
        return '''
        interleave axis=2 ${SOURCES[1]} |
        window min1=%g max1=%g | scale axis=1 |
        grey title="Interleaved (%s)" label1=Time unit1=s
        label2="Inline" unit2= labelfat=3 font=2 titlefat=3
        screenratio=0.6 screenht=8
        ''' % (tmin,tmax,title)

    def wiplot(title):
        return '''
        interleave axis=2 ${SOURCES[1]} |
        window min1=%g max1=%g | scale axis=1 |
        wiggle poly=y transp=y yreverse=y
        title="Interleaved (%s)"
        label1=Time unit1=s label2="In-line"
        ''' % (tmin,tmax,title)
        
    Plot(pp,plot('PP'))
    Flow(pp+'i',pp,ifreq)
    Plot(pp+'i',freqplot('PP Local Frequency'))

    Result(pp+'line',pp,'Overlay')

    PS = ('PS','SS')[ss]

    Plot(ps,
         '''
         window min1=%g max1=%g |
         grey title="%s" label1=Time unit1=s 
         labelfat=3 font=2 titlefat=3 label2=Trace unit2=
         ''' % (tmin,tmax*2.,PS))

    Flow(pp+'s0',pp,'spectra all=y')

    scanplot = '''
    window min1=%g max1=%g |
    byte gainpanel=all allpos=y |
    grey3 frame1=%d frame3=%d frame2=%d color=j flat=n
    label1=Time unit1=s label3="In-line" label2="Relative Gamma"
    wanttitle=n
    ''' % (tmin,tmax,frame1,(trace-o2)/j2,ng/2)

    simplot = '''
    window min1=%g max1=%g |
    grey title="%s" allpos=y 
    color=j clip=1
    label1="Time (s)" 
    ''' % (tmin,tmax,'%s')

    warpit = warping(niter,200,200)

    for i in range(iter):
        wrp = warp 
        
        #################
        # INITIAL WARPING
        #################

        def n(s):
            return '%s-%s-%d' % (name,s,i)

        psw = n('psw')
        Flow(psw,[ps,pp,wrp],warp0)
        Plot(psw,plot('Warped ' + PS))
        
        dif = n('dif')
        Plot(dif,[psw,pp],
             'add scale=1,-1 ${SOURCES[1]} | ' + plot('Difference'))

        gamma = n('gamma')
        Flow(gamma,wrp,warp2gamma(ss))
        Plot(gamma,vplot)

        Result(psw,[pp,psw,dif,gamma],'TwoRows')

        psw1 = n('psw1')
        pp1 = n('pp1')
        Flow(pp1,pp,'window n2=1 f2=286')
        Flow(psw1,psw,'window n2=1 f2=286')

        ppps = n('ppps')
        Result(ppps,[psw1,pp1],
               '''
               add scale=1,-1 ${SOURCES[1]} |
               cat ${SOURCES[0]} ${SOURCES[1]} axis=2 |
               dots gaineach=n Xscreenwd=9.225 Xscreenht=5.2 Xyyscale=0.8
               labels="Difference:SS warped:PP" label1=Time unit1=s
               title="LTF transform balancing"
               ''')

        ####################
        # SPECTRAL BALANCING
        ####################

        si = n('si')
        Flow(si,psw,ifreq)
        Plot(si,freqplot(PS + ' Local Frequency'))

##         msk = n('msk')
##         Flow(msk,[si,pp+'i'],getmask)

##         sr = n('sr')
##         pr = n('pr')

##         Flow(sr+'0',[msk,si,pp+'i'],psrect(frect))
##         Flow(sr,[psw,sr+'0',pp],balance)

##         Flow(pr+'0',[msk,si,pp+'i'],pprect(frect))
##         Flow(pr,[pp,pr+'0',pp],balance)

        ppltft = n('ppltft')
        Flow(ppltft,pp,
             '''
             ltft rect=%d verb=n | transp
             ''' % frect,split=[2,471],reduce="cat axis=3")
        Result(ppltft,
               '''
               transp | window max2=100 |
               math output="abs(input)" | real |
               byte allpos=y gainpanel=40 pclip=100 |
               grey3 transp=n yreverse=n color=j 
               title="PP before matched" label1=Time unit1=s
               frame1=512 frame2=100 frame3=235 point1=0.8 point2=0.3
               labelfat=3 font=2 titlefat=3
               ''')

        ppltftspe = n('ppltftspe')
        Flow(ppltftspe,ppltft,
             '''
             math output="abs(input)" | real 
             ''')
        pswltft = n('pswltft')
        Flow(pswltft,psw,
             '''
             ltft rect=%d verb=n | transp
             ''' % frect,split=[2,471],reduce="cat axis=3")
        Result(pswltft,
               '''
               transp | window max2=100 |
               math output="abs(input)" | real |
               byte allpos=y gainpanel=40 pclip=100 |
               grey3 transp=n yreverse=n color=j 
               title="PS warped before matched" label1=Time unit1=s
               frame1=512 frame2=100 frame3=235 point1=0.8 point2=0.3
               labelfat=3 font=2 titlefat=3
               ''')
        pswltftspe = n('pswltftspe')
        Flow(pswltftspe,pswltft,
             '''
             math output="abs(input)" | real 
             ''')
        
        pprick = n('pprick')
        Flow(pprick,ppltftspe,
             '''
             ricker niter=1000 ma=$TARGET verb=n m=40
             ''',stdout=0,split=[2,1024])
        pswrick = n('pswrick')
        Flow(pswrick,pswltftspe,
             '''
             ricker niter=1000 ma=$TARGET verb=n m=40
             ''',stdout=0,split=[2,1024])
        
        ppshape = n('ppshape')
        pswshape = n('pswshape')
        Flow([ppshape,pswshape],[ppltft,pswltft,pprick,pswrick],
             '''
             freshape in2=${SOURCES[1]} ma=${SOURCES[2]}
             ma2=${SOURCES[3]} out2=${TARGETS[1]}
             ''')
        
##         Result('ppmorphltft',
##                '''       
##                transp | window max2=100 |
##                math output="abs(input)" | real |
##                grey transp=n yreverse=n color=j
##                title="PP matched" label1=Time unit1=s
##                screenratio=0.35 screenht=5
##                labelfat=3 font=2 titlefat=3 parallel2=n scalebar=y
##                ''')     
##         Result('psmorphltft',
##                '''       
##                transp | window max2=100 |
##                math output="abs(input)" | real |
##                grey transp=n yreverse=n color=j
##                title="PS matched" label1=Time unit1=s
##                screenratio=0.35 screenht=5
##                labelfat=3 font=2 titlefat=3 parallel2=n scalebar=y
##                ''')     

        sr = n('sr')
        pr = n('pr')
        Flow(pr,ppshape,'transp | ltft inv=y verb=y')
        Flow(sr,pswshape,'transp | ltft inv=y verb=y')
        
##         Plot('ippmorphltft',
##              '''
##              graph plotfat=10 label1=Time unit1=s unit2=
##              wanttitle=n labelsz=10 label2="PP matched"
##              ''')
##         Plot('ipsmorphltft',
##              '''
##              graph plotfat=10 label1=Time unit1=s unit2=
##              wanttitle=n labelsz=10 label2="PS matched"
##              ''')
##         Result('match','ippmorphltft ipsmorphltft','OverUnderAniso')
        sr1 = n('sr1')
        pr1 = n('pr1')
        Flow(pr1,pr,'window n2=1 f2=286')
        Flow(sr1,sr,'window n2=1 f2=286')

        ltftppps = n('ltftppps')
        Result(ltftppps,[sr1,pr1],
               '''
               add scale=1,-1 ${SOURCES[1]} |
               cat ${SOURCES[0]} ${SOURCES[1]} axis=2 |
               dots gaineach=0 Xscreenwd=9.225 Xscreenht=5.2 Xyyscale=0.8
               labels="Difference:PS warped:PP" label1=Time unit1=s
               title="LTF transform balancing"
               ''')
        
##         Plot('ltftbefore','pp morph',
##              '''
##              cat axis=2 ${SOURCES[1]} |
##              spectra | window max1=100 |
##              graph dash=0,1 title="Before LTFT matching" max2=2 min2=0
##              label1=Frequency unit1=Hz label2= unit2=
##              ''')
##         Plot('ltftafter','ippmorphltft ipsmorphltft',
##              '''
##              cat axis=2 ${SOURCES[1]} |
##              spectra | window max1=100 |
##              graph dash=0,1 title="After LTFT matching" max2=2 min2=0
##              label1=Frequency unit1=Hz label2= unit2=
##              ''')
##         Result('ltftspectra','ltftbefore ltftafter','SideBySideAniso')
        pi = n('pi')
        Flow(pi,pr,ifreq)
        Flow(si+'2',sr,ifreq)
        
        Plot(si+'2',freqplot(PS + ' Local Frequency'))
        Plot(pi,freqplot('PP Local Frequency'))
        Result(si,[pp+'i',si,pi,si+'2'],'TwoRows')

        s0 = psw+'s0'
        Flow(s0,psw,'spectra all=y')
        Plot(s0,[pp+'s0',s0],specplot('Before'))
        
        s1 = psw+'s1'
        Flow(s1,sr,'spectra all=y')
        Flow(pr+'s1',pr,'spectra all=y')
        Plot(s1,[pr+'s1',s1],specplot('After'))

        Result(n('sp'),[s0,s1],'SideBySideIso')

        if i == 0:
            in0 = n('in0')
            Flow(pr+'in0',pr,interg)
            Flow(sr+'in0',sr,interg)
            Plot(in0,    [pr+'in0',sr+'in0'],giplot('Before'))

            Flow(pr+'in0w',pr,interw)
            Flow(sr+'in0w',sr,interw)
            Plot(in0+'w',[pr+'in0w',sr+'in0w'],wiplot('Before'))
            
            Result(in0,in0,'Overlay')
            Result(in0+'w',in0+'w','Overlay')

            sim0 = n('sim0')
            Flow(sim0,[pr,sr],simil(rect1,rect2))
            Result(sim0,simplot % 'Before')
           
        Plot(sr,plot('Warped and Balanced ' + PS))
        Plot(pr,plot('Balanced PP'))
        
        dif = dif+'2'
        Plot(dif,[sr,pr],
             'add scale=1,-1 ${SOURCES[1]} | ' + plot('Difference'))

        Result(sr,[pr,sr,dif,gamma],'TwoRows')        

        ############
        # GAMMA SCAN
        ############

        g1 = 2-g0
        warpscan2 = warpscan(ng,g0,g1,rect1,1,int(0.5+rect2/j2))
        
        sc = n('sc')

        Flow(sr+'2',sr,'window j2=%d' % j2)
        Flow(pr+'2',pr,'window j2=%d' % j2)
        
        Flow(sc,[sr+'2',pr+'2'],warpscan2)
        Result(sc,scanplot)
        
        pk = n('pk')

        if i==0:
            Flow(pk+'0',sc,pick(max(pmin,g0),min(pmax,g1),
                                rect1,4*rect2/j2,an=an))
        else:
            Flow(pk+'0',sc,pick(g0,g1,rect1,4*rect2/j2,an=an))

        Flow(pk,pk+'0',
             '''
             transp memsize=500 |
             spline n1=%d d1=1 o1=%g |
             transp memsize=500  |
             math output="(input-1)*x1"
             ''' % (nx,o2))

        #########
        # WARPING
        #########

        warp = n('wrp')
        Flow([warp,psw+'2'],[sr,pr,pk,wrp],warpit,stdout=-1)
        Plot(psw+'2',plot('Warped ' + PS))
        
        dif = n('dif2')
        Plot(dif,[psw+'2',pr],
             'add scale=1,-1 ${SOURCES[1]} | ' + plot('Difference'))

        gamma = n('gamma2')
        Flow(gamma,warp,warp2gamma(ss))
        Plot(gamma,vplot)

        Result(psw+'2',[pr,psw+'2',dif,gamma],'TwoRows')

        if i == iter-1:
            in1 = n('in1')
            Flow(pr+'in1',pr,interg)
            Flow(psw+'2in1',psw+'2',interg)
            Plot(in1,[pr+'in1',psw+'2in1'],giplot('After'))
            Flow(pr+'in1w',pr,interw)
            Flow(psw+'2in1w',psw+'2',interw)
            Plot(in1+'w',[pr+'in1w',psw+'2in1w'],wiplot('After'))
            Result(in1,in1,'Overlay')
            Result(in1+'w',in1+'w','Overlay')
            Result(in0+'1',[in0,in1],'SideBySideIso')
            Result(in0+'1w',[in0+'w',in1+'w'],'OverUnderAniso')

            sim1 = n('sim1')
            Flow(sim1,[pr,psw+'2'],simil(rect1,rect2))
            Result(sim1,simplot % 'After')

            Flow(psw+'1',[ps,pp,warp],warp0)
            Result(psw+'1',plot('Warped ' + PS))

            rt = n('rt')
            Flow(psw+'i',psw+'1',ifreq)            
            Flow(rt,psw+'i',
                 '''
                 math output="sqrt(1+12*(1/input^2-1/%g^2))" |
                 dd type=int
                 ''' % (fmax*2*math.pi*dt))

            dl = n('dl')
            Flow(dl,[psw+'1',rt],
                 '''
                 deblur rect=${SOURCES[1]}
                 verb=y niter=100 eps=0.04 nliter=1
                 ''')
            Result(dl,
                   '''
                   window min1=%g max1=%g |
                   grey title="Deblurred %s" label1="Time (s)"
                   ''' % (tmin,tmax,PS))

            Flow('e'+gamma,warp,warp2egamma(ss))
            Result(gamma,'e'+gamma,vplot)
        
        g0 = (g0+1)*0.5

def warp1(name,      # name prefix
          pp,ps,     # PP and PS images
          warp,      # initial warp
          tmax,      # maximum time for display
          tmin=0,    # minimum time for display
          gmin=1,    # minimum gamma
          gmax=4,    # maximum gamma
          niter=20,  # warping iterations
          dt=0.004,  # time sampling
          fmin=0,    # minimum frequency
          fmax=40,   # maximum frequency
          frect=12,  # frequency smoothing
          ng=101,    # number of gammas
          g0=0.96,   # first gamma
          pmin=0,    # minimum gamma for picking
          pmax=2,    # maximum gamma for picking
          an=0.5,    # anisotropy for picking
          rect1=50,  # vertical smoothing
          iter=2,    # number of iterations
          ss=0
          ):

    if version.old_version():
        return # think how to do it better

    graph = '''
    graph wanttitle=n min2=%g max2=%g min1=%g max1=%g
    wherexlabel=t wheretitle=b crowd=0.8 label2="Vp/Vs"
    ''' % (gmin,gmax,tmin,tmax)

    dplot ='''
    add scale=1,-1 ${SOURCES[1]} |
    cat ${SOURCES[0]} ${SOURCES[1]} axis=2 |
    window min1=%g max1=%g |
    dots gaineach=0
    labels="Difference:PS warped:PP" label1=Time unit1=s
    ''' % (tmin,tmax)

    def iphase(title):
        return '''
        cat axis=2 ${SOURCES[1]} |
        scale dscale=%g | 
        graph title="Local Frequency (%s)" label1="Time (s)"
        min2=%g max2=%g min1=%g max1=%g
        dash=0,1 label2="Frequency (Hz)"
        ''' % (0.5/(math.pi*dt),title,fmin,fmax,tmin,tmax)

    warpit = warping(niter,200)

    for i in range(iter):
        #################
        # INITIAL WARPING
        #################
        wrp = warp 
      
        def showpick(case):
            return '''
            graph transp=y min2=%g max2=%g min1=%g max1=%g
            yreverse=y plotcol=%d plotfat=%d 
            wantaxis=n wanttitle=n pad=n
            ''' % (g0,g1,tmin,tmax,(7,0)[case],(5,1)[case])

        def n(s):
            return '%s-%s-%d' % (name,s,i)

        gamma = n('gamma')
        Flow(gamma,wrp,warp2gamma(ss));
        Plot(gamma,graph)

        psw = n('psw')
        Flow(psw,[ps,pp,wrp],warp0)
        Plot(psw,[psw,pp],dplot)

        Result(psw,[gamma,psw],'OverUnderAniso')

        ####################
        # SPECTRAL BALANCING
        ####################

        ppft = n('ppft')
        psft = n('psft')

        ltft = 'ltft rect=%d | transp' % frect

        Flow(ppft,pp,ltft)
        Flow(psft,psw,ltft)

        Flow(ppft+'a',ppft,'math output="abs(input)" | real')
        Flow(psft+'a',psft,'math output="abs(input)" | real')

        ftplot = '''
        window min1=%g max1=%g max2=%g |
        grey allpos=y color=j labelfat=3 font=2 titlefat=3
        screenht=6. screenratio=0.45
        ''' % (fmin,fmax,tmax)

        Plot(ppft+'a',ftplot
             +'title="PP before balancing" scalebar=y maxval=0.014')
        Result(ppft+'a','Overlay')
        Plot(psft+'a',ftplot
             +'title="Warped SS before balancing" scalebar=y maxval=0.014')
        Result(psft+'a','Overlay')
        Result(n('ft0'),[ppft+'a',psft+'a'],'OverUnderAniso')

        pprick = n('pprick')
        psrick = n('psrick')

        Flow(pprick,ppft+'a',
             '''
             ricker niter=1000 ma=$TARGET verb=n m=40
             ''',stdout=0, split=[2,1024])
        Flow(psrick,psft+'a',
             '''
             ricker niter=1000 ma=$TARGET verb=n m=40
             ''',stdout=0, split=[2,1024])

        rickplot = '''
        cat axis=3 ${SOURCES[1]} | window n1=1 max2=%g | 
        math output="sqrt(input)" |
        graph title="Dominant Frequency" 
        label2=Frequency unit2=Hz min2=%g max2=%g
        ''' % (tmax,fmin,fmax)

        Result(n('rick'),[pprick,psrick],rickplot)
        
        Flow([ppft+'b',psft+'b'],[ppft,psft,pprick,psrick],
             '''
             freshape in2=${SOURCES[1]} ma=${SOURCES[2]}
             ma2=${SOURCES[3]} out2=${TARGETS[1]}
             ''')
        Flow(ppft+'c',ppft+'b','math output="abs(input)" | real')
        Flow(psft+'c',psft+'b','math output="abs(input)" | real')

        Plot(ppft+'c',ftplot
             +'title="PP after balancing" scalebar=y maxval=0.014')
        Result(ppft+'c','Overlay')
        Plot(psft+'c',ftplot
             +'title="Warped SS after balancing" scalebar=y maxval=0.014')
        Result(psft+'c','Overlay')
        Result(n('fta'),[ppft+'c',psft+'c'],'OverUnderAniso')

        sr = n('sr')
        pr = n('pr')

        Flow(pr,ppft+'b','transp | ltft inv=y')
        Flow(sr,psft+'b','transp | ltft inv=y')

        Plot(psw+'1',[sr,pr],dplot)
        Result(psw+'1',[gamma,psw+'1'],'OverUnderAniso')

        ############
        # GAMMA SCAN
        ############

        g1 = 2-g0
        
        warpscan1 = warpscan(2*ng,g0,g1,rect1)
        
        greyscan = '''
        window min1=%g max1=%g |
        grey title="Gamma scan" allpos=y 
        min2=%g max2=%g
        color=j pclip=100
        label1="Time (s)" label2="Gamma"
        ''' % (tmin,tmax,g0,g1)

        scn = n('scn')
        Flow(scn,[sr,pr],warpscan1)
        Plot(scn,greyscan)

        pik = n('pik')

        if i==0:
            Flow(pik+'0',scn,pick(max(pmin,g0),min(pmax,g1),2*rect1,an=an))
        else:
            Flow(pik+'0',scn,pick(g0,g1,2*rect1,an=an))

        Flow(pik,pik+'0','math output="(input-1)*x1" ')
        Plot(pik,pik+'0',showpick(0))
        Plot(pik+'0',showpick(1))
        Result(scn,[scn,pik,pik+'0'],'Overlay')

        #########
        # WARPING
        #########

        warp = n('wrp')

        Flow([warp,psw+'2'],[sr,pr,pik,wrp],warpit,stdout=-1)
        Flow(gamma+'2',warp,warp2gamma(ss))
        Plot(gamma+'2',graph)
        Plot(psw+'2',[psw+'2',pr],dplot)
        Result(psw+'2',[gamma+'2',psw+'2'],'OverUnderAniso')
        
        g0 = (g0+1)*0.5

        ############
        # Comparison
        ############

        Flow(psw+'line',psw,'math output=0.1')

        Flow(psw+'0c',[psw,pp],
             '''
             add scale=1,-1 ${SOURCES[1]} |
             cat ${SOURCES[0]} ${SOURCES[1]} axis=2
             ''')
        Flow(psw+'1c',[sr,pr],
             '''
             add scale=1,-1 ${SOURCES[1]} |
             cat ${SOURCES[0]} ${SOURCES[1]} axis=2
             ''')
        Flow(psw+'2c',[psw+'2',pr],
             '''
             add scale=1,-1 ${SOURCES[1]} |
             cat ${SOURCES[0]} ${SOURCES[1]} axis=2
             ''')
        
        Result(psw+'c',[psw+'2c',psw+'line',psw+'1c',psw+'line',psw+'0c'],
               '''
               cat ${SOURCES[1:5]} axis=2 |
               window min1=%g max1=%g |
               dots gaineach=0 labelfat=4 font=2 titlefat=4
               labels="Difference 3:SS registered:PP balanced:
               :Difference 2:SS balanced:PP balanced:
               :Difference 1:SS initial warped:PP"
               label1=Time unit1=s
               ''' % (tmin,tmax))
               
