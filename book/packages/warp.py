from rsfproj import *
import math

warp0 = '''
warp1 other=${SOURCES[1]} warpin=${SOURCES[2]}
verb=1 nliter=0
'''

dplot ='''
add scale=1,-1 ${SOURCES[1]} |
cat ${SOURCES[0]} ${SOURCES[1]} axis=2 |
dots gaineach=0
labels="Difference:PS warped:PP" label1="Time (s)"
''' 

getmask = 'add scale=1,-1 ${SOURCES[1]} | mask min=0 | dd type=float'

psrect = '''
math min=${SOURCES[2]} max=${SOURCES[1]}
output="sqrt(1+12*(1/min^2-1/max^2)*input)" | dd type=int
'''

pprect = '''
math min=${SOURCES[2]} max=${SOURCES[1]}
output="sqrt(1+12*(1/max^2-1/min^2)*(1-input))" | dd type=int
'''

balance = '''
nsmooth1 rect=${SOURCES[1]} |
abalance rect1=100 order=100 other=${SOURCES[2]}
'''

def warping(niter,rect1=1,rect2=1,rect3=1):
    return '''
    warp1 other=${SOURCES[1]} warpin=${SOURCES[2]}
    warpout=www.rsf
    verb=1 nliter=%d noamp=1 rect1=%d rect2=%d rect3=%d > ${TARGETS[1]} &&
    warpadd < ${SOURCES[3]} add=www.rsf > $TARGET &&
    rm www.rsf
    ''' % (niter,rect1,rect2,rect3)

def pick(rect1=1,rect2=1,rect3=1):
    return 'pick rect1=%d rect2=%d rect3=%d an=0.5 | window' % (rect1,rect2,rect3)

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

def warp2(name,       # name prefix
          pp,ps,      # PP and PS images
          warp,       # initial warp
          nx,         # number of traces
          tmax,       # maximum time for display
          tmin=0,     # minimum time for display
          trace=None, # seleted trace
          gmin=1,     # minimum gamma
          gmax=4,     # maximum gamma
          dt=0.004,   # time sampling
          fmin=0,     # minimum frequency
          fmax=40,    # maximum frequency
          frame1=10,  # time frame          
          ng=101,     # number of gammas
          g0=0.96,    # first gamma
          rect1=50,   # vertical smoothing
          rect2=50,   # lateral smoothing
          iter=2,     # number of iterations
          ss=0,       # PP-PS (0) or PP-SS (1)
          inter=1,    # interleaving
          clip=6      # display clip
          ):

    interg = 'pad n2=%d | put n2=%d n3=%d | stack' % ((nx/inter+1)*inter,inter,nx/inter+1)
    inter = 2*inter    
    interw = 'pad n2=%d | put n2=%d n3=%d | stack' % ((nx/inter+1)*inter,inter,nx/inter+1)
      
    if trace:
        for case in (pp,ps,warp):
            Flow(case+'1',case,'window n2=1 min2=%d' % trace)
        warp1(name+'t',pp+'1',ps+'1',warp+'1',tmax,tmin,
              gmin,gmax,dt,fmin,fmax,ng,g0,rect1,iter,ss)
    else:
        trace=10

    def plot(title):
        return '''
        window min1=%g max1=%g |
        grey title="%s" label1="Time (s)" clip=%g 
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
        dash=0,1
        ''' % (title,4*fmax)

    def giplot(title):
        return '''
        interleave axis=2 ${SOURCES[1]} |
        window min1=%g max1=%g |
        grey
        title="Interleaved (%s)"
        label1="Time (s)" label2="In-line"
        ''' % (tmin,tmax,title)

    def wiplot(title):
        return '''
        interleave axis=2 ${SOURCES[1]} |
        window min1=%g max1=%g |
        wiggle poly=y transp=y yreverse=y
        title="Interleaved (%s)"
        label1="Time (s)" label2="In-line"
        ''' % (tmin,tmax,title)
        
    Plot(pp,plot('PP'))
    Flow(pp+'i',pp,ifreq)
    Plot(pp+'i',freqplot('PP Local Frequency'))

    PS = ('PS','SS')[ss]

    Plot(ps,'grey title=%s label1="Time (s)" ' % PS)

    Flow(pp+'s0',pp,'spectra all=y')

    scanplot = '''
    window min1=%g max1=%g |
    byte gainpanel=all allpos=y |
    grey3 frame1=%d frame3=%d frame2=%d color=j flat=n
    label1="Time (s)" label3="In-line" label2="Relative Gamma"
    wanttitle=n
    ''' % (tmin,tmax,frame1,trace,ng/2)

    warpit = warping(2,200,200)

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

        ####################
        # SPECTRAL BALANCING
        ####################

        si = n('si')
        Flow(si,psw,ifreq)
        Plot(si,freqplot(PS + ' Local Frequency'))
        Result(si,[pp+'i',si],'OverUnderAniso')

        msk = n('msk')
        Flow(msk,[si,pp+'i'],getmask)

        sr = n('sr')
        pr = n('pr')

        Flow(sr+'0',[msk,si,pp+'i'],psrect)
        Flow(sr,[psw,sr+'0',pp],balance)

        Flow(pr+'0',[msk,si,pp+'i'],pprect)
        Flow(pr,[pp,pr+'0',pp],balance)

        si = si+'2'
        pi = n('pi')
        Flow(pi,pr,ifreq)
        Flow(si,sr,ifreq)
        
        Plot(si,freqplot(PS + ' Local Frequency'))
        Plot(pi,freqplot('PP Local Frequency'))
        Result(si,[pi,si],'OverUnderAniso')

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
        warpscan2 = warpscan(ng,g0,g1,rect1,1,rect2)
        
        sc = n('sc')
        Flow(sc,[sr,pr],warpscan2)
        Result(sc,scanplot)
        
        pk = n('pk')
        Flow(pk+'0',sc,pick(rect1,4*rect2))
        Flow(pk,pk+'0','math output="(input-1)*x1" ')

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
            Result(in0+'1',[in0,in1],'SideBySideAniso')
            Result(in0+'1w',[in0+'w',in1+'w'],'OverUnderAniso')

            Flow(psw+'1',[ps,pp,warp],warp0)
            Result(psw+'1',plot('Warped ' + PS))

            rt = n('rt')
            Flow(psw+'i',psw+'1',ifreq)            
            Flow(rt,psw+'i','math output="sqrt(1+12*(1/input^2-1/%g^2))" | dd type=int' % (fmax*2*math.pi*dt))

            dl = n('dl')
            Flow(dl,[psw+'1',rt],'deblur rect=${SOURCES[1]} verb=y niter=100 eps=0.04 nliter=1')
            Result(dl,'''
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
          dt=0.004,  # time sampling
          fmin=0,    # minimum frequency
          fmax=40,   # maximum frequency
          ng=101,    # number of gammas
          g0=0.96,   # first gamma
          rect1=50,  # vertical smoothing
          iter=2,    # number of iterations
          ss=0
          ):

    graph = '''
    graph wanttitle=n min2=%g max2=%g min1=%g max1=%g
    wherexlabel=t wheretitle=b crowd=0.8 label2="Vp/Vs"
    ''' % (gmin,gmax,tmin,tmax)

    iphase = '''
    cat axis=2 ${SOURCES[1]} |
    scale dscale=%g | 
    graph title="Local Frequency" label1="Time (s)"
    min2=%g max2=%g min1=%g max1=%g
    dash=0,1 label2="Frequency (Hz)"
    ''' % (0.5/(math.pi*dt),fmin,fmax,tmin,tmax)

    warpit = warping(2,200)

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

        phase = 'iphase rect1=%d' % (2*rect1)

        ifr = n('ifr')
        ppi = n('ppi')
        psi = n('psi')
        Flow(ppi,pp, phase)
        Flow(psi,psw,phase)
        Result(ifr,[ppi,psi],iphase)

        msk = n('msk')
        Flow(msk,[psi,ppi],getmask)

        sr = n('sr')
        pr = n('pr')

        Flow(sr+'0',[msk,psi,ppi],psrect)
        Flow(sr,[psw,sr+'0',pp],balance)
    
        Flow(pr+'0',[msk,psi,ppi],pprect)
        Flow(pr,[pp,pr+'0',pp],balance)

        Plot(psw+'1',[sr,pr],dplot)
        Result(psw+'1',[gamma,psw+'1'],'OverUnderAniso')

        for case in (sr,pr):
            Flow(case+'i',case,phase)
        Result(ifr+'1',[pr+'i',sr+'i'],iphase)

        ############
        # GAMMA SCAN
        ############

        g1 = 2-g0
        
        warpscan1 = warpscan(ng,g0,g1,rect1)
        
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
        Flow(pik+'0',scn,pick(2*rect1))
        Flow(pik,pik+'0','math output="(input-1)*x1" ')
        Plot(pik,pik+'0',showpick(0))
        Plot(pik+'0',showpick(1))
        Result(scn,[scn,pik,pik+'0'],'Overlay')

        #########
        # WARPING
        #########

        warp = n('wrp')

        Flow([warp,psw+'2'],[sr,pr,pik,wrp],warpit,stdout=-1)
        Flow(gamma+'2',wrp,warp2gamma(ss))
        Plot(gamma+'2',graph)
        Plot(psw+'2',[psw+'2',pr],dplot)
        Result(psw+'2',[gamma+'2',psw+'2'],'OverUnderAniso')
        
        g0 = (g0+1)*0.5
