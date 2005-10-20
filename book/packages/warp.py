from rsfproj import *
import math

warp2gamma = '''
math output="input+x1" |
smoothder |
math output="2*input-1" 
'''

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

inter = '''
pad end2=1 > shift.rsf &&
pad < $SOURCE beg2=1 | add shift.rsf | window f2=1 j2=2 > $TARGET &&
rm shift.rsf
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

def warp2(name,       # name prefix
          pp,ps,      # PP and PS images
          warp,       # initial warp
          trace=None, # seleted trace
          gmin=1,     # minimum gamma
          gmax=4,     # maximum gamma
          dt=0.004,   # time sampling
          fmin=0,     # minimum frequency
          fmax=40,    # maximum frequency
          ng=101,     # number of gammas
          g0=0.96,    # first gamma
          rect1=50,   # vertical smoothing
          rect2=50,   # lateral smoothing
          iter=2,     # number of iterations
          clip=6
          ):

    if trace:
        for case in (pp,ps,warp):
            Flow(case+'1',case,'window n2=1 min2=%d' % trace)
#        warp1(name+'t',pp+'1',ps+'1',warp+'1',gmin,gmax,dt,fmin,fmax,ng,g0,rect1,iter)

    def plot(title):
        return '''
        grey title="%s" label1="Time (s)" clip=%g
        ''' % (title,clip)

    vplot = plot('Vp/Vs') + '''
    clip=%g scalebar=y color=j bias=%g minval=%g maxval=%g
    ''' % (0.5*(gmax-gmin),0.5*(gmin+gmax),gmin,gmax)

    balance = '''
    nsmooth1 rect=${SOURCES[1]} |
    abalance rect1=%d rect2=%d order=100 other=${SOURCES[2]}
    ''' % (rect1,rect2)

    ifreq = 'iphase rect1=100 rect2=10 order=100'

    def freqplot(title):
        return '''
        scale scale dscale=%g |
        %s clip=%g bias=%g color=j scalebar=y barlabel="Frequency (Hz)"
        ''' % (0.5/(math.pi*dt),plot(title),(fmax-fmin)*0.5,(fmax+fmin)*0.5)

    def specplot(title):
        return '''
        cat axis=2 ${SOURCES[1]} |
        graph title="%s" max1=%g label1="Frequency (Hz)"
        dash=0,1
        ''' % (title,4*fmax)

    Plot(pp,plot('PP'))
    Flow(pp+'i',pp,ifreq)
    Plot(pp+'i',freqplot('PP Local Frequency'))

    Flow(pp+'s0',pp,'spectra all=y')

    g1 = 2-g0
    warpscan2 = warpscan(ng,g0,g1,rect1,1,rect2)

    scanplot = '''
    byte gainpanel=all allpos=y |
    grey3 frame1=750 frame3=10 frame2=25 color=j flat=n
    label1="Time (s)" label3="In-line" label2="Relative Gamma"
    wanttitle=n
    '''

    wrp = warp 
    for i in range(iter):
        #################
        # INITIAL WARPING
        #################

        def n(s):
            return '%s-%s-%d' % (name,s,i)

        psw = n('psw')
        Flow(psw,[ps,pp,wrp],warp0)
        Plot(psw,plot('Warped PS'))
        
        dif = n('dif')
        Plot(dif,[psw,pp],
             'add scale=1,-1 ${SOURCES[1]} | ' + plot('Difference'))

        gamma = n('gamma')
        Flow(gamma,wrp,warp2gamma)
        Plot(gamma,vplot)

        Result(psw,[pp,psw,dif,gamma],'TwoRows')

        ####################
        # SPECTRAL BALANCING
        ####################

        si = n('si')
        Flow(si,psw,ifreq)
        Plot(si,freqplot('PS Local Frequency'))
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
        
        Plot(si,freqplot('PS Local Frequency'))
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

        in0 = n('in0')
        Flow(pr+'in',pr,inter,stdout=-1)
        Flow(sr+'in',sr,inter,stdout=-1)
        Plot(in0,[pr+'in',sr+'in'],
             '''
             interleave axis=2 ${SOURCES[1]} |
             grey
             title='Interleaved (Before)'
             label1='Time (s)' label2='In-line'
             ''')
        Plot(in0+'w',[pr+'in',sr+'in'],
             '''
             interleave axis=2 ${SOURCES[1]} |
             wiggle poly=y transp=y yreverse=y
             title='Interleaved (Before)'
             label1='Time (s)' label2='In-line'
             ''')

        Plot(sr,plot('Warped and Balanced PS'))
        Plot(pr,plot('Balanced PP'))
        
        dif = dif+'2'
        Plot(dif,[sr,pr],
             'add scale=1,-1 ${SOURCES[1]} | ' + plot('Difference'))

        Result(sr,[pr,sr,dif,gamma],'TwoRows')        

        ############
        # GAMMA SCAN
        ############
        sc = n('sc')
        Flow(sc,[sr,pr],warpscan2)
        Result(sc,scanplot)
        


def warp1(name,      # name prefix
          pp,ps,     # PP and PS images
          warp,      # initial warp
          gmin=1,    # minimum gamma
          gmax=4,    # maximum gamma
          dt=0.004,  # time sampling
          fmin=0,    # minimum frequency
          fmax=40,   # maximum frequency
          ng=101,    # number of gammas
          g0=0.96,   # first gamma
          rect1=50,  # vertical smoothing
          iter=2     # number of iterations
          ):

    graph = '''
    graph wanttitle=n min2=%g max2=%g
    wherexlabel=t wheretitle=b crowd=0.8 label2="Vp/Vs"
    ''' % (gmin,gmax)

    iphase = '''
    cat axis=2 ${SOURCES[1]} |
    scale dscale=%g | 
    graph title="Local Frequency" label1="Time (s)"
    min2=%g max2=%g dash=0,1 label2="Frequency (Hz)"
    ''' % (0.5/(math.pi*dt),fmin,fmax)

    warp1 = warping(1,200)

    wrp = warp 
    for i in range(iter):
        #################
        # INITIAL WARPING
        #################

        g1 = 2-g0

        warpscan1 = warpscan(ng,g0,g1,rect1)

        greyscan = '''
        grey title="Gamma scan" allpos=y 
        min2=%g max2=%g
        color=j pclip=100
        label1="Time (s)" label2="Gamma"
        ''' % (g0,g1)
        
        def showpick(case):
            return '''
            graph transp=y min2=%g max2=%g
            yreverse=y plotcol=%d plotfat=%d 
            wantaxis=n wanttitle=n pad=n
            ''' % (g0,g1,(7,0)[case],(5,1)[case])

        def n(s):
            return '%s-%s-%d' % (name,s,i)

        gamma = n('gamma')
        Flow(gamma,wrp,warp2gamma);
        Plot(gamma,graph)

        psw = n('psw')
        Flow(psw,[ps,pp,wrp],warp0)
        Plot(psw,[psw,pp],dplot)

        Result(psw,[gamma,psw],'OverUnderAniso')

        ####################
        # SPECTRAL BALANCING
        ####################

        ifr = n('ifr')
        ppi = n('ppi')
        psi = n('psi')
        Flow(ppi,pp, 'iphase rect1=100')
        Flow(psi,psw,'iphase rect1=100')
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
            Flow(case+'i',case,'iphase rect1=100')
        Result(ifr+'1',[pr+'i',sr+'i'],iphase)

        ############
        # GAMMA SCAN
        ############

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

        wrp = n('wrp')

        Flow([wrp,psw+'2'],[sr,pr,pik,warp],warp1,stdout=-1)
        Flow(gamma+'2',wrp,warp2gamma)
        Plot(gamma+'2',graph)
        Plot(psw+'2',[psw+'2',pr],dplot)
        Result(psw+'2',[gamma+'2',psw+'2'],'OverUnderAniso')
        
        g0 = (g0+1)*0.5
