from rsf.proj import *
import math

pi2 = 2*math.pi

def NAR(data,
        dt,           # time sampling
        tmax,         # maximum time
        fmax,         # maximum frequency
        rect1,        # smoothing radius for frequency estimation
        rect2,        # smoothing radius for amplitude estimation
        ns=8,         # number of components
        list=[5,6,7], # significant components
        amin=0.2,     # minimum aplitude for significant components
        ):
    'Non-stationary auto-regressions'

    def n(name):
        return '%s-%s' % (data,name)

    Flow(n('shift'),data,'shift1 ns=%d' % ns)

    # analytical trace
    Flow(n('itrace'),data,'envelope hilb=y')
    Flow(n('ctrace'),[data,n('itrace')],'cmplx ${SOURCES[1]}')

    Flow(n('ishift'),n('shift'),'envelope hilb=y')
    Flow(n('cshift'),[n('shift'),n('ishift')],'cmplx ${SOURCES[1]}')

    Flow([n('cpef'),n('cpre')],[n('cshift'),n('ctrace')],
         'clpf match=${SOURCES[1]} rect1=%d pred=${TARGETS[1]}' % rect1)
    Flow(n('cerr'),[n('cpre'),n('ctrace')],'add scale=-1,1 ${SOURCES[1]}')

    Result(n('cerr'),'real | graph title="Nonstationary Deconvolution" ')

    # Find complex polynomial roots

    Flow(n('cpoly'),n('cpef'),
         'window n2=1 | math output=-1 | cat axis=2 $SOURCE')
    Flow(n('croots'),n('cpoly'),
         'transp | roots verb=n niter=100 sort=r | transp')

    Flow(n('group'),n('croots'),
         'math output="-arg(input)/%g" | real ' % (pi2*dt))
    Result(n('group'),
       '''
       graph title=Frequencies yreverse=y pad=n
       wanttitle=n scalebar=n bartype=v
       plotfat=3 grid=y label2=Frequency unit2=Hz 
       min2=0 max2=%g plotcol=6,5,4,6,5,4,6,5,4
       screenratio=0.8 crowd1=0.75 crowd2=0.3 labelsz=6
        ''' % fmax)

    Flow(n('freqs'),n('group'),
         'causint | math output="input*%g/(x1+%g)" ' % (dt,dt))

    Result(n('freqs'),
       '''
       graph title=Frequencies yreverse=y pad=n wanttitle=n
       scalebar=n bartype=v
       plotfat=3 grid=y label2=Frequency unit2=Hz 
       min2=0 max2=%g
       ''' % fmax)

    Flow(n('comps'),n('freqs'),
         'rtoc | math output="exp(I*input*x1*%g)" ' % pi2 )

    Flow([n('cwht'),n('cfit')],[n('comps'),n('ctrace')],
         'clpf match=${SOURCES[1]} rect1=%d pred=${TARGETS[1]}' % rect2)

    Flow(n('cdif'),[n('cfit'),n('ctrace')],'add scale=1,-1 ${SOURCES[1]}')

    Result(n('fit'),[n('cdif'),n('cfit'),n('ctrace')],
           '''
           cat axis=2 ${SOURCES[1:3]} | real | 
           dots labels=Difference:Fit:Original gaineach=n
           ''')

    # mask low amplitudes

    Flow(n('mask'),n('cwht'),
         'math output="abs(input)" | real | mask min=%g' % amin)

    for emd in list:
        mask = n('mask%d' % emd)
        Flow(mask,n('mask'),'window n2=1 f2=%d' % emd)

        grup = n('group%d' % emd)
        Flow(grup,[n('group'),mask],
             '''
             window n2=1 f2=%d | 
             rtoc | math output="x1+I*input" | 
             transp plane=12 | headerwindow mask=${SOURCES[1]} | window
             ''' % emd)
        Plot(grup,
             '''
             graph title=Frequencies yreverse=y pad=n
             wanttitle=n scalebar=n bartype=v
             plotfat=5 plotcol=%d grid=%d label2=Frequency unit2=Hz
             min1=0 max1=%g min2=0 max2=150 
             ''' % (emd-1,emd==list[0],tmax))
        group.append(grup)

    Result(n('mgroup'),group,'Overlay')
