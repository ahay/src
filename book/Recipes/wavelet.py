try:    from rsf.cluster import *
except: from rsf.proj    import *
import wplot
import math

# ------------------------------------------------------------
def makeCarrier(fpeak,par):
    return '''
    math n1=%d o1=%g d1=%g
    output="cos(2*3.14*(%g)*(x1+(%g)))" |
    put label1=%s unit1='%s'
    '''%(par['nt'],par['ot'],par['dt'],
         fpeak, -par['kt']*par['dt'],  # peak freq
         par['lt'],par['ut'])

def gauEnvelope(fband,par):
    return '''
    math n1=%d o1=%g d1=%g
    output="exp(-2*3.14^2*(%g)^2*(x1+(%g))^2)" |
    put label1=%s unit1='%s'
    '''%(par['nt'],par['ot'],par['dt'],
         fband/4.0, -par['kt']*par['dt'], # freq std dev
         par['lt'],par['ut'])

def boxEnvelope(fband,par):
    return '''
    math n1=%d o1=%g d1=%g
    output="sin(3.14*(%g)*(x1-%g))/(3.14*(%g)*(x1-%g+1e-6))" |
    put label1=%s unit1='%s'
    '''%(par['nt'],par['ot'],par['dt'],
         fband, par['kt']*par['dt'],
         fband, par['kt']*par['dt'],
         par['lt'],par['ut'])

def makeWavelet(w,c,e,par):
    Flow(w,[c,e],'add mode=p ${SOURCES[1]}')

# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
def makeBoxSpec(s,fPEAK,fBAND,par):
    par['df'] = 0.5/par['tmax']
    par['of'] = 0.0
    par['nf'] = 2*par['nt']-1

    boxS = fBAND/(2*math.sqrt(3));
    boxL = fPEAK - math.sqrt(3) * boxS;
    boxH = fPEAK + math.sqrt(3) * boxS;
    kf = 1 + boxL/par['df'] 
    lf = 1 + boxH/par['df']
    
    Flow(s,None,
    'spike nsp=1 mag=0 n1=%(nf)d o1=%(of)g d1=%(df)g |'%par +          
    '''
    spike nsp=1 mag=1 k1=%d l1=%d |
    smooth rect1=5 | 
    scale axis=123 
    '''%(kf,lf))

def makeGauSpec(s,fPEAK,fBAND,par):
    par['df'] = 0.5/par['tmax']
    par['of'] = 0.0
    par['nf'] = 2*par['nt']-1

    gauS = fBAND/(2*math.sqrt(3));
    gauC = fPEAK;

    Flow(s,None,
    'spike nsp=1 mag=0 n1=%(nf)d o1=%(of)g d1=%(df)g |'%par +                 
    '''
    math output="exp(-0.5*((x1-%g)/%g)^2)" |
    scale axis=123
    '''%(gauC,gauS))

def makeTimeWavelet(t,f,par):
    Flow(t,f,
    '''
    rtoc | fft3d cnt=y axis=0 | real |
    window f1=%d n1=%d | put o1=%f d1=%f |
    scale axis=123
    '''%(par['nt']-par['kt'],par['nt'],par['ot'],par['dt']))
    
    
# ------------------------------------------------------------
def plotWavelet(plot,wvl,env,par):
    Plot(wvl,wplot.waveplot('plotcol=5',par))
    Plot(env,wplot.waveplot('plotcol=3',par))
    Result(plot,[wvl,env],'Overlay')
    
def plotSpectrum(plot,wvl,par):
    Flow(wvl+'spec',wvl,
        'spectra | put d1=%g'%(1./(par['nt']*par['dt'])))
    Result(plot,wvl+'spec',
        wplot.specplot('label1=%(lf)s unit1=%(uf)s plotfat=10 max1=30'%par,par))
