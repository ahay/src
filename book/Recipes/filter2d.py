try:
    from rsf.cluster import *
except:
    from rsf.proj import *
import fdmod,math 

'''
    This module should help to do 2d filtering
    
    I will add in the future some filter for upgoing/downgoing
    wavefield separation.
    
    Esteban Diaz,

    Center for Wave Phenomena
'''


def fourier2d(four2d,fin,a1=1,a2=2):
    # Fourier transform along axis 1 and 2
    # Input is real, output is complex
    Flow(four2d,fin,
      '''
      rtoc |
      fft3 axis=%d | 
      fft3 axis=%d  
      '''%(a1,a2))

def invfourier2d(invfour2d,fin):
    # Fourier transform along axis 1 and 2
    # Input is complex (fourier domain), output is real
    Flow(four2d,fin,
      '''
      fft3 axis=1 inv=y| 
      fft3 axis=2 inv=y|
      real
      ''')

def spectra2d(spectra2d,fin):
    # 2D fourier spectra 
    # output is real
    Flow(spectra2d,fin,
      '''
      math output="sqrt(input*conj(input))" | 
      real |
      window min1=0 min2=0
      ''')

def plot2d (custom):
    return '''
    grey allpos=yes color=j  %s
    '''%custom

def f2dv1(filt,inp,slp=0.1,smooth="rect1=1 rect2=1 repeat=1",sym=1):
    # 2D filter based on spectra slopes,
    # slope=d1/d2  (change on axis1 over change on axis2
    # inp: complex file with the same grid of the fourier2d()
    # Default slope=0.1    
    # Filt is a complex file, with the same filter in
    # both real and imaginary parts
    # smooth: custom smoothing for tapering the filter
    tmp='tmp'
    if (sym==1) :
        Flow(filt,inp,
            '''
            real |
            math output="1*(%g*abs(x2) -abs(x1))/sqrt(%g*%g +1)" |
            mask min=0 | 
            dd type=float | 
            smooth %s |
            rtoc |
            math output="input -sqrt(-1)*input"
            '''%(slp,slp,slp,smooth))
    else:
        Flow(filt,inp,
            '''
            real |
            math output="1*(%g*x2 -x1)/sqrt(%g*%g +1)" |
            mask min=0 | 
            dd type=float | 
            smooth %s |
            rtoc |
            math output="input -sqrt(-1)*input"
            '''%(slp,slp,slp,smooth))
                                   


 
def f2dv1v2(filt,inp,slp1=0.1,slp2=0.4,smooth="rect1=1 rect2=1 repeat=1",sym=0,reject=1):
    # 2D filter based on spectra slopes,
    # slope=d1/d2  (change on axis1 over change on axis2
    # slp2 >slp1
    # fan filter
    # inp: complex file with the same grid of the fourier2d()
    # Default slope=0.1    
    # Filt is a complex file, with the same filter in
    # both real and imaginary parts
    # smooth: custom smoothing for tapering the filter
    tslp1=slp1
    tslp2=slp2
    if (abs(slp2)<abs(slp1)): 
       tslp2=slp1
       tslp1=slp2
    slp1=tslp1
    slp2=tslp2  
    f2dv1(filt+'tmp1',inp,slp1,smooth="rect1=1 rect2=1 repeat=1",sym=sym)
 
    f2dv1(filt+'tmp2',inp,slp2,smooth="rect1=1 rect2=1 repeat=1",sym=sym)

    if reject==1:
        Flow(filt,[filt+'tmp2',filt+'tmp1'],
          '''
          add mode=a scale=1,-1 ${SOURCES[1]}|
          real |
          math output="abs(input)" | 
          smooth %s |
          scale axis=123 |
          rtoc | 
          math output="input- sqrt(-1)*input"
          '''%smooth)
    else:
        Flow(filt,[filt+'tmp2',filt+'tmp1'],
          '''
          add mode=a scale=1,-1 ${SOURCES[1]}|
          real |
          math output="abs(input)" | 
          smooth %s |
          scale axis=123 |
          rtoc | 
          math output="1-(input- sqrt(-1)*input)"
          '''%smooth)
        


def kzfilt(filt,inp,smooth="rect1=1 rect2=1 repeat=1",forw=1.0):
    # Filter designed to keep positive kz
    forw=forw/abs(forw)

    Flow(filt,inp,
       '''
       real |
       math output="abs(sign(%g*x1)-1.0)*0.5" |
       smooth %s |
       rtoc |
       math output="input -sqrt(-1)*input"
       '''%(forw,smooth))


def kxfilt(filt,inp,smooth="rect1=1 rect2=1 repeat=1",forw=1.0):
    # Filter designed to keep positive kz
    forw=forw/abs(forw)

    Flow(filt,inp,
       '''
       real |
       math output="abs(sign(%g*x2)-1.0)*0.5" |
       smooth %s |
       rtoc |
       math output="input -sqrt(-1)*input"
       '''%(forw,smooth))


def highpass2d(filt,inp,smooth="rect1=1 rect2=1 repeat=1",cutoff=0.3):
    # Filter designed to keep positive kz

    cut=cutoff*math.sqrt(2)

    Flow(filt,inp,
       '''
       real |
       math output="sqrt((x1)^2 +x2^2)" |
       mask min=%g |
       dd type=float |
       smooth %s |
       rtoc |
       math output="input -sqrt(-1)*input"
       '''%(cut,smooth))

def lowpass2d(filt,inp,smooth="rect1=1 rect2=1 repeat=1",cutoff=0.3):
    # Filter designed to keep positive kz

    cut=cutoff*math.sqrt(2)

    Flow(filt,inp,
       '''
       real |
       math output="sqrt((x1)^2 +x2^2)" |
       mask min=%g |
       dd type=float |
       math output="1.0-input"|
       smooth %s |
       rtoc |
       math output="input -sqrt(-1)*input"
       '''%(cut,smooth))



def filtinvert(filtdata,fourier2d,filt2d):
    Flow(filtdata,[fourier2d,filt2d],
        '''
        add mode=p ${SOURCES[1]} |
        fft3 axis=1 inv=y |
        fft3 axis=2 inv=y | 
        real
        ''')

    

