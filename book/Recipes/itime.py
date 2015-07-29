from rsf.proj import *
import math


def itimec(out,data, cltf=' ',div=' '):
   '''Instantaneous traveltimes as a function of t and f (CLTFT).
      cltf: parameters for sfcltft.
      div:  parameters for smooth division.'''
   
   Flow(data+'_cltf',  data,'rtoc | cltft %s ' %cltf)
   Flow(data+'_dcltf', data,
       'rtoc | math output="input*x1*I" | cltft %s ' %cltf)

   Flow(out, [data+'_dcltf', data+'_cltf'], 'cdivn den=${SOURCES[1]} %s | imag'%div)


def itime1(out,data, ltft=' ',div=' '):
   '''Instantaneous traveltimes as a function of t and f (LTFT).
      ltft: parameters for sfltft.
      div:  parameters for smooth division.'''
   
   Flow(data+'_ltft1',  data,'ltft %s ' %ltft)
   Flow(data+'_dltft1', data,
       'math output="input*x1" | ltft %s | math output="input*I" ' %ltft)

   Flow(out, [data+'_dltft1', data+'_ltft1'], 'cdivn den=${SOURCES[1]} %s | imag'%div)


def itime2(out,data, ltft=' ',div=' '):
   '''Instantaneous traveltimes as a function of t and f (LTFT).
      ltft: parameters for sfltft.
      div:  parameters for smooth division.'''
   
   Flow(data+'_ltft2', data,'ltft %s ' %ltft)
   Flow(data+'_dltft2', data,
        '''math output="input*x1" | ltft %s | math output="I*input" ''' %ltft)
   Flow(out, [data+'_dltft2', data+'_ltft2'],
        'cdivn den=${SOURCES[1]} %s | imag' %div)


def itimef(out,data, ft=' ', div=' '):
   '''Instantaneous traveltimes as a function of t and f (FT)
      ft:  fft1 parameters.
      div: parameters for smooth division.'''
   
   Flow(data+'_ft',  data,'fft1 %s' %ft )
   Flow(data+'_dft', data,
        'math output="input*x1" | fft1 %s | math output="input*I"' %ft)
   Flow(out, [data+'_dft', data+'_ft'],
        'cdivn den=${SOURCES[1]} %s | imag | put label2="tau"' %div)


"""
def itime1(out,data, ltft=' ',div=' ', eps=0.2):
   '''Instantaneous traveltimes as a function of t and f (LTFT).
      ltf: parameters for sfltft.
      div: parameters for smooth division.
      eps: smoothness parameter for division.'''
   
   Flow(data+'_ltft1', data,'ltft %s ' %ltft)

   Flow(data+'_imag_der1', data+'_ltft1',
       'imag | transp plane=12 | smoothder eps=%g | transp plane=12 | rtoc' %eps)
   Flow(out, [data+'_ltft1', data+'_imag_der1'],
        '''real | transp plane=12 | smoothder eps=%g | transp plane=12 |
           rtoc | math img=${SOURCES[1]} output="input + I*img" |    
           cdivn den=${SOURCES[1]} %s | imag
        '''%(eps,div))

        #   cmplx ${SOURCES[1]} | cdivn den=${SOURCES[1]} %s | imag
        #'''%(eps,div))
"""


