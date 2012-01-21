from rsf.proj import *
import math

def vlines(out,inp,max=1.,axes='',plot=''):
   ''' Plot with vertical lines on specific positions'''
   
   max = abs(max);
   Plot(out,inp,
        '''spray axis=2 n=2 o=%g d=%g | rtoc | math output="input+I*x2" | transp | 
           graph %s wanttitle=n wantaxis=n %s''' %(-max,2*max,axes,plot))


def hlines(out,inp,max=1.,axes='',plot=''):
   ''' Plot with vertical lines on specific positions'''
   
   max = abs(max);
   Plot(out,inp,
        '''spray axis=2 n=2 o=%g d=%g | rtoc | math output="x2 + I*input" | transp | 
           graph %s wanttitle=n wantaxis=n %s transpose=y''' %(-max,2*max,axes,plot))


# def hlines2(out,inp,axes:w

