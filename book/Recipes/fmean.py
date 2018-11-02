from rsf.proj import *
import math

   
def Freqmean(file,r=10,std='n'):
   '''Mean and std of frequency for t-f data.
      r:   smoothing radius (for division).
      std: if 'y', return the std as well.'''
   
   Flow(file+'_s1', file, 
        'window min2=0. squeeze=n | math output="abs(input)*x2" | real | stack')
   Flow(file+'_s0',  file, 
        'window min2=0. squeeze=n | math output="abs(input)" | real | stack')
   Flow(file+'_fmn', '%s_s1 %s_s0' %(file,file), 
        'divn rect1=%d den=${SOURCES[1]} verb=n' %r)

   if std=='y':
      Flow(file+'_s2', file,
           'window min2=0. squeeze=n | math output="abs(input)*x2*x2" | real | stack')

      Flow(file+'_fstd','%s_fmn %s_s0 %s_s1 %s_s2' %(file,file,file,file),
           '''math s0=${SOURCES[1]} s1=${SOURCES[2]} s2=${SOURCES[3]} 
              output="s2-2*s1*input+s0*input*input" | 
              divn den=${SOURCES[1]} rect1=%d verb=n | clip2 lower=0 | 
              math output="sqrt(input)"''' %r)

def Freqplots(out, inp,axes=' ', scalebar='n',
              plotfat=10,plotcol=1, dash=1, transp='y', yreverse='y'):
   Plot(out,inp,
        '''graph wanttitle=n wantaxis=n scalebar=%s transp=%s yreverse=%s 
           %s plotfat=%d plotcol=%d dash=%d
        '''%(scalebar,transp,yreverse,axes,plotfat,plotcol,dash))

def iscmplx(file):
   Flow('%s.par'%file, file, 'get parform=n data_format')
   print('ahem')
   txt = open('%s.par'%file, 'r')
   line = txt.readlines()
   txt.close()
   
   if line[-7:]=='complex':
      return True
   else:
      return False



