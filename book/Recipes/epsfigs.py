import os

def Convert(src,clr):   
   src1,src2 = os.path.split(src)
   
   if os.path.isfile(src) and src[-4:]=='.vpl':
      trg = src[:-3]+'eps'
      print 'Converting Fig/%s to Fig/%s.eps...' %(src2,src2[:-4]),
                      
      fail = os.system( 'vpconvert %s color=%s %s\n'%(src,clr,trg) )
      if fail:
         print '\nFailed to convert "%s" to eps!\n' %src2 
         sys.exit(1)
      print ' Done.'
   else:
      print 'File "Fig/%s" not found' %src2



base = os.getcwd()

def Epsfigs(files=None,color=False):
   """Convert all vpl files in Fig to eps. This is useful, if you use Madagascar 
for the experiments but not for the latex compilation.
   
USAGE
   In the SConstruct file:
   
      from rsf.recipes.epsfigs import Epsfigs
      Epsfigs(files=None,color=False)
   
   If files is None, all the vpl files in Fig/ are converted.
   
NOTE: You need to run scons twice (cannot find dependencies).
   """
   if color=='y' or color==1 or color==True:
      clr='y'
   elif color=='n' or color==0 or color==False:
      clr='n'

   figdir  = os.path.join(base,'Fig')
   if not os.path.isdir(figdir):
      return None
   
   if files is None:      
      figlist = os.listdir(figdir)
      for fig in figlist:
         if len(fig)<=4 or fig[-4:]!='.vpl':
            figlist.remove(fig)
   else:
      figlist = files.split()
      for i in range(len(figlist)):
         figlist[i] += '.vpl'
          
   for fig in figlist:
      src = os.path.join(figdir,fig)
      Convert(src,clr)





