import sys, os, string, re, commands, types, py_compile
from os.path import join
import SCons


from sfconfigure import ToolCreator

from SCons.Script import *

# CONSTANTS -- DO NOT CHANGE
context_success = 1
context_failure = 0
py_success = 0 # user-defined
unix_failure = 1


# Make sure error messages stand out visually
def stderr_write(message):
    sys.stderr.write('\t%s\n' % message)

def need_pkg(env,type,fatal=True):
    mypkg = env['PKG'][type].get( env['PLAT']['distro'])
    if mypkg:
        stderr_write('Optional package: ' + mypkg)
        
    if fatal:
        stderr_write('Fata Error: Needed package: ' + type)
        raise Exception
        sys.exit(unix_failure)

def CreateFalseTool( name, context ):
    newtool = ToolCreator( name, dest=context.env['tool_dest'])
    newtool.Exists(False)
    newtool.CreateTool( context.env )

def AllChecks( context ):
    return {
    'CheckCC': cc,
    'CheckCCDebug':ccdebug,
    'CheckAR':ar,
    'CheckLIBS':libs,
    'CheckC99':c99 ,
    'CheckX11':x11,
    'CheckPPM':ppm,
    'CheckJPEG':jpeg,
    'CheckOPENGL':opengl,
    'CheckBLAS':blas,
    'CheckMPI':mpi ,
    'CheckOMP':omp ,
    'CheckCXX':cxx,
    'CheckF77':f77,
    'CheckF90':f90,
    'CheckMATLAB':matlab,
    'CheckOCTAVE':octave,
    'CheckPYTHON':python,
    }
    


def check_all(context):
    

    # FDNSI = Failure Does Not Stop Installation
#    identify_platform(context)
    context.env.Tool('rsfroot')
    context.env.Tool('platform_ext')
    
    cctool = ToolCreator('rsfcc', dest=context.env['tool_dest'])
    cctool.Exists(True)
    
    cxxtool = ToolCreator('rsfcxx',dest=context.env['tool_dest'])
    f77tool = ToolCreator('rsff77',dest=context.env['tool_dest'])
    f90tool = ToolCreator('rsff90',dest=context.env['tool_dest'])
    mattool = ToolCreator('rsfmatlab',dest=context.env['tool_dest'])
    octtool = ToolCreator('rsfoctave',dest=context.env['tool_dest'])
    pyttool = ToolCreator('rsfpython',dest=context.env['tool_dest'])
    numpytool = ToolCreator('numpytool',dest=context.env['tool_dest'])
    
    context.env['RSFCC'] = cctool
    context.env['RSFCXX']  = cxxtool
    context.env['RSFF77']  = f77tool
    context.env['RSFF90']  = f90tool
    context.env['RSFMATLAB']  = mattool
    context.env['RSFOCTAVE']  = octtool
    context.env['RSFPYTHON']  = pyttool
    context.env['NUMPYTOOL']  = numpytool

    cc  (context)
    ccdebug( context )

    ar  (context)
    libs(context)
    c99 (context) # FDNSI
    
    x11 (context) # FDNSI
    ppm (context) # FDNSI
    jpeg(context) # FDNSI
    opengl(context) # FDNSI
    blas(context) # FDNSI
    mpi (context) # FDNSI
    omp (context) # FDNSI
#    api = api_options(context)
    
    disabled = GetOption('disabled')
    if 'c++' in disabled:
        cxxtool.Exists(False)
        context.Message("checking c++ API ... ")
        context.Result('no (explicitly disabled)')
    else:
        cxx(context)
        
    if 'f77' in disabled:
        f77tool.Exists(False)
        context.Message("checking f77 API ... ")
        context.Result('no (explicitly disabled)')
    else:
        f77(context)
        
    if 'f90' in disabled:
        f90tool.Exists(False)
        context.Message("checking f90 API ... ")
        context.Result('no (explicitly disabled)')

    else:
        f90(context)
        
    if 'matlab' in disabled:
        mattool.Exists(False)
        context.Message("checking matlab API ... ")
        context.Result('no (explicitly disabled)')

    else:
        matlab(context)
        
    if 'octave' in disabled:
        octtool.Exists(False)
        context.Message("checking octave API ... ")
        context.Result('no (explicitly disabled)')

    else:
        octave(context)
        
    if 'python' in disabled:
        pyttool.Exists(False)
        context.Message("checking python API ...")
        context.Result('no (explicitly disabled)')

    else:
        python(context)

    cctool.CreateTool( context.env )
    cxxtool.CreateTool( context.env )
    f77tool.CreateTool( context.env )
    f90tool.CreateTool( context.env )
    mattool.CreateTool( context.env )
    octtool.CreateTool( context.env )
    pyttool.CreateTool( context.env )
    numpytool.CreateTool( context.env )

def ccdebug(context):
    dbtool = ToolCreator('rsfcc_debug', dest=context.env['tool_dest'])
    
    CCFLAGS = context.env['CCFLAGS']
    CCFLAGS.append( '-g')
    
    text = '''
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    context.Message("checking if '-g' works ... ")
    
    res = context.TryLink(text,'.c')
    context.Result(res)
    if res:
        dbtool.Append(CCFLAGS='-g')
    
    dbtool.Exists(True)    
    
    dbtool.CreateTool( context.env )
    
    CCFLAGS.pop()

# A C compiler is needed by most Madagascar programs
# Failing this test stops the installation.
def cc(context):
    cctool = context.env['RSFCC']

    plat = context.env['PLAT']
    pkg = context.env['PKG']
    
    pkg['gcc'] = {'fedora':'gcc'}
    
    context.Message("checking for C compiler ... ")
    CC = context.env.get('CC',WhereIs('gcc'))
    if CC:
        context.Result(CC)
    else:
        context.Result(context_failure)
        need_pkg( context.env, 'gcc')
    text = '''
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    pkg['libc'] = {'fedora':'glibc',
                   'generic':'libc6-dev'}

    context.Message("checking if %s works ... " % CC)
    res = context.TryLink(text,'.c')
    context.Result(res)
    if not res:
        need_pkg( context.env, 'libc')
    if string.rfind(CC,'gcc') >= 0:
        
        oldflag = context.env.get('CCFLAGS',[])
        oldflag = Split( oldflag )
        context.env['CCFLAGS']= oldflag
#        print "oldflag",oldflag
        for flag in (['-std=gnu99','-Wall','-pedantic'],
                     ['-std=gnu9x', '-Wall','-pedantic'],
                     ['-Wall','-pedantic']):
            
            context.Message("checking if gcc accepts '%s' ... " % ' '.join(flag))
            
            context.env['CCFLAGS'] = oldflag + flag
            
            res = context.TryCompile(text,'.c')
            context.Result(res)
            if res:
                break
        if not res:
            context.env['CCFLAGS'] = oldflag
        else:
            
            cctool.Append( CCFLAGS=flag )
        # large file support
        (status,lfs) = commands.getstatusoutput('getconf LFS_CFLAGS')
        if not status and lfs:
            oldflag = context.env.get('CCFLAGS')
            context.Message("checking if gcc accepts '%s' ... " % lfs)
            context.env['CCFLAGS'] = oldflag + Split(lfs)
            res = context.TryCompile(text,'.c')
            context.Result(res)
            if not res:
                context.env['CCFLAGS'] = oldflag
            else:
                cctool.Append( CCFLAGS=lfs )
        # if Mac OS X and fink installed, update CPPPATH and LIBPATH
        if plat['OS'] == 'darwin' and os.path.isdir('/sw'):
            context.env['CPPPATH'] = context.env.Append(CPPPATH=['/sw/include',])
            
            context.env['LIBPATH'] = context.env.get('LIBPATH',[]) + \
                                         ['/sw/lib',]
                                         
            context.env['LINKFLAGS'] = context.env.get('LINKFLAGS','') + \
                                         ' -framework Accelerate'
            
            cctool.Append( CPPPATH=['/sw/include'] ,LIBPATH=['/sw/lib'],
                           FRAMEWORKS=['Accelerate'])
#                           LINKFLAGS=' -framework Accelerate')
            
    elif plat['OS'] == 'sunos':
        if '-O2' in context.env.get('CCFLAGS',[]):
            context.env['CCFLAGS'].remove( '-O2' )
            context.env.Append(CCFLAGS=['-xO2'])
            cctool.Append( CCFLAGS=['-x02'] )
            
    # remove duplicates
    CCFLAGS = dict.fromkeys( context.env.get('CCFLAGS',[]) ).keys()
    context.env['CCFLAGS'] = CCFLAGS
    cctool.UpdateExists(True)
    

# Used for building libraries.
def ar(context):
    pkg = context.env['PKG']

    pkg['ar']={'fedora':'binutils'}

    context.Message("checking for ar ... ")
    AR = context.env.get('AR',WhereIs('ar'))
    if AR:
        context.Result(AR)
        context.env['AR'] = AR
    else:
        context.Result(context_failure)
        need_pkg( context.env, 'ar')
        
    cctool = context.env['RSFCC']
    cctool.UpdateExists(bool(AR))

# Failing this check stops the installation.
def libs(context):
    plat = context.env['PLAT']

    cctool = context.env['RSFCC']

    context.Message("checking for libraries ... ")
#    LIBS = context.env.get('LIBS','m')
#    if type(LIBS) is not types.ListType:
#        LIBS = string.split(LIBS)
    LIB = None
    
    if plat['OS'] == 'sunos' or plat['OS'] == 'hp-ux' or plat['OS'] == 'hpux':
        LIB = 'nsl'
    elif plat['OS'] == 'cygwin':
        LIBS.append('rpc')
    elif plat['OS'] == 'darwin':
        LIB = 'mx'
    elif plat['OS'] == 'interix':
        LIB = 'rpclib'
        
    text = '''
    #include <rpc/types.h>
    #include <rpc/xdr.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''
    
    
    if LIB:
        LIBS = [LIB,'m']
    else:
        LIBS = ['m']

    cctool.Append( LIBS=LIBS )
    context.env.Append( LIBS=LIBS )

    res = context.TryLink(text,'.c')
    if res:
        context.Result(str(LIB))
        context.env['LIBS'] = LIBS
        
        cctool.UpdateExists(True)
        
    else:
        context.Result(context_failure)
        need_pkg( context.env, 'libs')
        
        cctool.UpdateExists(False)


# Complex number support according to ISO C99 standard
def c99(context):
    
    cctool = context.env['RSFCC']

    context.Message("checking complex support ... ")
    text = '''
    #include <complex.h>
    #include <math.h>
    int main(int argc,char* argv[]) {
    float complex c;
    float f;
    f = cabsf(ccosf(c));
    return (int) f;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
    else:
        context.env.Append(CCFLAGS=['-DNO_COMPLEX'] )
        cctool.Append(CCFLAGS=['-DNO_COMPLEX'] )
        context.Result(context_failure)
        need_pkg( context.env, 'c99', fatal=False)

# The two lists below only used in the x11 check
xinc = [
    '/usr/X11/include',
    '/usr/X11R6/include',
    '/usr/X11R5/include',
    '/usr/X11R4/include',
    '/usr/include/X11',
    '/usr/include/X11R6',
    '/usr/include/X11R5',
    '/usr/include/X11R4',
    '/usr/local/X11/include',
    '/usr/local/X11R6/include',
    '/usr/local/X11R5/include',
    '/usr/local/X11R4/include',
    '/usr/local/include/X11',
    '/usr/local/include/X11R6',
    '/usr/local/include/X11R5',
    '/usr/local/include/X11R4',
    '/usr/X386/include',
    '/usr/x386/include',
    '/usr/XFree86/include/X11',
    '/usr/include',
    '/usr/local/include',
    '/usr/unsupported/include',
    '/usr/athena/include',
    '/usr/local/x11r5/include',
    '/usr/lpp/Xamples/include',
    '/usr/openwin/include',
    '/usr/openwin/share/include'
    ]

xlib = [
    '/usr/X11/lib64',
    '/usr/X11/lib',
    '/usr/X11R6/lib64',
    '/usr/X11R6/lib',
    '/usr/X11R5/lib',
    '/usr/X11R4/lib',
    '/usr/lib/X11',
    '/usr/lib/X11R6',
    '/usr/lib/X11R5',
    '/usr/lib/X11R4',
    '/usr/local/X11/lib',
    '/usr/local/X11R6/lib',
    '/usr/local/X11R5/lib',
    '/usr/local/X11R4/lib',
    '/usr/local/lib/X11',
    '/usr/local/lib/X11R6',
    '/usr/local/lib/X11R5',
    '/usr/local/lib/X11R4',
    '/usr/X386/lib',
    '/usr/x386/lib',
    '/usr/XFree86/lib/X11',
    '/usr/lib',
    '/usr/local/lib',
    '/usr/unsupported/lib',
    '/usr/athena/lib',
    '/usr/local/x11r5/lib',
    '/usr/lpp/Xamples/lib',
    '/lib/usr/lib/X11',
    '/usr/openwin/lib',
    '/usr/openwin/share/lib'
    ]

# If this check is failed
# you may not be able to display .vpl images on the screen
def x11(context):
    text = '''
    #include <X11/Intrinsic.h>
    #include <X11/Xaw/Label.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''
    plat = context.env['PLAT']
    
    x11tool  = ToolCreator('x11', dest=context.env['tool_dest'])
    x11tool.Exists(True)
    
    context.Message("checking for X11 headers ... ")
    INC = context.env.get('XINC','')
    if type(INC) is not types.ListType:
        INC = string.split(INC)

    oldpath = context.env.get('CPPPATH',[])

    if not oldpath:
        oldpath = []        
    res = None
    for path in filter(lambda x:
                       os.path.isfile(os.path.join(x,'X11/Xaw/Label.h')),
                       INC+xinc):
#        import pdb;
#        pdb.set_trace()
        context.env['CPPPATH'] = oldpath + [path,] 
        res = context.TryCompile(text,'.c')

        if res:
            context.Result(path)
            context.env['XINC'] = [path,]
            break

    if not res:
        context.Result(context_failure)
        x11tool.Exists(False)
        stderr_write('xtpen (for displaying .vpl images) will not be built.')
        need_pkg( context.env, 'xaw', fatal=False)
        context.env['XINC'] = None
        return
    else:
        x11tool.Replace( XINC=[path])
        x11tool.Append( CPPPATH=[path])

    context.Message("checking for X11 libraries ... ")
    LIB = context.env.get('XLIBPATH','')
    if type(LIB) is not types.ListType:
        LIB = string.split(LIB)

    oldlibpath = context.env.get('LIBPATH',[])
    oldlibs = context.env.get('LIBS',[])

    XLIBS = context.env.get('XLIBS')
    if XLIBS:
        if type(XLIBS) is not types.ListType:
            XLIBS = string.split(XLIBS)
    else:
        if  plat['OS'] == 'interix':
            XLIBS =  ['Xaw','Xt','Xmu','X11','Xext','SM','ICE']
        elif plat['OS'] == 'linux' or plat['OS'] == 'posix':
            XLIBS = ['Xaw','Xt']
        else:
            XLIBS = ['Xaw','Xt','X11']

    res = None
    for path in filter(os.path.isdir,LIB+xlib):
        context.env['LIBPATH'] = oldlibpath + [path,] 
        res = context.TryLink(text,'.c')

        if res:
            context.Result(path)
            context.env['XLIBPATH'] = [path,]
            context.env['XLIBS'] = XLIBS
            break
    if not res:
        x11tool.Exists(False)
        context.Result(context_failure)
        context.env['XLIBPATH'] = None
    else:
        x11tool.Replace(XLIBPATH = [path] )
        x11tool.Replace(XLIBS = XLIBS )
        x11tool.Append( LIBS = XLIBS )
        x11tool.Append( LIBPATH = [path] )
#        context.env.Append( LIBS=XLIBS )
        

    context.env['CPPPATH'] = oldpath
    context.env['LIBPATH'] = oldlibpath
    context.env['LIBS'] = oldlibs
    
    x11tool.CreateTool( context.env )
    


def ppm(context):
    
    ppmtool  = ToolCreator('ppm', dest=context.env['tool_dest'])
    ppmtool.Exists(True)

    context.Message("checking for ppm ... ")
    LIBS = context.env.get('LIBS',['m'])
#    if type(LIBS) is not types.ListType:
    LIBS = Split(LIBS)

    text = '''
    #include <ppm.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''
    
    for ppm in [context.env.get('PPM','netpbm'),'netpbm.10']:
        
	LIBS.append(ppm)
	res = context.TryLink(text,'.c')
	
	if res:
	    context.Result(res)
	    context.env['PPM'] = ppm
        
	    break
	else:

	    LIBS.pop()

    if res:
        ppmtool.Append(LIBS=ppm)
        ppmtool.Exists(True)
        ppmtool.Replace( PPM=ppm)
        LIBS.pop()
    else:
        
        ppmtool.Exists(False)
        ppmtool.Replace( PPM=None)
        
        context.Result(context_failure)
        need_pkg( context.env, 'netpbm', fatal=False)
        context.env['PPM'] = None

    ppmtool.CreateTool( context.env )

# If this test is failed, no writing to jpeg files
def jpeg(context):
    
    jpgtool  = ToolCreator('jpeg', dest=context.env['tool_dest'])
    jpgtool.Exists(True)

    context.Message("checking for jpeg ... ")
    
    LIBS = context.env.get('LIBS',['m'])
    LIBS = Split(LIBS)
    
    jpeg = context.env.get('JPEG','jpeg')
    LIBS.append(jpeg)
    text = '''
    #include <stdio.h>
    #include <jpeglib.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['JPEG'] = jpeg
        jpgtool.Exists(True)
        jpgtool.Replace( JPEG=jpeg)
        jpgtool.Append( LIBS=jpeg)
        
    else:
        context.Result(context_failure)
        stderr_write('sfbyte2jpg will not be built.')
        need_pkg( context.env, 'jpeg', fatal=False)
        context.env['JPEG'] = None
        jpgtool.Exists(False)
        jpgtool.Replace( JPEG=None)

    LIBS.pop()
    jpgtool.CreateTool( context.env )


# If this test is failed, no opengl programs
def opengl(context):
    
    
    plat = context.env['PLAT']
    pkg = context.env['PKG']
    
    pkg['opengl'] = {'generic':'mesa-libGL-devel',
                     'fedora': 'mesa-libGL-devel + freeglut + freeglut-devel'}

    ogltool  = ToolCreator('opengl', dest=context.env['tool_dest'])
#    ogltool.Exists(True)

    context.Message("checking for OpenGL ... ")
    
    LIBS = context.env.get('LIBS',['m'])
    LIBS = Split(LIBS)
    
    LINKFLAGS = context.env.get('LINKFLAGS',[])
    
    if plat['OS'] == 'darwin':
        
        FRAMEWORKS = ['AGL','OpenGL','GLUT']
        context.env.Append(FRAMEWORKS=FRAMEWORKS)
#        oglflags = [' -framework AGL -framework OpenGL -framework GLUT']
#        context.env['LINKFLAGS'] = LINKFLAGS + oglflags
        ogl = []
        oglflags = []
    else:
        oglflags = None
        ogl = context.env.get('OPENGL')
        if ogl:
            ogl = Split(ogl)
        else:
            ogl = []
    context.env['LIBS'] = LIBS + ogl

    text = '''
    #ifdef __APPLE__
    #include <OpenGL/glu.h>
    #include <GLUT/glut.h>
    #else
    #include <GL/glu.h>
    #include <GL/glut.h>
    #endif
    int main(int argc,char* argv[]) {
    glutInit (&argc, argv);
    return 0;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)    
        context.env['OPENGL'] = ogl 
        context.env['OPENGLFLAGS'] = oglflags
        ogltool.Exists(True)
        ogltool.Append( LIBS=ogl )
        ogltool.Append( LINKFLAGS=oglflags )
        
        ogltool.Replace( OPENGL=ogl )
        ogltool.Replace( OPENGLFLAGS=oglflags )
    else:
        context.env['OPENGL'] = None 
        context.env['OPENGLFLAGS'] = None
        context.Result(context_failure)
        need_pkg( context.env, 'opengl', fatal=False)

        ogltool.Exists(False)
        
        ogltool.Replace( OPENGL=None )
        ogltool.Replace( OPENGLFLAGS=None )

    if res:
        glew(context,LIBS,ogl)
    else:
        glewtool  = ToolCreator('glew', dest=context.env['tool_dest'])
        glewtool.Exists(False)
        glewtool.CreateTool( context.env )

    context.env['LIBS'] = LIBS
    context.env['LINKFLAGS'] = LINKFLAGS
    
    ogltool.CreateTool( context.env )


# If this test is failed, no GLEW programs
def glew(context,LIBS,ogl):
    context.Message("checking for GLEW ... ")

    text = '''
    #include <GL/glew.h>
    #ifdef __APPLE__
    #include <GLUT/glut.h>
    #else
    #include <GL/glut.h>
    #endif
    int main(int argc,char* argv[]) {
    GLenum err;
    glutInit(&argc, argv);
    err = glewInit();
    return 0;
    }\n'''

    glewtool  = ToolCreator('glew', dest=context.env['tool_dest'])

    GLEW = context.env.get('GLEW','GLEW')
    context.env['LIBS'] =  LIBS + [GLEW] + ogl 
        
    res = context.TryLink(text,'.c')

    if res:
        glewtool.Exists(True)
        glewtool.Replace( GLEW=GLEW )
        glewtool.Append( LIBS=GLEW )
        context.Result(res)
        context.env['GLEW'] = GLEW
    else:
        glewtool.Exists(False)
        context.Result(context_failure)
        need_pkg( context.env, 'glew', fatal=False)

    glewtool.CreateTool( context.env )

def blas(context):
    
    plat = context.env['PLAT']
    pkg = context.env['PKG']
    pkg['blas'] = {'fedora':'blas + blas-devel + atlas + atlas-devel'}
    
    blastool  = ToolCreator('blas', dest=context.env['tool_dest'])

    context.Message("checking for BLAS ... ")
    
    LIBS = context.env.get('LIBS',['m'])
    LIBS = Split(LIBS)
    
    blas = context.env.get('BLAS','blas')
    
    LIBS.append(blas)
    text = '''
    #ifdef __APPLE__
    #include <vecLib/vBLAS.h>
    #else
    #include <cblas.h>
    #endif
    int main(int argc,char* argv[]) {
    float d, x[]={1.,2.,3.}, y[]={3.,2.,1.};
    d = cblas_sdot(3,x,1,y,1);
    return 0;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['LIBS'] = LIBS
        context.env['BLAS'] = blas

        context.env['RSFCC'].Append( LIBS=blas )
        context.env['RSFCXX'].Append( LIBS=blas )
        context.env['RSFMATLAB'].Append( LIBS=blas )
        context.env['RSFF77'].Append( LIBS=blas )
        context.env['RSFF90'].Append( LIBS=blas )

        blastool.Append( LIBS=blas)
        blastool.Replace( BLAS=blas)
        
        if plat['OS'] == 'cygwin':
            context.env['ENV']['PATH'] = context.env['ENV']['PATH'] + \
                                         ':/lib/lapack'
    else:
        # some systems require cblas and atlas
        LIBS.append('cblas')
        LIBS.append('atlas')
        res = context.TryLink(text,'.c')
        
        if res:
            blastool.Append( LIBS=[blas,'cblas','atlas'])
            blastool.Replace( BLAS=blas)
            blastool.Exists(True)
            context.Result(res)
            context.env['LIBS'] = LIBS
            context.env['BLAS'] = blas
        else:
            context.Result(context_failure)
            context.env.Append(CCFLAGS='-DNO_BLAS')
            context.env.Append(CXXFLAGS= ['-DNO_BLAS'] )
            LIBS.pop()
            LIBS.pop()
            LIBS.pop()
            context.env['BLAS'] = None
            need_pkg( context.env, 'blas', fatal=False)

            context.env['RSFCC'].Append( CPPDEFINES='NO_BLAS' )
            context.env['RSFCXX'].Append( CPPDEFINES='NO_BLAS' )
            context.env['RSFMATLAB'].Append( CPPDEFINES='NO_BLAS' )
            context.env['RSFF77'].Append( CPPDEFINES='NO_BLAS' )
            context.env['RSFF90'].Append( CPPDEFINES='NO_BLAS' )

            blastool.Exists( False )
    
    blastool.CreateTool( context.env )

def mpi(context):
    pkg = context.env['PKG']
    pkg['mpi'] = {'fedora':'openmpi, openmpi-devel, openmpi-libs'}
    
    
    mpitool  = ToolCreator('mpi', dest=context.env['tool_dest'])

    context.Message("checking for MPI ... ")
    mpicc = context.env.get('MPICC',WhereIs('mpicc'))
    if mpicc:
        context.Result(mpicc)
        context.Message("checking if %s works ... " % mpicc)
        # Try linking with mpicc instead of cc
        text = '''
        #include <mpi.h>
        int main(int argc,char* argv[]) {
        MPI_Init(&argc,&argv);
        MPI_Finalize();
        }\n'''
        cc = context.env.get('CC')
        context.env['CC'] = mpicc
        res = context.TryLink(text,'.c')
        context.env['CC'] = cc
    else: # mpicc not found
        context.Result(context_failure)
        res = None
    if res:
        context.Result(res)
        context.env['MPICC'] = mpicc
        
        mpitool.Exists(True)
        mpitool.Append(CC=mpicc)
        
    else:
        context.Result(context_failure)
        need_pkg( context.env, 'mpi', fatal=False)
        context.env['MPICC'] = None

        mpitool.Exists(False)
    
    mpitool.CreateTool( context.env )

def ncpus(env):
    'Detects number of CPUs'
    plat = env['PLAT']
    if plat['OS'] in ('linux','posix'):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            nr_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if type(nr_cpus) is int:
                if nr_cpus > 0:
                    return nr_cpus
    elif plat['OS'] == 'darwin':
        nr_cpus = int(os.popen2('sysctl -n hw.ncpu')[1].read())
        if nr_cpus > 0:
            return nr_cpus
    return 2 # default number of processors


def omp(context):
    
    omptool  = ToolCreator('omp', dest=context.env['tool_dest'])

    pkg = context.env['PKG']
    pkg['omp'] = {'fedora':'libgomp'}
    
    if ncpus(context.env) == 1:
        context.env['OMP'] = False
        omptool.Exists(False)
        omptool.CreateTool(context.env)
        return # only 1 cpu. OMP not needed
    
    context.Message("checking for OpenMP ... ")
    LIBS  = context.env.get('LIBS',[])
    CC    = context.env.get('CC','gcc')
    flags = context.env.get('CCFLAGS','')
    gcc = (string.rfind(CC,'gcc') >= 0)
    icc = (string.rfind(CC,'icc') >= 0)
    if gcc:
        LIBS.append('gomp')
        CCFLAGS = flags + [' -fopenmp']
    elif icc:
        LIBS.append('guide')
        LIBS.append('pthread')
        CCFLAGS = flags + [' -openmp', '-D_OPENMP']
    else:
        CCFLAGS = flags

    text = '''
    #include <omp.h>
    int main(void) {
    int nt;
#pragma omp parallel
{
    nt = omp_get_num_threads();
}
    return 0;
    }
    '''
    context.env['LIBS'] = LIBS
    context.env['CCFLAGS'] = CCFLAGS
    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['OMP'] = True
    else:
        context.Result(context_failure)
        need_pkg( context.env, 'omp', fatal=False)
        if gcc:
            LIBS.pop()
        if icc:
            LIBS.pop()
            LIBS.pop()
        context.env['LIBS'] = LIBS
        context.env['CCFLAGS'] = flags
        context.env['OMP'] = False

def api_options(context):
    
    cppapitool  = ToolCreator('cppapi', dest=context.env['tool_dest'])
    f77apitool  = ToolCreator('f77api', dest=context.env['tool_dest'])
    f90apitool  = ToolCreator('f90api', dest=context.env['tool_dest'])
    pythonapitool  = ToolCreator('pythonapi', dest=context.env['tool_dest'])
    matlabapitool  = ToolCreator('matlabapi', dest=context.env['tool_dest'])
    octaveapitool  = ToolCreator('octaveapi', dest=context.env['tool_dest'])

    context.Message("checking API options ... ")
    api = context.env.get('API')
    if api:
        api = Split(api)
        api = map(string.lower, api)
    else:
        api = []

    valid_api_options = ['','c++', 'fortran', 'f77', 'fortran-90',
                         'f90', 'python', 'matlab', 'octave']

    for option in api:
        if not option in valid_api_options:
            api.remove(option)

    # Make tests for fortrans in API easy
    for i in range(len(api)):
        if api[i] == 'fortran':
            api[i] = 'f77'
        elif api[i] == 'fortran-90':
            api[i] = 'f90'

    # Eliminate duplicates if user was redundant
    try: # sets module was introduced in Py 2.3
        import sets
        api = list(sets.Set(api))
        del sets
    except:
        pass # Not a big deal if this is not done

    # Improve output readability
    if api == ['']:
        context.Result('none')
    else:
        context.Result(str(api))
    context.env['API'] = api

    if 'c++' in api:
        cppapitool.Exists(True)
    if 'f77' in api:
        f77apitool.Exists(True)
    if 'f90' in api:
        f90apitool.Exists(True)
    if 'python' in api:
        pythonapitool.Exists(True)
    if 'matlab' in api:
        matlabapitool.Exists(True)
    if 'octave' in api:
        octaveapitool.Exists(True)
    
    cppapitool.CreateTool(context.env)
    f77apitool.CreateTool(context.env)
    f90apitool.CreateTool(context.env)
    pythonapitool.CreateTool(context.env)
    matlabapitool.CreateTool(context.env)
    octaveapitool.CreateTool(context.env)

    return api


# For the C++ API
def cxx(context):
    
    cxxtool = context.env['RSFCXX']
    
    context.Message("checking for C++ compiler ... ")
    CXX = context.env.get('CXX')
    
    if CXX:
        cxxtool.Exists(True)
        cxxtool.Replace( CXX=CXX )
        
        
        context.Result(CXX)
    else:
        cxxtool.Exists(False)
        context.Result(context_failure)
        need_pkg( context.env, 'c++')
        
    oldflag = context.env.get('CXXFLAGS',[])
#    import pdb; pdb.set_trace()
    
    oldflag = Split(oldflag)
    context.env['CXXFLAGS'] = oldflag
    
    context.Message("checking if %s works ... " % CXX)
    text = '''
    #include <valarray>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''
    res = context.TryLink(text,'.cc')
    context.Result(res)
    if not res:
        del context.env['CXX']
        sys.exit(unix_failure)
        
    if CXX[-3:]=='g++':
        for flag in [['-Wall','-pedantic'],['-pedantic']]:
            context.Message("checking if g++ accepts '%s' ... " % 
                            ' '.join(flag))
            context.env['CXXFLAGS'] = oldflag + flag
            res = context.TryCompile(text,'.cc')
            context.Result(res)
            if res:
                cxxtool.Append(CXXFLAGS=flag )
                break
            
        if not res:
            context.env['CXXFLAGS'] = oldflag
    
    CXXFLAGS = dict.fromkeys( context.env.get('CXXFLAGS',[]) ).keys()
    context.env['CXXFLAGS'] = CXXFLAGS

    cxxtool.CreateTool( context.env )
    
# Used in checks for both f77 and f90
fortran = {'g77':'f2cFortran',
           'f77':'f2cFortran',
           'gfortran':'NAGf90Fortran',
           'gfc':'NAGf90Fortran',
           'f2c':'f2cFortran'}

def f77(context):
    
#    import pdb;pdb.set_trace()
    f77tool  = context.env['RSFF77']

    context.Message("checking for F77 compiler ... ")
    F77 = context.env.get('F77')
    if not F77:
        compilers = ['gfortran','g77','f77','f90','f95','xlf90','pgf90',
                     'ifort','ifc','pghpf','gfc']
        F77 = context.env.Detect(compilers)
        if not F77:
            for comp in compilers:
                F77 = WhereIs(comp)
                if F77:
                    break
        context.env['F77'] = F77
        
        if F77:
            f77tool.Replace(F77=F77)
            context.Result(F77)
        else:
            context.Result(context_failure)
            need_pkg( context.env, 'f77',False)
            
            f77tool.Exists(False)
            return
    else:
        context.Result(F77)
                
    if os.path.basename(F77) == 'ifc' or os.path.basename(F77) == 'ifort':
        intel(context)
        context.env.Append(F77FLAGS=' -Vaxlib')
        f77tool.Append(F77FLAGS=' -Vaxlib')
    
    text = '''      program Test
      stop
      end
      '''
      
    context.Message("checking if %s works ... " % F77)
    oldlink = context.env.get('LINK')
    context.env['LINK'] = F77
    res = context.TryLink(text,'.f')
    context.env['LINK'] = oldlink
    context.Result(res)
    if not res:
        del context.env['F77']
#        sys.exit(unix_failure)
        f77tool.Exists(False)
#        f77tool.CreateTool( context.env)
        return
    else:
        f77tool.Exists(True)
        f77tool.Replace(LINK=F77)
        
    F77base = os.path.basename(F77)
    if F77base[:3] == 'f77' and plat['OS'] == 'sunos':
        cfortran = 'sunFortran'
    else:
        cfortran = fortran.get(F77base,'NAGf90Fortran')
    context.env['CFORTRAN'] = cfortran
    f77tool.Replace( CFORTRAN=cfortran )
    context.Message("checking %s type ... " % F77)
    context.Result(cfortran)
    
#    f77tool.CreateTool( context.env)

def f90_write_autofile(extension,content):
    filename=os.path.join('api','f90','ptr_sz.'+extension)
    handle = open(filename,'w')
    handle.write(content)
    handle.close()

def f90_write_ptr_sz(plat):
    'Tell poor Fortran whether platform is 32 or 64 bit'
    if plat['arch'] == '32bit':
        str_insert_f90 = '9'
        str_insert_c   = '32'
    elif plat['arch'] == '64bit':
        str_insert_f90 = '12'
        str_insert_c   = '64'
    msg = 'File created by config'
    f90_write_autofile('h',
        '/* %s */\n#define RSF%sBIT\n' % (msg,str_insert_c) )
    f90_write_autofile('f90',
        '! %s\ninteger, parameter :: PTRKIND=selected_int_kind(%s)\n' % (msg,str_insert_f90) )

def f90(context):
    f90tool = context.env['RSFF90']
    
    context.Message("checking for F90 compiler ... ")
    F90 = context.env.get('F90')
    if not F90:
        compilers = ['gfortran','gfc','f90','f95','xlf90','pgf90',
                     'ifort','ifc','pghpf']
        F90 = context.env.Detect(compilers)
        if not F90:
            for comp in compilers:
                F90 = WhereIs(comp)
                if F90:
                    break
        context.env['F90'] = F90
        if F90:
            f90tool.Replace(F90=F90,LINK=F90)
            context.Result(F90)
        else:
            
            context.Result(context_failure)
            need_pkg( context.env, 'f90',False)
            
            f90tool.Exists(False)
            return
        
    if os.path.basename(F90) == 'ifc' or os.path.basename(F90) == 'ifort':
        intel(context)
        context.env.Append(F90FLAGS=' -Vaxlib')
        f90tool.Append(F90FLAGS=' -Vaxlib')

    main = '''program Test
    end program Test
    '''
    module = '''module testf90
    end module testf90
    '''
    context.Message("checking if %s works ... " % F90)
    
    oldlink = context.env.get('LINK')
    context.env['LINK'] = F90
    res1 = context.TryCompile(module,'.f90')
    res2 = context.TryLink(main,'.f90')
    
    context.env['LINK'] = oldlink
    
    context.Result(res1 and res2)
    
    if not res1 or not res2:
        del context.env['F90']
        f90tool.Exists(False)
#        f90tool.CreateTool(context.env)
        
        return
    
        
    base = os.path.basename(F90)
    context.Message("checking %s type ... " % base)
    cfortran = fortran.get(base,'NAGf90Fortran')
    context.env['CFORTRAN90'] = cfortran
    
    f90tool.Replace( CFORTRAN90=cfortran)

    context.Result(cfortran)
    context.Message("checking F90 module extension ... ")
    f90module = re.compile(r'(?:testf90|TESTF90)(\.\w+)$')
    suffix = ''
    here = os.getcwd()
    for file in os.listdir(here):
        gotit = f90module.match(file)
        if gotit:
            suffix = gotit.group(1)
            os.remove(file)
            break
        
    context.env['F90MODSUFFIX'] = suffix
    
    f90tool.Replace( F90MODSUFFIX=suffix)
    
    context.Result(suffix)
    f90_write_ptr_sz(context.env['PLAT'])
    f90tool.Exists(True)
#    f90tool.CreateTool(context.env)

def matlab(context):
    
#    mtool  = ToolCreator('rsfmatlab', dest=context.env['tool_dest'])
    mtool  = context.env['RSFMATLAB']
    
    context.Message("checking for Matlab ... ")
    matlab = WhereIs('matlab')
    if matlab:
        context.Result(matlab)
        MATLABPATH = os.environ.get('MATLABPATH')
        if MATLABPATH:
            MATLABPATH += ':${lib_prefix}'
        else:
            MATLABPATH = '${lib_prefix}'
            
        MATCOM = 'MATLABPATH=%s %s -nosplash -nojvm -nodesktop' % \
            (MATLABPATH,matlab)
        context.env['MATLAB'] = MATCOM
        mtool.Replace(MATLAB=MATCOM)
    else:
        context.Result(context_failure)
#        stderr_write('Please install Matlab.')
        context.env['MATLAB'] = None
        mtool.Replace(MATLAB=None)
        mtool.Exists(False)
#        mtool.CreateTool(context.env)
        return
#        sys.exit(unix_failure)

    context.Message("checking for mex ... ")
    mex = WhereIs('mex')
    if mex:
        context.Result(mex)
        context.env['MEX'] = mex
        mtool.Replace(MEX=mex)
    else:
        context.Result(context_failure)
        mtool.Replace(MEX=None)
        context.env['MEX'] = None
        
        mtool.Exists(False)
#        mtool.CreateTool(context.env)
        
        return

    # See http://www.mathworks.com/access/helpdesk/help/techdoc/ref/mex.html
    plat = context.env['PLAT']
    
    if plat['OS'] == 'linux' or plat['OS'] == 'posix':
        if plat['arch'] == '32bit':
            suffix = 'glx'
        else:
            suffix = 'a64'
    elif plat['OS'] == 'sunos':
        suffix = 'sol'
    elif plat['OS'] == 'darwin':
        if plat['cpu'] == 'i386':
            suffix = 'maci'
        else:
            suffix = 'mac'
    else:
        suffix = 'glx'
    context.env['MEXSUFFIX'] = '.mex' + suffix

    mtool.Replace(MEXSUFFIX='.mex' + suffix)

#    mtool.CreateTool(context.env)

def octave(context):
    
#    otool  = ToolCreator('rsfoctave', dest=context.env['tool_dest'])
    otool = context.env['RSFOCTAVE']
    
    context.Message("checking for Octave ... ")
    octave = WhereIs('octave')
    if octave:
        context.Result(octave)
        context.env['OCTAVE'] = octave
        otool.Replace(OCTAVE=octave)
        
        context.Message("checking for mkoctfile ... ")
        mkoctfile = WhereIs('mkoctfile')
        
        if mkoctfile:
            context.Result(mkoctfile)
            context.env['MKOCTFILE'] = mkoctfile
            otool.Replace(MKOCTFILE=mkoctfile)
            otool.Exists(True)
        else:
            context.Result(context_failure)
            stderr_write('Please install mkoctfile.')
            need_pkg( context.env, 'mkoctfile',False)
            otool.Exists(False)
            
    else: # octave not found
        context.Result(context_failure)
#        stderr_write('Please install Octave.')
        need_pkg( context.env, 'octave',False)
        otool.Exists(False)


#    otool.CreateTool(context.env)

def python(context):
    
#    pytool  = ToolCreator('rsfpython', dest=context.env['tool_dest'])
    pytool = context.env['RSFPYTHON']

    context.Message("checking for SWIG ... ")
    
    if not Tool('swig').exists(context.env):
        pytool.Exists(False)
#        pytool.CreateTool( context.env )
        context.Result( False )
        return
    else:
        context.Result( 'yes' )
    
    context.Message("checking for numpy ... ")
    
    numpytool = context.env['NUMPYTOOL'] 
    
    try:
        import numpy
    except ImportError:
        
        context.Result(context_failure)
        need_pkg( context.env, 'numpy',False)

        pytool.Exists(False)
#        pytool.CreateTool( context.env )
        context.Result( False )
        numpytool.Exists( False )
        return

    else:
        context.Result(context_success)
        context.env['PYMODULES'] = ['numpy']
        pytool.Append(PYMODULES='numpy')
        numpytool.Append( CPPPATH=join( numpy.__path__[0], 'core/include/' ))
        numpytool.Exists(True)

    context.Message("checking for scipy ... ")
    try:
        import scipy
    except ImportError:
        context.Result(context_failure)
        need_pkg( context.env, 'scipy', fatal=False)
    else:
        context.Result(context_success)
        context.env.Append(PYMODULES='scipy')
        pytool.Append(PYMODULES='scipy')

    context.Message("checking for pyct ... ")
    try:
        import pyct
    except ImportError:
        context.Result(context_failure)
    else:
        context.Result(context_success)
        context.env.Append(PYMODULES='pyct')
        pytool.Append(PYMODULES='pyct')
        
    pytool.Exists(True)
        

def intel(context):
    '''Trying to fix weird intel setup.'''
    libdirs = string.split(os.environ.get('LD_LIBRARY_PATH',''),':')
    libs = filter (lambda x: re.search('intel',x) and os.path.isdir(x),
                   libdirs)
    context.env.Append(ENV={'LD_LIBRARY_PATH':string.join(libs,':')})
    for key in ('INTEL_FLEXLM_LICENSE','INTEL_LICENSE_FILE','IA32ROOT'):
        license = os.environ.get(key)
        if license:
            context.env.Append(ENV={key:license})

def options(opts):
    opts.Add('ENV','SCons environment')
    opts.Add('AR','Static library archiver')
    opts.Add('JPEG','The libjpeg library')
    opts.Add('OPENGL','OpenGL libraries','GL GLU glut')
    opts.Add('OPENGLFLAGS','Flags for linking OpenGL libraries')
    opts.Add('GLEW','GLEW library','GLEW')
    opts.Add('MPICC','MPI C compiler')
    opts.Add('OMP','OpenMP support')
    opts.Add('BLAS','The BLAS library')
    opts.Add('PPM','The netpbm library')
    opts.Add('CC','The C compiler')
    opts.Add('CCFLAGS','General options that are passed to the C compiler',
             ['-O2'])
    opts.Add('CPPPATH',
             'The list of directories that the C preprocessor will search',[])
    opts.Add('LIBPATH',
             'The list of directories that will be searched for libraries')
    opts.Add('LIBS',
             'The list of libraries that will be linked with executables')
    opts.Add('LINKFLAGS','General options that are passed to the linker')
    opts.Add('XLIBPATH','Location of X11 libraries')
    opts.Add('XLIBS','X11 libraries')
    opts.Add('XINC','Location of X11 headers')
    opts.Add('PROGPREFIX','The prefix used for executable file names','sf')
    opts.Add('API','Support for additional languages. Possible values: c++, fortran or f77, fortran-90 or f90, matlab, octave, python')
    opts.Add('CXX','The C++ compiler')
    opts.Add('CXXFLAGS','General options that are passed to the C++ compiler',
             ['-O2'])
    opts.Add('F77','The Fortran-77 compiler')
    opts.Add('F77FLAGS','General options that are passed to the F77 compiler',
             ['-O2'])
    opts.Add('CFORTRAN','Type of the Fortran-77 compiler (for cfortran.h)')
    opts.Add('F90','The Fortran-90 compiler')
    opts.Add('F90FLAGS','General options that are passed to the F90 compiler',
             ['-O2'])
    opts.Add('CFORTRAN90','Type of the Fortran-90 compiler (for cfortran.h)')
    opts.Add('F90MODSUFFIX','Suffix of Fortran-90 module interface files')
    opts.Add('MEXSUFFIX','Suffix for mex files')
    opts.Add('MEX','Matlab function compiler')
    opts.Add('MATLAB','Matlab interpreter')
    opts.Add('OCTAVE','Octave interpreter')
    opts.Add('MKOCTFILE','Octave function compiler')
    opts.Add('PYMODULES','List of Python modules available')

