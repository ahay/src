# Copyright (C) 2005-2010 University of Texas at Austin #
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys, os, glob, string, re, commands, types

try: # The subprocess module was introduced in Python 2.4
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False

import SCons

# The following adds all SCons SConscript API to the globals of this module.
version = map(int,string.split(SCons.__version__,'.')[:3])
if version[0] >= 1 or version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
    from SCons.Script import *
else:  # old style
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

# CONSTANTS -- DO NOT CHANGE
context_success = 1
context_failure = 0
unix_failure = 1

def escape_seq(keyword):
    'A few ANSI escape sequences for highlighting warnings and errors'
    # Bold for warnings, and yellow on red background for fatal errors
    # These are independent from user-selected terminal emulator background
    # Details at http://en.wikipedia.org/wiki/ANSI_escape_code
    ansi_codes = {
        'bold'  :1,
        'redbg' :41,
        'yellow':93,
        'end'   :0
        }
    return '\033[%dm' % ansi_codes[keyword]

def stderr_write(message, highlight_mode=None):
    'Make sure error messages stand out visually'

    if highlight_mode != None:
        if highlight_mode == 'bold':
            prependix = escape_seq('bold')
        elif highlight_mode == 'yellow_on_red':
            prependix = escape_seq('redbg') + escape_seq('yellow')
        message =  prependix + message + escape_seq('end')

    sys.stderr.write('\n  %s\n' % message)

plat = {'OS': 'unknown',
        'distro': 'unknown',
        'version': 'unknown',
        'arch': 'unknown',
        'cpu': 'unknown'}
pkg = {}

def need_pkg(pkgtype,fatal=True):
    global pkg, plat
    pkgnm = 'unknown'
    if pkg.has_key(pkgtype):
        if pkg[pkgtype].has_key(plat['distro']):
            pkgnm = pkg[pkgtype].get(plat['distro'])
            if fatal:
                stderr_write('Needed package: ' + pkgnm,'yellow_on_red')
            else:
                stderr_write('Optional package: ' + pkgnm,'bold')
    if fatal:
        stderr_write('Fatal missing dependency','yellow_on_red')
        sys.exit(unix_failure)

def check_all(context):

    # FDNSI = Failure Does Not Stop Installation
    identify_platform(context)
    cc  (context)
    ar  (context)
    libs(context)
    c99 (context) # FDNSI
    x11 (context) # FDNSI
    opengl(context) # FDNSI
    sfpen(context) # FDNSI
    ppm (context) # FDNSI
    tiff (context) # FDNSI
    gd  (context) # FDNSI
    plplot  (context) # FDNSI
    ffmpeg  (context) # FDNSI
    cairo(context) # FDNSI
    jpeg(context) # FDNSI
    blas(context) # FDNSI
    lapack(context) # FDNSI
    mpi (context) # FDNSI
    pthreads (context) # FDNSI
    omp (context) # FDNSI
    petsc(context) # FDNSI
    cuda(context) # FDNSI
    api = api_options(context)
    if 'c++' in api:
        cxx(context)
    if 'f77' in api:
        f77(context)
    if 'f90' in api:
        f90(context)
    if 'matlab' in api:
        matlab(context)
    if 'octave' in api:
        octave(context)
    if 'python' in api:
        python(context)
    if 'java' in api:
        java(context)

def identify_platform(context):
    global plat
    context.Message("checking platform ... ")
    plat['OS'] = context.env.get('PLATFORM',sys.platform)
    try:
        from platform import architecture, uname
        if sys.version_info[:2] < (2, 6): # Python < 2.6
            from platform import dist
        else:
            from platform import linux_distribution as dist

        plat['arch'] = architecture()[0]
        del architecture

        name = uname()[2].split('.')[-2]
        if plat['OS'] in ('linux', 'posix', 'linux2'):            
            if dist()[0].lower() == 'fedora':
                plat['OS'] = 'linux'
                plat['distro'] = 'fedora'
                plat['version'] = dist()[1]
            elif dist()[0].lower() == 'redhat':
                plat['OS'] = 'linux'
                plat['distro'] = 'rhel' # Red Hat Enterprise Linux
                plat['version'] = dist()[1]
            elif dist()[0].lower() == 'ubuntu':
                plat['OS'] = 'linux'
                plat['distro'] = 'ubuntu'
                plat['version'] = dist()[1]
            elif dist()[0].lower() == 'debian':
                plat['OS'] = 'linux'
                plat['distro'] = 'debian'
                plat['version'] = dist()[1]   
            elif dist()[0] == 'SuSE':
                plat['OS'] = 'linux'
                plat['distro'] = 'suse'        
            elif name[-7:] == 'generic':
                plat['OS'] = 'linux'
                plat['distro'] = 'ubuntu'
        elif plat['OS'] == 'sunos':
            if name[:2] == '10':
                plat['distro'] = '10' # Solaris 10
        elif plat['OS'] == 'darwin':
             plat['distro'] = uname()[2]
             plat['cpu'] = uname()[5] # i386 / powerpc
        elif plat['OS'] == 'irix':
             plat['distro'] = uname()[2]
        elif plat['OS'] in ('hp-ux', 'hpux'):
             plat['OS'] = 'hpux'
             plat['distro'] = uname()[2].split('.')[-2]
        del uname, dist
    except: # "platform" not installed. Python < 2.3
        # For each OS with Python < 2.3, should use specific
        # commands through os.system to find distro/version
        # Not known if what follows works everywhere:
        plat_nm = os.uname()[4]
        if plat_nm == 'x86_64':
            plat['arch'] = '64bit'
        elif plat_nm == 'i686':
            plat['arch'] = '32bit'
    context.Result('%(OS)s [%(distro)s]' % plat)

pkg['gcc'] = {'fedora':'gcc'}
pkg['libc'] = {'fedora':'glibc',
               'ubuntu':'libc6-dev'}

# A C compiler is needed by most Madagascar programs
# Failing this test stops the installation.
def cc(context):
    context.Message("checking for C compiler ... ")
    CC = context.env.get('CC',WhereIs('gcc'))
    if CC:
        context.Result(CC)
    else:
        context.Result(context_failure)
        need_pkg('gcc')
    text = '''
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    if string.rfind(CC,'icc') >= 0:
        intel(context)

    context.Message("checking if %s works ... " % CC)
    res = context.TryLink(text,'.c')
    context.Result(res)
    if not res:
        need_pkg('libc')
    if string.rfind(CC,'gcc') >= 0 and \
           string.rfind(CC,'pgcc') < 0:
        oldflag = context.env.get('CCFLAGS')
        for flag in ('-std=gnu99 -Wall -pedantic',
                     '-std=gnu99 -Wall -pedantic',
                     '-std=gnu9x -Wall -pedantic',
                     '-Wall -pedantic'):
            context.Message("checking if gcc accepts '%s' ... " % flag)
            context.env['CCFLAGS'] = oldflag + ' ' + flag
            res = context.TryCompile(text,'.c')
            context.Result(res)
            if res:
                break        
        if not res:
            context.env['CCFLAGS'] = oldflag
        # large file support
        (status,lfs) = commands.getstatusoutput('getconf LFS_CFLAGS')
        if not status and lfs:
            oldflag = context.env.get('CCFLAGS')
            context.Message("checking if gcc accepts '%s' ... " % lfs)
            context.env['CCFLAGS'] = oldflag + ' ' + lfs
            res = context.TryCompile(text,'.c')
            context.Result(res)
            if not res:
                context.env['CCFLAGS'] = oldflag

        # Mac OS X include path, library path, and link flags
        if plat['OS'] == 'darwin':
            context.env['LINKFLAGS'] = context.env.get('LINKFLAGS','') + \
                                       ' -framework Accelerate ' 
            if os.path.isdir('/opt'):   # paths for MacPorts
                context.env['CPPPATH'] = context.env.get('CPPPATH',[]) + \
                                         ['/opt/local/include',]
                context.env['LIBPATH'] = context.env.get('LIBPATH',[]) + \
                                         ['/opt/local/lib',]
            if os.path.isdir('/sw'):    # paths for Fink
                context.env['CPPPATH'] = context.env.get('CPPPATH',[]) + \
                                         ['/sw/include',]
                context.env['LIBPATH'] = context.env.get('LIBPATH',[]) + \
                                         ['/sw/lib',]
  

    elif plat['OS'] == 'sunos':
        context.env['CCFLAGS'] = string.replace(context.env.get('CCFLAGS',''),
                                                '-O2','-xO2')

pkg['ar']={'fedora':'binutils'}

# Used for building libraries.
def ar(context):
    context.Message("checking for ar ... ")
    AR = context.env.get('AR',WhereIs('ar'))
    if AR:
        context.Result(AR)
        context.env['AR'] = AR
    else:
        context.Result(context_failure)
        need_pkg('ar')

pkg['libs'] = {'fedora':'glibc-headers',
               'cygwin':'sunrpc (Setup...Net)'}

# Failing this check stops the installation.
def libs(context):
    context.Message("checking for libraries ... ")
    LIBS = context.env.get('LIBS','m')
    DYNLIB = context.env.get('DYNLIB')
    if not DYNLIB or DYNLIB[0].lower() == 'n':
        context.env['DYNLIB'] = ''
    else:
        context.env['DYNLIB'] = 'd'
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    if plat['OS'] in ('sunos', 'hpux'):
        LIBS.append('nsl')
    elif plat['OS'] == 'cygwin':
        LIBS.append('tirpc')
    elif plat['OS'] == 'darwin':
        LIBS.append('mx')
    elif plat['OS'] == 'interix':
        LIBS.append('rpclib')
    text = '''
    #include <rpc/types.h>
    #include <rpc/xdr.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(str(LIBS))
        context.env['LIBS'] = LIBS
    else:
        context.Result(context_failure)
        need_pkg('libs')

pkg['c99'] = {'fedora':'glibc-headers'}

# Complex number support according to ISO C99 standard
def c99(context):
    context.Message("checking complex support ... ")
    text = '''
    #include <complex.h>
    #include <math.h>
    int main(int argc,char* argv[]) {
    float complex c;
    float f;
    c = 0.0;
    f = cabsf(ccosf(c));
    return (int) f;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
    else:
        context.env['CPPDEFINES'] = \
            context.env.get('CPPDEFINES',[])+['NO_COMPLEX']
        context.Result(context_failure)
        need_pkg('c99', fatal=False)

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

pkg['xaw']={'fedora':'libXaw-devel',
            'ubuntu':'libxaw7-dev'}

# If this check is failed
# you may not be able to display .vpl images on the screen
def x11(context):
    text = '''
    #include <X11/Intrinsic.h>
    #include <X11/Xaw/Label.h>
    int main(int argc,char* argv[]) {
    return 0;
    }\n'''

    context.Message("checking for X11 headers ... ")
    INC = context.env.get('XINC','')
    if type(INC) is not types.ListType:
        INC = string.split(INC)

    oldpath = context.env.get('CPPPATH',[])

    res = None
    for path in filter(lambda x:
                       os.path.isfile(os.path.join(x,'X11/Xaw/Label.h')),
                       INC+xinc):
        context.env['CPPPATH'] = oldpath + [path,] 
        res = context.TryCompile(text,'.c')

        if res:
            context.Result(path)
            context.env['XINC'] = [path,]
            break

    if not res:
        context.Result(context_failure)
        need_pkg('xaw', fatal=False)
        context.env['XINC'] = None
        return

    context.Message("checking for X11 libraries ... ")
    LIB = context.env.get('XLIBPATH','')
    if type(LIB) is not types.ListType:
        LIB = string.split(LIB)

    oldlibpath = context.env.get('LIBPATH',[])
    oldlibs = context.env.get('LIBS',[])

    XLIBS = context.env.get('XLIBS')
    if XLIBS:
        if type(XLIBS) is not types.ListType:
            XLIBS = string.split(XLIBS,',')
    else:
        if  plat['OS'] == 'interix':
            XLIBS =  ['Xaw','Xt','Xmu','X11','Xext','SM','ICE']
        elif (plat['OS'] == 'linux' or plat['OS'] == 'posix') and \
                plat['distro'] != 'fedora':
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
        context.Result(context_failure)
        context.env['XLIBPATH'] = None

    context.env['CPPPATH'] = oldpath
    context.env['LIBPATH'] = oldlibpath
    context.env['LIBS'] = oldlibs

def check_pen(env,pen):
    if pen == 'xtpen' and (env.get('XINC') and env.get('XLIBPATH')):
        return 1
    if pen == 'oglpen' and (env.get('OPENGL') or env.get('OPENGLFLAGS')):
        return 1
    return 0

def sfpen(context):
    context.Message("checking for sfpen ... ")
    sfpen = context.env.get('SFPEN')
    
    if not sfpen:
        if plat['OS'] == 'cygwin' or plat['OS'] == 'darwin':
            pens = ('oglpen','xtpen')
        else:
            pens = ('xtpen','oglpen')
        for pen in pens:
            if check_pen(context.env,pen):
                sfpen = pen
                break
    else:
        if not check_pen(context.env,sfpen):
            sfpen = None

    if not sfpen:
        context.Result(context_failure)
        stderr_write('sfpen (for displaying .vpl images) will not be built.',
                     'bold')
    else:
        context.Result(sfpen)
        context.env['SFPEN'] = sfpen

pkg['netpbm'] = {'cygwin':'libnetpbm-devel (Setup...Devel)',
                 'darwin':'netpbm (fink)',
                 'fedora':'netpbm-devel',
                 'suse'  :'libnetpbm-devel',
                 'ubuntu':'libnetpbm10-dev'}

def ppm(context):
    context.Message("checking for ppm ... ")

    oldpath = context.env.get('CPPPATH',[])
    ppmpath = context.env.get('PPMPATH','/usr/include/netpbm')
    if os.path.isfile(os.path.join(ppmpath,'ppm.h')):
        context.env['CPPPATH'] = oldpath + [ppmpath]
    else:
        ppmpath = None

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <ppm.h>
    int main(int argc,char* argv[]) {
    ppm_init(&argc,argv);
    return 0;
    }\n'''
    for ppm in [context.env.get('PPM','netpbm'),'netpbm.10']:
        LIBS.append(ppm)
        res = context.TryLink(text,'.c')

        if res:
            context.Result(res)
            context.env['PPM'] = ppm
            context.env['PPMPATH'] = ppmpath    
            break
        else:
            LIBS.pop()

    if res:
        LIBS.pop()
    else:
        context.Result(context_failure)
        need_pkg('netpbm', fatal=False)
        context.env['PPM'] = None
    context.env['CPPPATH'] = oldpath

pkg['libtiff'] = {'suse':'libtiff-devel', 'fedora':'libtiff-devel'}

def tiff(context):
    context.Message("checking for tiff ... ")

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <tiffio.h>
    int main(int argc,char* argv[]) {
    TIFF *tiffout;
    tiffout = TIFFOpen("test.tif","wb");
    return 0;
    }\n'''
    tiff = context.env.get('TIFF','tiff')
    LIBS.append(tiff)
    res = context.TryLink(text,'.c')

    if res:
        context.Result(res)
        context.env['TIFF'] = tiff   
    else:
        context.Result(context_failure)
        context.env['TIFF'] = None
        stderr_write('sfbyte2tif, sftif2byte, and tiffpen will not be built.',
                     'bold')
        need_pkg('libtiff', fatal=False)

    LIBS.pop()

pkg['libgd'] = {'suse':'gd-devel',
                'ubuntu':'libgd2-xpm-dev'}

def gd(context):
    context.Message("checking for GD (PNG) ... ")

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <gd.h>
    int main(int argc,char* argv[]) {
    gdImagePtr image;
    image = gdImageCreate(10,10);
    gdImageSetClip(image,0,0,1,1);
    return 0;
    }\n'''
    gd = context.env.get('GD','gd')
    LIBS.append(gd)
    res = context.TryLink(text,'.c')

    if res:
        context.Result(res)
        context.env['GD'] = gd

        context.Message("checking for GD (GIF) ... ")
        text = '''
        #include <stdio.h>
        #include <gd.h>
        int main(int argc,char* argv[]) {
        gdImagePtr image;
        image = gdImageCreate(10,10);
        gdImageGifAnimBegin(image, stdout, 1, 0);
        return 0;
        }\n'''        
        res2 = context.TryLink(text,'.c')

        if res2:
            context.Result(res2)
            context.env['GIFANIM'] = True
        else:
            context.Result(context_failure)
            context.env['GIFANIM'] = None
    else:
        context.Result(context_failure)
        need_pkg('libgd', fatal=False)
        context.env['GD'] = None
    LIBS.pop()

pkg['plplot'] = {'fedora':'plplot-devel',
                 'suse':'libplplot-devel',
                 'ubuntu':'libplplot-dev'}

def plplot(context):
    context.Message("checking for plplot ... ")

    oldpath = context.env.get('CPPPATH',[])
    plplotpath = context.env.get('PLPLOTPATH')
    oldlibpath = context.env.get('LIBPATH',[])
    plplotlibpath = context.env.get('PLPLOTLIBPATH')
    if plplotpath and os.path.isfile(os.path.join(plplotpath,'plplot.h')):
        context.env['CPPPATH'] = oldpath + [plplotpath]
    else:
        for top in ('/usr/include','/usr/local/include',
                    '/sw/include','/opt/local/include'):
            plplotpath = os.path.join(top,'plplot')
            if os.path.isfile(os.path.join(plplotpath,'plplot.h')):
                context.env['CPPPATH'] = oldpath + [plplotpath]
                break
    if plplotlibpath:
        context.env['LIBPATH'] = oldlibpath + [plplotlibpath]

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <plplot.h>
    #include <plplotP.h>
    #include <plstrm.h>
    int main (int argc, char *argv[]) {
    plinit ();
    pladv (0);
    return 0;
    }\n'''
    plplot = context.env.get('PLPLOT','plplotd')
    LIBS.append(plplot)
    LIBS.append('ltdl')
    res = context.TryLink(text,'.c')

    if res:
        context.Result(res)
        context.env['PLPLOT'] = plplot
        context.env['PLPLOTPATH'] = plplotpath
        context.env['PLPLOTLIBPATH'] = plplotlibpath
    else:
        context.Result(context_failure)
        need_pkg('plplot', fatal=False)
        context.env['PLPLOT'] = None
    context.env['CPPPATH'] = oldpath
    context.env['LIBPATH'] = oldlibpath    
    LIBS.pop()
    LIBS.pop()

pkg['ffmpeg'] = {'fedora':'ffmpeg-devel',
                 'suse':'ffmpeg-devel',
                 'ubuntu':'libavcodec-dev'}

def ffmpeg(context):
    context.Message("checking for ffmpeg ... ")

    oldpath = context.env.get('CPPPATH',[])
    ffmpegpath = context.env.get('FFMPEGPATH')
    if ffmpegpath and os.path.isfile(os.path.join(ffmpegpath,'avcodec.h')):
        context.env['CPPPATH'] = oldpath + [ffmpegpath]
    else:
        for top in ('/usr/include','/usr/local/include',
                    '/sw/include','/opt/local/include'):
            ffmpegpath = os.path.join(top,'ffmpeg')
            if os.path.isfile(os.path.join(ffmpegpath,'avcodec.h')):
                context.env['CPPPATH'] = oldpath + [ffmpegpath]
                break
            ffmpegshortpath = ffmpegpath
            ffmpegpath = os.path.join(ffmpegshortpath,'libavcodec')
            if os.path.isfile(os.path.join(ffmpegpath,'avcodec.h')):
                context.env['CPPPATH'] = oldpath + [ffmpegshortpath,ffmpegpath]
                ffmpegpath = [ffmpegshortpath,ffmpegpath]
                break
            ffmpegpath = os.path.join(top,'libavcodec')
            if os.path.isfile(os.path.join(ffmpegpath,'avcodec.h')):
                context.env['CPPPATH'] = oldpath + [ffmpegpath]
                break

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <avcodec.h>
    int main (int argc, char *argv[]) {
    avcodec_init ();
    avcodec_register_all ();
    return 0;
    }\n'''
    ffmpeg = context.env.get('FFMPEG','avcodec')
    LIBS.append(ffmpeg)
    res = context.TryLink(text,'.c')

    if res:
        context.Result(res)
        context.env['FFMPEG'] = ffmpeg
        context.env['FFMPEGPATH'] = ffmpegpath
    else:
        context.Result(context_failure)
        need_pkg('ffmpeg', fatal=False)
        context.env['FFMPEG'] = None
    context.env['CPPPATH'] = oldpath    
    LIBS.pop()


pkg['cairo'] = {'suse':'cairo-devel',
                'ubuntu':'libcairo2-dev'}

def cairo(context):
    context.Message("checking for cairo (PNG) ... ")

    oldpath = context.env.get('CPPPATH',[])

    cairopath = context.env.get('CAIROPATH')
    if cairopath and os.path.isfile(os.path.join(cairopath,'cairo.h')):
        context.env['CPPPATH'] = oldpath + [cairopath]
    else:
        for top in ('/usr/include','/usr/local/include','/sw/include'):
            cairopath = os.path.join(top,'cairo')
            if os.path.isfile(os.path.join(cairopath,'cairo.h')):
                context.env['CPPPATH'] = oldpath + [cairopath]
                break

    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    text = '''
    #include <cairo.h>
    int main(int argc,char* argv[]) {
    cairo_surface_t *surface;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 10, 10);
    return 0;
    }\n'''
    cairo = context.env.get('CAIRO','cairo')
    LIBS.append(cairo)
    res = context.TryLink(text,'.c')

    if res:
        context.Result(res)
        context.env['CAIROPNG'] = cairo
        context.env['CAIROPATH'] = cairopath   
    else:
        context.Result(context_failure)
        need_pkg('cairo', fatal=False)
        context.env['CAIROPNG'] = None

    if res:
        for format in ('svg','pdf'):
            cap = format.upper()
            context.Message("checking for cairo (%s) ... " % cap)
            text = '''
            #include <cairo.h>
            #include <cairo-%s.h>
            int main(int argc,char* argv[]) {
            cairo_surface_t *surface;
            surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 10, 10);
            return 0;
            }\n''' % format
            res2 = context.TryLink(text,'.c')

            if res2:
                context.Result(res2)
                context.env['CAIRO' + cap] = cairo
            else:
                context.Result(context_failure)
                context.env['CAIRO' + cap] = None

    context.env['CPPPATH'] = oldpath    
    LIBS.pop()

pkg['jpeg'] = {'fedora':'libjpeg-devel',
               'ubuntu':'libjpeg62-dev'}

# If this test is failed, no writing to jpeg files
def jpeg(context):
    context.Message("checking for jpeg ... ")
    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    jpeg = context.env.get('JPEG','jpeg')
    LIBS.append(jpeg)
    text = '''
    #include <stdio.h>
    #include <jpeglib.h>
    int main(int argc,char* argv[]) {
    struct jpeg_compress_struct jpeg;
    struct jpeg_error_mgr jpeg_err;
    jpeg.err = jpeg_std_error(&jpeg_err);
    return 0;
    }\n'''

    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['JPEG'] = jpeg
    else:
        context.Result(context_failure)
        stderr_write('sfbyte2jpg, sfjpg2byte, and jpegpen will not be built.',
                     'bold')
        need_pkg('jpeg', fatal=False)
        context.env['JPEG'] = None

    LIBS.pop()

pkg['opengl'] = {'fedora':'mesa-libGL-devel + freeglut + freeglut-devel',
                 'suse'  :'freeglut-devel',
                 'ubuntu':'freeglut3-dev',
                 'cygwin':'opengl (Setup...Graphics)'}

# If this test is failed, no opengl programs
def opengl(context):
    global plat
    context.Message("checking for OpenGL ... ")
    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    LINKFLAGS = context.env.get('LINKFLAGS','')
    CPPPATH = context.env.get('CPPPATH',[])
    oglflags = None
    oglpath = None

    if plat['OS'] == 'darwin':
        ogl = []
        oglflags = ' -framework AGL -framework OpenGL -framework GLUT'
        context.env['LINKFLAGS'] = LINKFLAGS + oglflags
    elif plat['OS'] == 'cygwin':
        ogl = context.env.get('OPENGL','opengl32 glu32 glut32')
        oglpath = '/usr/include/opengl'
        context.env['CPPPATH'] = CPPPATH + [oglpath]
    else:
        ogl = context.env.get('OPENGL','GL GLU glut')

    if type(ogl) is not types.ListType:
        ogl = string.split(ogl)
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
        context.env['OPENGLPATH'] = oglpath
    else:
        context.env['OPENGL'] = None 
        context.env['OPENGLFLAGS'] = None
        context.env['OPENGLPATH'] = None
        context.Result(context_failure)
        need_pkg('opengl', fatal=False)

    context.env['LIBS'] = LIBS
    context.env['LINKFLAGS'] = LINKFLAGS
    context.env['CPPPATH'] = CPPPATH

pkg['blas'] = {'fedora':'blas + blas-devel + atlas + atlas-devel',
               'ubuntu':'libblas-dev'}

def blas(context):
    context.Message("checking for BLAS ... ")
    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
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
        if plat['OS'] == 'cygwin':
            context.env['ENV']['PATH'] = context.env['ENV']['PATH'] + \
                                         ':/lib/lapack'
    else:
        # some systems require cblas and atlas
        LIBS.append('cblas')
        LIBS.append('atlas')
        res = context.TryLink(text,'.c')
        if res:
            context.Result(res)
            context.env['LIBS'] = LIBS
            context.env['BLAS'] = blas
        else:
            context.Result(context_failure)
            context.env['CPPDEFINES'] = \
                context.env.get('CPPDEFINES',[]) + ['NO_BLAS']
            LIBS.pop()
            LIBS.pop()
            LIBS.pop()
            context.env['BLAS'] = None
            need_pkg('blas', fatal=False)

def lapack(context):
    context.Message("checking for LAPACK ... ")
    LIBS = context.env.get('LIBS','m')
    if type(LIBS) is not types.ListType:
        LIBS = string.split(LIBS)
    blas = context.env.get('BLAS')
    if not blas:
        blas = 'blas'
    lapack = context.env.get('LAPACK','lapack')
    mylibs = [lapack,blas]
    LIBS.extend(mylibs)
    text = '''
    int main(int argc,char* argv[]) {
    int i=1;
    return 0;
    }\n'''
    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['LAPACK'] = mylibs
        if plat['OS'] == 'cygwin':
            context.env['ENV']['PATH'] = context.env['ENV']['PATH'] + \
                                         ':/lib/lapack'
    else:
        # some systems require cblas and atlas
        mylibs.extend(['cblas','atlas'])
        LIBS.extend(['cblas','atlas'])
        res = context.TryLink(text,'.c')
        if res:
            context.Result(res)
            context.env['LAPACK'] = mylibs
        else:
            context.Result(context_failure)
            context.env['LAPACK'] = None
            need_pkg('lapack', fatal=False)
        LIBS.pop()
        LIBS.pop()
    LIBS.pop()
    LIBS.pop()

pkg['mpi'] = {'fedora':'openmpi, openmpi-devel, openmpi-libs'}

def mpi(context):
    context.Message("checking for MPI ... ")
    path = os.environ['PATH']
    if plat['OS'] == 'linux':
        if plat['distro'] == 'fedora':
            path += ':/usr/lib64/openmpi/bin/'
    mpicc = context.env.get('MPICC',WhereIs('mpicc', path))
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
        context.env.Append(ENV={'MPICH_CC':cc,'MPICH_CLINKER':cc})
        res = context.TryLink(text,'.c')
        context.env['CC'] = cc
    else: # mpicc not found
        context.Result(context_failure)
        res = None
    if res:
        context.Result(res)
        context.env['MPICC'] = mpicc
    else:
        context.Result(context_failure)
        need_pkg('mpi', fatal=False)
        context.env['MPICC'] = None

def cuda(context):
    context.Message("checking for CUDA ... ")

    CUDA_TOOLKIT_PATH = context.env.get('CUDA_TOOLKIT_PATH',
                                        os.environ.get('CUDA_TOOLKIT_PATH'))
    if CUDA_TOOLKIT_PATH:
        context.Result(CUDA_TOOLKIT_PATH)
        if os.path.isdir(CUDA_TOOLKIT_PATH):
            context.env['CUDA_TOOLKIT_PATH'] = CUDA_TOOLKIT_PATH
        else:
            stderr_write('Set CUDA_TOOLKIT_PATH to the location of CUDA toolkit',
                         'bold')
            context.env['CUDA_TOOLKIT_PATH'] = None
    else:
        home = os.environ.get('HOME', '')
        paths = [
            os.path.join(home, 'NVIDIA_CUDA_TOOLKIT'),
            os.path.join(home, 'Apps', 'NVIDIA_CUDA_TOOLKIT'),
            os.path.join(home, 'Apps', 'CudaToolkit'),
            os.path.join(home, 'Apps', 'CudaTK'),
            os.path.join('/usr', 'local', 'NVIDIA_CUDA_TOOLKIT'),
            os.path.join('/usr', 'local', 'CUDA_TOOLKIT'),
            os.path.join('/usr', 'local', 'cuda_toolkit'),
            os.path.join('/usr', 'local', 'CUDA'),
            os.path.join('/usr', 'local', 'cuda'),
            os.path.join('/Developer', 'NVIDIA CUDA TOOLKIT'),
            os.path.join('/Developer', 'CUDA TOOLKIT'),
            os.path.join('/Developer', 'CUDA')
            ]
        for path in paths:
            if os.path.isdir(path):
                CUDA_TOOLKIT_PATH = path
                break
        if CUDA_TOOLKIT_PATH:
            context.Result(CUDA_TOOLKIT_PATH)
            context.env['CUDA_TOOLKIT_PATH'] = CUDA_TOOLKIT_PATH
        else:
            context.Result(context_failure)
            context.env['CUDA_TOOLKIT_PATH'] = None
            return

    path = ':'.join([os.environ['PATH'],os.path.join(CUDA_TOOLKIT_PATH,'bin')])
    nvcc = context.env.get('NVCC',WhereIs('nvcc',path))
    if nvcc:
        context.Message("checking if %s works ... " % nvcc)
        # Try compiling with nvcc instead of cc
        text = '''
        #include <cuda.h>
        #include <cuda_runtime_api.h>
        int main(int argc,char* argv[]) {
        cuInit (0); /* Use first GPU device */
        cudaSetDevice (0);
        }\n'''
        cc = context.env.get('CC')
        cflags = context.env.get('CCFLAGS')
        libs = context.env.get('LIBS')
        libpath = context.env.get('LIBPATH')
        linkflags = context.env.get('LINKFLAGS')
        context.env['CC'] = nvcc
        context.env['CCFLAGS'] = '--x=cu'
        context.env['LIBS'] = ['cuda','cudart']
        context.env['LIBPATH'] = [os.path.join(CUDA_TOOLKIT_PATH,'lib64'),
                                  os.path.join(CUDA_TOOLKIT_PATH,'lib')]
        context.env['LINKFLAGS'] = ''
        res = context.TryLink(text,'.c')
        context.env['CC'] = cc
        context.env['CCFLAGS'] = cflags
        context.env['LIBS'] = libs
        context.env['LIBPATH'] = libpath
        context.env['LINKFLAGS'] = linkflags
    else:
        context.Result(context_failure)
        res = None
    if res:
        context.Result(res)
        context.env['NVCC'] = nvcc
    else:
        context.Result(context_failure)
        context.env['NVCC'] = None

pkg['petsc'] = {'ubuntu':'petsc-dev',
                'fedora':'petsc-devel'}

def petsc(context):
    petscdir = context.env.get('PETSCDIR',os.environ.get('PETSC_DIR'))

    if not petscdir:
        return
    testdir = os.path.join(os.getcwd(),'user/petsc')
    
    # Run make in order to catch PETSc compilation options
    if have_subprocess: # use subprocess.Popen() if possible, for Py 2.4 and up
        popen = subprocess.Popen('make PETSC_DIR=%s options' % petscdir, 
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 cwd=testdir)
        if popen.wait() != 0:
            return
        makeout = popen.stdout.read()
    else: # otherwise use os.popen2(), deprecated in Py 2.6
        makeout = os.popen2('make PETSC_DIR=%s -C %s options' % 
                            (petscdir,testdir))[1].read()
        
    context.Message("checking for PETSc ... ")

    # Compiler name
    compiler = re.compile(r'^\S*')
    petsccc = compiler.findall(makeout)
    petsccc = WhereIs(petsccc[0])
    # All include paths (-I/...)
    includes = re.compile(r'\s\-I(\S*)')
    petscpath = includes.findall(makeout)
    # All lib paths (-L/...)
    libpaths = re.compile(r'\s\-L(\S*)')
    petsclibpath = libpaths.findall(makeout)
    # All libs (-l...)
    libs = re.compile(r'\s\-l(\S*)')
    petsclibs = libs.findall(makeout)

    oldcc = context.env.get('CC')
    oldpath = context.env.get('CPPPATH',[])
    oldlibpath = context.env.get('LIBPATH',[])
    oldlibs = context.env.get('LIBS',[])

    context.env['CC'] = petsccc
    context.env['CPPPATH'] = oldpath + [petscpath,] 
    context.env['LIBPATH'] = oldlibpath + [petsclibpath,]
    context.env['LIBS'] = oldlibs + [petsclibs,]

    text = '''
    #include <petscksp.h>
    int main(int argc,char* argv[]) {
    KSP Solver; Mat A;
    MPI_Comm comm = MPI_COMM_WORLD;
    PetscInitialize (&argc, &argv, 0, 0);
    MatCreate (comm, &A);
    KSPCreate (comm, &Solver);
    KSPSetType (Solver, KSPGMRES);
    PetscFinalize ();
    }\n'''
    res = context.TryLink(text,'.c')

    context.env['CPPPATH'] = oldpath
    context.env['LIBPATH'] = oldlibpath
    context.env['LIBS'] = oldlibs
    context.env['CC'] = oldcc

    if res:
        context.Result(res)
        context.env['PETSCDIR'] = petscdir
        context.env['PETSCPATH'] = petscpath
        context.env['PETSCLIBPATH'] = petsclibpath
        context.env['PETSCLIBS'] = petsclibs
        context.env['PETSCCC'] = petsccc
    else:
        context.Result(context_failure)
        need_pkg('petsc', fatal=False)

def ncpus():
    'Detects number of CPUs'
    if plat['OS'] in ('linux','posix'):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            nr_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if type(nr_cpus) is int:
                if nr_cpus > 0:
                    return nr_cpus
    elif plat['OS'] == 'darwin':
        if have_subprocess: # use subprocess.Popen() if possible, for Py 2.4 and up
            pipe = subprocess.Popen('sysctl -n hw.ncpu', shell=True, stdout=subprocess.PIPE).stdout
            nr_cpus = int(pipe.read())
        else: # otherwise use os.popen2(), deprecated in Py 2.6
            nr_cpus = int(os.popen2('sysctl -n hw.ncpu')[1].read())
        if nr_cpus > 0:
            return nr_cpus
    return 2 # default number of processors

pkg['omp'] = {'fedora':'libgomp'}

def omp(context):
    if ncpus() == 1:
        context.env['OMP'] = False
        return # only 1 cpu. OMP not needed
    context.Message("checking for OpenMP ... ")
    LIBS  = context.env.get('LIBS',[])
    CC    = context.env.get('CC','gcc')
    flags = context.env.get('CCFLAGS','')
    lflags = context.env.get('LINKFLAGS','')
    pgcc =  (string.rfind(CC,'pgcc') >= 0)
    gcc = (string.rfind(CC,'gcc') >= 0)
    icc = (string.rfind(CC,'icc') >= 0)
    if pgcc:
        CCFLAGS = flags + ' -mp'
        LINKFLAGS = lflags + ' -mp'
    elif gcc:
        LIBS.append('gomp')
        CCFLAGS = flags + ' -fopenmp'
        LINKFLAGS = lflags
    elif icc:
        CCFLAGS = flags + ' -openmp -D_OPENMP'
        LINKFLAGS = lflags + ' -openmp' 
    else:
        CCFLAGS = flags
        LINKFLAGS = lflags

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
    context.env['LINKFLAGS'] = LINKFLAGS
    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['OMP'] = True
    else:
        context.Result(context_failure)
        need_pkg('omp', fatal=False)
        if gcc:
            LIBS.pop()
        context.env['LIBS'] = LIBS
        context.env['CCFLAGS'] = flags
        context.env['LINKFLAGS'] = lflags
        context.env['OMP'] = False

def pthreads(context):
    context.Message("checking for Posix threads ... ")

    flags = context.env.get('LINKFLAGS','')
    LIBS  = context.env.get('LIBS',[])
    CC    = context.env.get('CC','gcc')
    pgcc =  (string.rfind(CC,'pgcc') >= 0)
    gcc = (string.rfind(CC,'gcc') >= 0)
    icc = (string.rfind(CC,'icc') >= 0)
    if icc or pgcc:
        LIBS.append('pthread')
    elif gcc and plat['OS'] != 'darwin' and plat['OS'] != 'cygwin':
        context.env.Append(LINKFLAGS='-pthread')

    text = '''
    #include <pthread.h>
    #include <stdlib.h>
    int main(void) {
    pthread_t peers;
    pthread_create(&peers, NULL, NULL, NULL);
    return 0;
    }
    '''
    context.env['LIBS'] = LIBS
    res = context.TryLink(text,'.c')
    if res:
        context.Result(res)
        context.env['PTHREADS'] = True
    else:
        context.Result(context_failure)
        need_pkg('pthreads', fatal=False)
        if icc or pgcc:
            LIBS.pop()
        context.env['LIBS'] = LIBS
        context.env.Replace(LINKFLAGS=flags)
        context.env['PTHREADS'] = False

def api_options(context):
    context.Message("checking API options ... ")
    api = context.env.get('API')
    if api:
        api = string.split(string.lower(api),',')
    else:
        api = []

    valid_api_options = ['','c++', 'fortran', 'f77', 'fortran-90',
                         'f90', 'python', 'matlab', 'octave', 'java']

    for option in api:
        if not option in valid_api_options:
            api.remove(option)

    # Make tests for fortrans in API easy
    for i in range(len(api)):
        if api[i] == 'fortran':
            api[i] = 'f77'
        elif api[i] == 'fortran-90':
            api[i] = 'f90'

    # Eliminate duplicates if user was redundant.
    # For Py 2.4 and up this can be done more elegantly with:
    # api = list(set(api))
    api_dict = {}
    for i in api: api_dict[i] = 0
    api = api_dict.keys()

    # Improve output readability
    if api == ['']:
        context.Result('none')
    else:
        context.Result(str(api))
    context.env['API'] = api

    return api

pkg['c++'] = {'fedora':'gcc-c++',
              'suse'  :'gcc-c++',
              'ubuntu':'g++'}

# For the C++ API
def cxx(context):
    context.Message("checking for C++ compiler ... ")
    CXX = context.env.get('CXX')
    if CXX:
        context.Result(CXX)
    else:
        context.Result(context_failure)
        need_pkg('c++')
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
        oldflag = context.env.get('CXXFLAGS')
        for flag in ['-Wall -pedantic']:
            context.Message("checking if g++ accepts '%s' ... " % flag)
            context.env['CXXFLAGS'] = oldflag + ' ' + flag
            res = context.TryCompile(text,'.cc')
            context.Result(res)
            if res:
                break
        if not res:
            context.env['CXXFLAGS'] = oldflag

# Used in checks for both f77 and f90
fortran = {'g77':'f2cFortran',
           'f77':'f2cFortran',
           'ifort':'f2cFortran',
           'gfortran':'NAGf90Fortran',
           'gfc':'NAGf90Fortran',
           'f2c':'f2cFortran'}

pkg['f77'] = {'fedora':'gcc-gfortran',
              'ubuntu':'g77'}

def f77(context):
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
        context.Result(F77)
        context.env['F77'] = F77
    else:
        context.Result(context_failure)
        need_pkg('f77')
    if os.path.basename(F77) == 'ifc' or os.path.basename(F77) == 'ifort':
        intel(context)
        context.env.Append(F77FLAGS=' -Vaxlib')
        context.env['FORTRAN'] = F77
    
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
        sys.exit(unix_failure)
    F77base = os.path.basename(F77)
    if F77base[:3] == 'f77' and plat['OS'] == 'sunos':
        cfortran = 'sunFortran'
    else:
        cfortran = fortran.get(F77base,'NAGf90Fortran')
    context.env['CFORTRAN'] = cfortran 
    context.Message("checking %s type ... " % F77)
    context.Result(cfortran)

pkg['f90'] = {'fedora':'gcc-gfortran',
              'suse'  :'gcc-fortran',
              'ubuntu':'gfortran'}

def f90_write_autofile(extension,content):
    filename=os.path.join('api','f90','ptr_sz.'+extension)
    handle = open(filename,'w')
    handle.write(content)
    handle.close()

def f90_write_ptr_sz():
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
        context.Result(F90)
    else:
        context.Result(context_failure)
        need_pkg('f90')
    if os.path.basename(F90) == 'ifc' or os.path.basename(F90) == 'ifort':
        intel(context)
        context.env.Append(F90FLAGS=' -Vaxlib')
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
        sys.exit(unix_failure)
    base = os.path.basename(F90)
    context.Message("checking %s type ... " % base)
    cfortran = fortran.get(base,'NAGf90Fortran')
    context.env['CFORTRAN90'] = cfortran
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
    context.Result(suffix)
    f90_write_ptr_sz()

def matlab(context):
    context.Message("checking for Matlab ... ")
    matlab = WhereIs('matlab')
    if matlab:
        context.Result(matlab)
        RSFROOT_lib = os.path.join(context.env.get('RSFROOT'),'lib')
        MATLABPATH = os.environ.get('MATLABPATH')
        if MATLABPATH:
            MATLABPATH += ':' + RSFROOT_lib
        else:
            MATLABPATH = RSFROOT_lib
        context.env['MATLAB'] = 'MATLABPATH=%s %s ' \
                                '-nosplash -nojvm -nodesktop' %(MATLABPATH,
                                                                matlab)
    else:
        context.Result(context_failure)
        stderr_write('Please install Matlab.','yellow_on_red')
        context.env['MATLAB'] = None
        sys.exit(unix_failure)

    context.Message("checking for mex ... ")
    mex = WhereIs('mex')
    if mex:
        context.Result(mex)
        context.env['MEX'] = mex
        if plat['OS'] == 'darwin':
            maci64 = os.environ.get('MACI64',0)                    
            context.env.Append(ENV={'MACI64':maci64})
    else:
        context.Result(context_failure)
        stderr_write('Please install mex.','yellow_on_red')
        context.env['MEX'] = None
        sys.exit(unix_failure)

    # See http://www.mathworks.com/access/helpdesk/help/techdoc/ref/mex.html
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

pkg['octave'] = {'fedora':'octave',
                 'ubuntu':'octave'}

pkg['mkoctave'] = {'fedora':'octave-devel',
                   'ubuntu':'octave-headers'}

def octave(context):
    context.Message("checking for Octave ... ")
    octave = WhereIs('octave')
    if octave:
        context.Result(octave)
        context.env['OCTAVE'] = octave
        context.Message("checking for mkoctfile ... ")
        mkoctfile = WhereIs('mkoctfile')
        if mkoctfile:
            context.Result(mkoctfile)
            context.env['MKOCTFILE'] = mkoctfile
        else:
            context.Result(context_failure)
            stderr_write('Please install mkoctfile.','bold')
            need_pkg('mkoctfile')
    else: # octave not found
        context.Result(context_failure)
        stderr_write('Please install Octave.','bold')
        need_pkg('octave')

pkg['swig'] = {'fedora':'swig',
               'suse'  :'swig',
               'ubuntu':'swig'}
pkg['numpy'] = {'fedora':'numpy',
                'ubuntu':'python-scipy, python-numpy-dev'}

def python(context):
    context.Message("checking for SWIG ... ")
    if 'swig' in Environment().get('TOOLS'):
        context.Result( WhereIs('swig') )
    else:
        context.Result(context_failure)
        need_pkg('swig')

    context.Message("checking for numpy ... ")
    try:
        import numpy
        context.Result(context_success)
    except:
        context.Result(context_failure)
        need_pkg('numpy')


pkg['java-devel'] = {'ubuntu':'openjdk-6-jdk'}
pkg['minesjtk'] = {}

def java(context):
    context.Message("checking for javac ... ")
    JAVAC = context.env.get('JAVAC',WhereIs('javac'))
    if JAVAC:
        context.Result(JAVAC)
        context.env['JAVAC'] = JAVAC
    else:
        context.Result(context_failure)
        need_pkg('java-devel')
    
    context.Message("checking for JAVA_HOME ... ")
    JAVA_HOME = context.env.get('JAVA_HOME',
                                os.environ.get('JAVA_HOME',
                                               '/usr/lib/jvm/java-6-openjdk'))
    if not JAVA_HOME:  # Check for JAVA_SDK as well, Mac fix
        JAVA_HOME = context.env.get('JAVA_SDK',os.environ.get('JAVA_SDK'))
    if JAVA_HOME:
        context.Result(JAVA_HOME)
        context.env['JAVA_HOME'] = JAVA_HOME
    else:
        context.Result(context_failure)
        need_pkg('java-devel')

    context.Message("Checking for Mines JTK ...")
    pkg['minesjtk'][plat['distro']] = 'Mines JTK http://inside.mines.edu/~dhale/jtk/\n\tSet MINESJTK to the location of edu_mines_jtk.jar'

    MINESJTK = context.env.get('MINESJTK',os.environ.get('MINESJTK'))
    if MINESJTK:
        context.Result(MINESJTK)
        if os.path.isdir(MINESJTK):
            MINESJTK = os.path.join(MINESJTK, 'edu_mines_jtk.jar')
        if os.path.isfile(MINESJTK) and \
                os.path.basename(MINESJTK) == 'edu_mines_jtk.jar':            
            context.env['MINESJTK'] = MINESJTK
        else:
            stderr_write('Set MINESJTK to the location of edu_mines_jtk.jar',
                         'bold')
            context.env['MINESJTK'] = None
    else:
        context.Result(context_failure)
        need_pkg('minesjtk', fatal=False)

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

def set_options(env,my_opts=None):
    'get options from config file'
    from rsf.prog import RSFROOT
    config = os.path.join(RSFROOT, 'etc', 'madagascar','config.py')
    if not os.path.isfile(config):
        return
    opts = options(config)
    if my_opts:
        for opt in my_opts.keys():
            opts.Add(opt,my_opts[opt])
    opts.Update(env)

def options(file):
    global version

    if version[0] < 1 or (version[0] == 1 and version[1] < 2):
        opts=Options(file)
    else:
        opts=Variables(file)

    opts.Add('ENV','SCons environment')
    opts.Add('RSFROOT','Top Madagascar installation directory')
    opts.Add('AR','Static library archiver')
    opts.Add('JPEG','The libjpeg library')
    opts.Add('OPENGL','OpenGL libraries')
    opts.Add('OPENGLFLAGS','Flags for linking OpenGL libraries')
    opts.Add('OPENGLPATH','Path to OpenGL headers')
    opts.Add('MPICC','MPI C compiler')
    opts.Add('PETSCDIR','Portable, Extensible Toolkit for Scientific computation - installation directory')
    opts.Add('PETSCPATH','PETSc - path to headers')
    opts.Add('PETSCLIBPATH','PETSc - path to libraries')
    opts.Add('PETSCLIBS','PETSc - libraries')
    opts.Add('PETSCCC','PETSc - compiler')
    opts.Add('OMP','OpenMP support')
    opts.Add('PTHREADS','Posix threads support')
    opts.Add('BLAS','The BLAS library')
    opts.Add('LAPACK','The LAPACK library')
    opts.Add('PPM','The netpbm library')
    opts.Add('TIFF','The libtiff library')
    opts.Add('GD','The GD library')
    opts.Add('GIFANIM','Support for GIF animation in the GD library')
    opts.Add('PLPLOT','PLPLOT library')
    opts.Add('PLPLOTPATH','Path to PLPLOT headers')
    opts.Add('PLPLOTLIBPATH','Path to PLPLOT library')
    opts.Add('FFMPEG','The ffmpeg library')
    opts.Add('FFMPEGPATH','Path to ffmpeg codec headers')
    opts.Add('CAIRO','The cairo library')
    opts.Add('CAIROPNG','The cairo library (PNG support)')
    opts.Add('CAIROSVG','The cairo library (SVG support)')
    opts.Add('CAIROPDF','The cairo library (PDF support)')
    opts.Add('CAIROPATH','Path to cairo header files')
    opts.Add('PPMPATH','Path to netpbm header files')
    opts.Add('SFPEN','Preference for sfpen')
    opts.Add('DYNLIB','Compiling with dynamic libraries')
    opts.Add('CC','The C compiler')
    opts.Add('CCFLAGS','General options that are passed to the C compiler',
             '-O2')
    opts.Add('CPPPATH',
             'The list of directories that the C preprocessor will search')
    opts.Add('CPPDEFINES',
             'List of defines for the C preprocessor')
    opts.Add('LIBPATH',
             'The list of directories that will be searched for libraries')
    opts.Add('LIBS',
             'The list of libraries that will be linked with executables')
    opts.Add('LINKFLAGS','General options that are passed to the linker')
    opts.Add('XLIBPATH','Location of X11 libraries')
    opts.Add('XLIBS','X11 libraries')
    opts.Add('XINC','Location of X11 headers')
    opts.Add('PROGPREFIX','The prefix used for executable file names','sf')
    opts.Add('API','Support for additional languages. Possible values: c++, fortran or f77, fortran-90 or f90, matlab, octave, python, java')
    opts.Add('CXX','The C++ compiler')
    opts.Add('CXXFLAGS','General options that are passed to the C++ compiler',
             '-O2')
    opts.Add('F77','The Fortran-77 compiler')
    opts.Add('FORTRAN','The generic Fortran compiler')
    opts.Add('F77FLAGS','General options that are passed to the F77 compiler',
             '-O2')
    opts.Add('CFORTRAN','Type of the Fortran-77 compiler (for cfortran.h)')
    opts.Add('F90','The Fortran-90 compiler')
    opts.Add('F90FLAGS','General options that are passed to the F90 compiler',
             '-O2')
    opts.Add('CFORTRAN90','Type of the Fortran-90 compiler (for cfortran.h)')
    opts.Add('F90MODSUFFIX','Suffix of Fortran-90 module interface files')
    opts.Add('MEXSUFFIX','Suffix for mex files')
    opts.Add('MEX','Matlab function compiler')
    opts.Add('MATLAB','Matlab interpreter')
    opts.Add('OCTAVE','Octave interpreter')
    opts.Add('MKOCTFILE','Octave function compiler')
    opts.Add('JAVAC','The Java compiler')
    opts.Add('JAVA_HOME','Location of jdk')
    opts.Add('MINESJTK','Location of edu_mines_jtk.jar')
    opts.Add('CUDA_TOOLKIT_PATH','Location of CUDA toolkit')
    opts.Add('NVCC','NVIDIA C compiler')
    
    return opts
