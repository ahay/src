
import sys

pkg = {}

pkg['libc'] = {'fedora':'glibc',
               'generic':'libc6-dev'}
pkg['ar']={'fedora':'binutils'}
pkg['libs'] = {'fedora':'glibc-headers',
               'cygwin':'sunrpc (Setup...Libs)'}
pkg['c99'] = {'fedora':'glibc-headers'}
pkg['xaw']={'fedora':'libXaw-devel',
            'generic':'libxaw7-dev'}

pkg['netpbm'] = {'fedora':'netpbm-devel',
                 'generic':'libnetpbm10-dev',
                 'darwin':'netpbm (fink)',
                 'cygwin':'libnetpbm-devel (Setup...Devel)'}
pkg['jpeg'] = {'fedora':'libjpeg-devel',
               'generic':'libjpeg62-dev'}
pkg['glew'] = {'generic':'libglew + libglew-dev',
               'fedora': 'glew + glew-devel'}


pkg['c++'] = {'fedora':'gcc-c++',
              'generic':'g++'}
# Used in checks for both f77 and f90
fortran = {'g77':'f2cFortran',
           'f77':'f2cFortran',
           'gfortran':'NAGf90Fortran',
           'gfc':'NAGf90Fortran',
           'f2c':'f2cFortran'}

pkg['f77'] = {'fedora':'gcc-gfortran',
              'generic':'g77'}

pkg['f90'] = {'fedora':'gcc-gfortran',
              'generic':'gfortran'}
pkg['octave'] = {'fedora':'octave',
                 'generic':'octave'}

pkg['mkoctave'] = {'fedora':'octave-devel',
                   'generic':'octave-headers'}
pkg['swig'] = {'fedora':'swig',
               'generic':'swig'}
pkg['numpy'] = {'fedora':'numpy',
                'generic':'python-scipy, python-numpy-dev'}
pkg['scipy'] = {'fedora':'scipy'}


def generate( env ):
    
    plat = {'OS': 'unknown',
            'distro': 'unknown',
            'arch': 'unknown',
            'cpu': 'unknown'}
    
    env['PLAT'] = plat
    env['PKG'] = pkg
    
    plat['OS'] = env.get('PLATFORM',sys.platform)

    # Check for distributions / OS versions
    try:
        from platform import architecture, uname
    except ImportError:
        
     # "platform" not installed. Python < 2.3
        # For each OS with Python < 2.3, should use specific
        # commands through os.system to find distro/version
        # Not known if what follows works everywhere:
        plat_nm = os.uname()[4]
        if plat_nm == 'x86_64':
            plat['arch'] = '64bit'
        elif plat_nm == 'i686':
            plat['arch'] = '32bit'
            
    else:
        
        plat['arch'] = architecture()[0]
        name = uname()[2].split('.')[-1]
        if plat['OS'] in ('linux', 'posix'):
            if name[:2] == 'fc':
                plat['distro'] = 'fedora'
            elif name[:2] == 'EL' or name[:2] == 'el':
                plat['distro'] = 'RedHat EL' # Redhat Enterprise
            elif name[-7:] == 'generic':
                plat['distro'] = 'generic' # Ubuntu
        elif plat['OS'] == 'sunos':
            if name[:2] == '10':
                plat['distro'] = '10' # Solaris 10
        elif plat['OS'] == 'darwin':
             plat['distro'] = uname()[2]
             plat['cpu'] = uname()[5] # i386 / powerpc
        elif plat['OS'] == 'irix':
             plat['distro'] = uname()[2]
        elif plat['OS'] in ('hp-ux', 'hpux'):
             plat['distro'] = uname()[2].split('.')[-2]
        del architecture, uname

def exists( env ):
    return 1



