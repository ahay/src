import sys
sys.path.append('../../framework')
from bldutil import *

class HuiSconsTargets:
    '''Simple wrapper for convinent building'''
    docs = []
    has_lapack = False
    has_nvcc = False
    def __init__(self, cfiles=None, ccfiles=None, cufiles=None, cmpifiles=None):
        self.c = cfiles
        self.cc = ccfiles
        self.cu = cufiles
        self.c_mpi = cmpifiles

# -----------------------------------------------------------------------------
    @classmethod
    def install_docs(cls,env,libdir,glob_build):
        if glob_build and cls.docs:
            env.Depends(cls.docs, '#/framework/rsf/doc.py')
            user = os.path.basename(os.getcwd())
            main = 'sf%s.py' % user
            doc = env.RSF_Docmerge(main,cls.docs)
            env.Install(libdir,doc)

# -----------------------------------------------------------------------------

    def build_c(self, env, glob_build, srcroot, bindir, libdir, pkgdir):
        bldroot = '../..'
        if not glob_build:
            bldroot = env.get('RSFROOT', os.environ.get('RSFROOT',sys.prefix))
            # SConscript(os.path.join(srcroot, 'api', 'c', 'SConstruct'))
        if not self.c:
            docs_c = None
        else:
            docs_c = build_install_c(env, self.c, srcroot, bindir, libdir, glob_build, bldroot)
        if glob_build:
            HuiSconsTargets.docs.append(docs_c)

# -----------------------------------------------------------------------------

    def build_cc(self, env, glob_build, srcroot, bindir, libdir, pkgdir):
        if not self.cc:
            docs_cc = None
        else:
            docs_cc = self.build_install_cc(env, self.cc, srcroot, bindir, libdir, glob_build)
        if glob_build:
            HuiSconsTargets.docs.append(docs_cc)

# -----------------------------------------------------------------------------
    def build_cu(self, env, glob_build, srcroot, bindir, libdir, pkgdir):
        # Assume using C API rather than CPP API
        if not self.cu:
            docs_cu = None
        else:
            docs_cu = self.build_install_cu(env, self.cu, srcroot, bindir, libdir, glob_build)
        if glob_build:
            HuiSconsTargets.docs.append(docs_cu)

# -----------------------------------------------------------------------------
    def build_c_mpi(self, env, glob_build, srcroot, bindir, libdir, pkgdir):
        bldroot = '../..'
        if not glob_build:
            bldroot = env.get('RSFROOT', os.environ.get('RSFROOT',sys.prefix))
        if not self.c_mpi:
            docs_c_mpi = None
        else:
            docs_c_mpi = build_install_c_mpi(env, self.c_mpi, srcroot, bindir, glob_build, bldroot)
        if glob_build:
            HuiSconsTargets.docs.append(docs_c_mpi)

# -----------------------------------------------------------------------------

    # @classmethod
    def build_install_cc(self,env, progs_cc, srcroot, bindir, libdir, glob_build):
        'Build and install CC programs'
        dynlib = env.get('DYNLIB','')
        env.Prepend(CPPPATH=[os.path.join(srcroot,'include')],
                    LIBPATH=[os.path.join(srcroot,'lib')],
                    LIBS=[dynlib+'rsf'])
        if 'c++' in env.get('API',[]):
            env.Prepend(LIBS=[dynlib+'rsf++'])
        mains_cc = Split(progs_cc)
        for prog in mains_cc:
            sources = ['M' + prog]
            if self.has_lapack:
                print HuiSconsTargets.has_lapack
                prog = env.Program(prog,map(lambda x: x + '.cc',sources))
            else:
                prog = env.RSF_Place('sf'+prog,None,var='LAPACK',package='lapack')

            if glob_build:
                install = env.Install(bindir,prog)
                if dynlib and env['PLATFORM'] == 'darwin':
                    env.AddPostAction(install,
                    'install_name_tool -change '
                    'build/api/c/libdrsf.dylib '
                    'build/api/c++/libdrsf++.dylib '
                    '%s/libdrsf.dylib %s' % (libdir,install[0]))
        if glob_build:
            docs_cc = map(lambda prog: env.Doc(prog,'M'+prog+'.cc'),mains_cc)
        else:
            docs_cc = None

        return docs_cc

# -----------------------------------------------------------------------------

    # @classmethod
    def build_install_cu(self,env, progs_cu, srcroot, bindir, libdir, glob_build):
        'Build and install CUDA programs'
        dynlib = env.get('DYNLIB','')

        env.Prepend(CPPPATH=[os.path.join(srcroot,'include')],
                    LIBPATH=[os.path.join(srcroot,'lib')],
                    LIBS=[dynlib+'rsf'])

        mains_cu = Split(progs_cu)
        for prog in mains_cu:
            if not glob_build:
                chk_exists(prog,'cu')

            sources = env.Command(prog+'.cpp','M'+prog+'.cu','cp $SOURCE $TARGET')
            if self.has_nvcc:
                prog = env.Program(prog,sources)
            else:
                prog = env.RSF_Place('sf'+prog,None,var='CUDA Toolkit',package='CUDA Toolkit')

            if glob_build:
                install = env.Install(bindir,prog)
                if dynlib and env['PLATFORM'] == 'darwin':
                    env.AddPostAction(install,
                    'install_name_tool -change '
                    'build/api/c/libdrsf.dylib '
                    '%s/libdrsf.dylib %s' % (libdir,install[0]))
        if glob_build:
            docs_cu = map(lambda prog: env.Doc(prog,'M'+prog+'.cu'),mains_cu)
        else:
            docs_cu = None

        return docs_cu

# -----------------------------------------------------------------------------


