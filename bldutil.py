import glob, os, re, string, py_compile

try:
    import configure
except:
    import rsfconf as configure

# The following adds all SCons SConscript API to the globals of this module.
import SCons
version = map(int,string.split(SCons.__version__,'.')[:3])
if version[0] == 1  or \
   version[1] >= 97 or \
   (version[1] == 96 and version[2] >= 90):
    from SCons.Script import *
else:  # old style
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

################################################################################

# Constants used in multiple functions

__py_success = 0 # user-defined
__include = re.compile(r'#include\s*\"([^\"]+)\.h\"')

################################################################################

def __docmerge(target=None,source=None,env=None):
    outfile = target[0].abspath
    out = open(outfile,'w')
    out.write('import rsfdoc\n\n')
    for src in map(str,source):
        inp = open(src,'r')
        for line in inp.readlines():
                out.write(line)
        inp.close()
    alias = env.get('alias',{})
    for prog in alias.keys():
        out.write("rsfdoc.progs['%s']=%s\n" % (prog,alias[prog]))
    out.close()
    py_compile.compile(outfile,outfile+'c')
    return __py_success

def __pycompile_emit(target, source, env):
    target.append(str(target[0])+'c')
    return target,source 

Docmerge = Builder(action=Action(__docmerge,varlist=['alias']),
                   emitter=__pycompile_emit)

################################################################################

def Debug():
    'Environment for debugging'
    env = Environment()
    srcroot = os.environ.get('RSFSRC', '../..')
    opts = configure.options(os.path.join(srcroot,'config.py'))
    opts.Update(env)
    env['CCFLAGS'] = env.get('CCFLAGS','').replace('-O2','-g')
    if  env['PLATFORM'] == 'sunos':
        env['CCFLAGS'] = string.replace(env.get('CCFLAGS',''),'-xO2','-g')
    env['F90FLAGS'] = string.replace(env.get('F90FLAGS',''),'-O2','-g')
    env.SConsignFile(None)
    env.Append(BUILDERS={'RSF_Include':configure.Header,
                         'RSF_Place':configure.Place},
               SCANNERS=[configure.Include])
    return env

################################################################################

def depends(env,list,file):
    'Find dependencies for C'
    filename = env.File(file+'.c').abspath
    # replace last occurence of build/
    last = filename.rfind('build/')
    if last >= 0:
        filename = filename[:last] + filename[last+6:]
    fd = open(filename,'r')
    for line in fd.readlines():
        for inc in __include.findall(line):
            if inc not in list and inc[0] != '_':
                list.append(inc)
                depends(env,list,inc)
    fd.close()

################################################################################

def depends90(env,list,file):
    'Find dependencies for Fortran-90'

    include90 = re.compile(r'^[^!]*use\s+(\S+)')
    filename = env.File(file+'.f90').abspath
    # replace last occurence of build/
    last = filename.rfind('build/')
    if last >= 0:
        filename = filename[:last] + filename[last+6:]
    fd = open(filename,'r')
    for line in fd.readlines():
        for inc in include90.findall(line):
            if inc not in list and inc != 'rsf':
                list.append(inc)
                depends90(env,list,inc)
    fd.close()

################################################################################

def chk_exists(prog, ext='c', mainprog=True):
    'Check if file corresponding to program name exists'

    ext = ext.lstrip('.') # In case the user put a dot
    prognm = prog + '.' + ext
    if mainprog:
        prognm = 'M' + prognm
    try:
        assert os.path.isfile(prognm)
    except:
        msg = 'Missing file: ' + os.path.join(os.getcwd(),prognm)
        configure.stderr_write(msg,'yellow_on_red')
        sys.exit(configure.unix_failure)

################################################################################

def build_install_c(env, progs_c, bindir, glob_build, bldroot):
    'Build and install C programs'

    env.Prepend(CPPPATH=[os.path.join(bldroot,'include')],
            LIBPATH=[os.path.join(bldroot,'lib')],
            LIBS=['rsf'])

    if glob_build:
        dir = string.replace(os.getcwd(),'/build','') # aka RSFSRC/user/$USER
        src = map(os.path.basename,glob.glob(os.path.join(dir,'[a-z]*.c')))
    else:
        src = glob.glob('[a-z]*.c')

    for source in src:
        inc = env.RSF_Include(source,prefix='')
        obj = env.StaticObject(source)
        env.Depends(obj,inc)

    mains_c = Split(progs_c)
    for prog in mains_c:
        if not glob_build:
            chk_exists(prog)
        sources = ['M' + prog]
        depends(env, sources, 'M'+prog)
        prog = env.Program(prog, map(lambda x: x + '.c',sources))
        if glob_build:
            env.Install(bindir,prog)

    if glob_build:
        docs_c = map(lambda prog: env.Doc(prog,'M'+prog),mains_c)
    else:
        docs_c = None

    return docs_c

################################################################################

def build_install_f90(env, progs_f90, bindir, api, bldroot, glob_build):
    'Build and install Fortran90 programs'

    mains_f90 = Split(progs_f90)

    if 'f90' in api:

        F90 = env.get('F90')
        assert F90 != None # The configure step should have found the compiler
        F90base = os.path.basename(F90)
        if F90base[:8] == 'gfortran' or F90base[:3] == 'gfc':
            env.Append(F90FLAGS=' -J${SOURCE.dir}')
        elif F90base == 'ifort':
            env.Append(F90FLAGS=' -module ${SOURCE.dir}')

        env.Prepend(LIBS='rsff90', # order matters when linking
                    F90PATH=os.path.join(bldroot,'include'))

        for prog in mains_f90:
            if not glob_build:
                chk_exists(prog, 'f90')
            obj_dep = []
            sources = ['M' + prog]
            depends90(env,sources,'M'+prog)
            for f90_src in sources:
                obj = env.StaticObject(f90_src+'.f90')
                # SCons mistakenly treats ".mod" files as ".o" files, and
                # tries to build them the same way (which fails). So we
                # explicitly keep just the ".o" files as dependencies:
                for fname in obj:
                    if os.path.splitext(fname.__str__())[1] == '.o':
                        obj_dep.append(fname)
            # Using obj_dep instead of the list of sources because when two
            # mains used the same module, object files for the module were 
            # created in both places, hence endless "double-define" warnings
            prog = env.Program(prog, obj_dep, LINK=F90)
            if glob_build:
                env.Install(bindir,prog)

    else: # Put in a placeholder
        for prog in mains_f90:
            prog = env.RSF_Place('sf'+prog,None,package='Fortran90+API=F90')
        if glob_build:
            env.Install(bindir,prog)

    if glob_build:
        docs_f90 = map(lambda prog: env.Doc(prog,'M'+prog+'.f90',lang='f90'),
               mains_f90)
    else:
        docs_f90 = None

    return docs_f90

################################################################################

def install_py_mains(env, progs_py, bindir):
    'Copy Python programs to bindir, generate list of self-doc files'

    mains_py = Split(progs_py)
    for prog in mains_py:
        env.InstallAs(os.path.join(bindir,'sf'+prog),'M'+prog+'.py')

    # Self-doc
    user = os.path.basename(os.getcwd())
    main = 'sf%s.py' % user
    docs_py = map(lambda prog: env.Doc(prog,'M'+prog+'.py',lang='python'),
           mains_py)

    return docs_py

################################################################################

def install_py_modules(env, py_modules, libdir):
    'Compile Python modules and install to libdir/rsfuser'

    rsfuser = os.path.join(libdir,'rsfuser')
    for module in Split(py_modules):
        env.RSF_Pycompile(module+'.pyc',module+'.py')
        env.Install(rsfuser,module+'.pyc')

################################################################################

def install_self_doc(env, libdir, docs_c=None, docs_py=None, docs_f90=None):

    docs = []
    if docs_c != None:
        docs += docs_c
    if docs_py != None:
        docs += docs_py
    if docs_f90 != None:
        docs += docs_f90

    env.Depends(docs,'#/framework/rsfdoc.py')	

    user = os.path.basename(os.getcwd())
    main = 'sf%s.py' % user
    doc = env.RSF_Docmerge(main,docs)
    env.Install(libdir,doc)

################################################################################

def add_ext_static_lib(env, libnm, root=os.environ.get('RSFROOT')):
    
    env['LIBS'].append(File(os.path.join(root,'lib','lib'+libnm+'.a')))
    env['CPPPATH'].append(os.path.join(root,'include'))
