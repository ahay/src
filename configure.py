import sys, os
import string, re

# The following adds all SCons SConscript API to the globals of this module.
import SCons.Script.SConscript
globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

toheader = re.compile(r'\n\n((?:[^\n]|\n[^\n])+)\n'
                      '\/\*(\^|\<(?:[^>]|\>[^*]|\>\*[^/])*\>)\*\/')

def header(target=None,source=None,env=None):
    inp = open(str(source[0]),'r')
    text = string.join(inp.readlines(),'')
    inp.close()
    file = str(target[0])
    prefix = env.get('prefix','')
    define = prefix + string.translate(os.path.basename(file),
                                       string.maketrans('.','_'))
    out = open(file,'w')    
    out.write('#ifndef _' + define + '\n')
    out.write('#define _' + define + '\n\n')
    for extract in toheader.findall(text):
        if extract[1] == '^':
            out.write(extract[0]+'\n\n')
        else:
            out.write(extract[0]+';\n')
            out.write('/*'+extract[1]+'*/\n\n')
    out.write('#endif\n')
    out.close()
    return 0

Header = Builder (action = Action(header,varlist=['prefix']),
                  src_suffix='.c',suffix='.h')

include = re.compile(r'#include\s*\"([^\"]+)\.h\"')

def depends(env,list,file):
    filename = string.replace(env.File(file+'.c').abspath,'build/','',1)
    fd = open(filename,'r')
    for line in fd.readlines():
        for inc in include.findall(line):
            if inc not in list and inc[0] != '_':
                list.append(inc)
                depends(env,list,inc)
    fd.close()

def included(node,env,path):
    file = os.path.basename(str(node))
    file = re.sub('\.[^\.]+$','',file)
    contents = node.get_contents()
    includes = include.findall(contents)
    if file in includes:
        includes.remove(file)
    return map(lambda x: x + '.h',includes)

Include = Scanner(name='Include',function=included,skeys=['.c'])

def check_all(context):
    cc(context)
    ar(context)
    libs(context)
    api = string.split(string.lower(context.env.get('API','')),',')
    if 'c++' in api:
        cxx(context)
    if 'fortran' in api:
        f77(context)
    if 'fortran-90' in api:
        f90(context)

def libs(context):
    context.Message("checking libraries ... ")
    LIBS = string.split(context.env.get('LIBS','m'))
    if sys.platform[:5] == 'sunos':
        LIBS.append('nsl')
    elif sys.platform[:6] == 'cygwin':
        LIBS.append('rpc')
    text = '''
    #include <rpc/types.h>
    #include <rpc/xdr.h>
    int main(int argc,char* argv[]) {
    return 0;
    }
    '''
    res = context.TryRun(text,'.c')
    if res[0]:
        context.Result(str(LIBS))
        context.env['LIBS'] = LIBS
    else:
        context.Result(0)
        context.Message("Please install RPC libraries.\n")
        sys.exit(1)

def ar(context):
    context.Message("checking ar ... ")
    AR = context.env.get('AR',WhereIs('ar'))
    if AR:
        context.Result(AR)
        context.env['AR'] = AR
    else:
        context.Result(0)
        sys.exit(1)

def cc(context):
    context.Message("checking C compiler ... ")
    CC = context.env.get('CC',WhereIs('gcc'))
    if CC:
        context.Result(CC)   
    else:
        context.Result(0)
        sys.exit(1)
    text = '''
    int main(int argc,char* argv[]) {
    return 0;
    }
    '''
    context.Message("checking if %s works ... " % CC)
    res = context.TryRun(text,'.c')
    context.Result(res[0])
    if not res[0]:
        sys.exit(1)
    if CC[-3:]=='gcc':
        oldflag = context.env.get('CCFLAGS')
        for flag in ('-std=gnu99 -Wall -pedantic',
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
    elif sys.platform[:5] == 'sunos':
        context.env['CCFLAGS'] = string.replace(context.env.get('CCFLAGS',''),
                                                '-O2','-xO2')

def cxx(context):
    context.Message("checking C++ compiler ... ")
    CXX = context.env.get('CXX')
    if CXX:
        context.Result(CXX)   
    else:
        context.Result(0)
        return
    text = '''
    #include <valarray>
    int main(int argc,char* argv[]) {
    return 0;
    }
    '''
    context.Message("checking if %s works ... " % CXX)
    res = context.TryRun(text,'.cc')
    context.Result(res[0])
    if not res[0]:
        del context.env['CXX']
        return
    if CXX == 'g++':
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

fortran = {'g77':'f2cFortran',
           'f2c':'f2cFortran'}

def f77(context):
    context.Message("checking F77 compiler ... ")
    F77 = context.env.get('F77')
    if F77:
        context.Result(F77)
    else:
        context.Result(0)
        return
    if os.path.basename(F77) == 'ifc':
        intel(context)
        context.env.Append(F77FLAGS=' -Vaxlib')
    text = '''      program Test
      stop
      end
      '''
    context.Message("checking if %s works ... " % F77)
    oldlink = context.env.get('LINK')
    context.env['LINK'] = F77
    res = context.TryRun(text,'.f')
    context.env['LINK'] = oldlink
    context.Result(res[0])
    if not res[0]:
        sys.stderr.write("No working F77 compiler detected.\n")
        del context.env['F77']
        sys.exit(1)
    cfortran = fortran.get(os.path.basename(F77),'NAGf90Fortran')
    context.env['CFORTRAN'] = cfortran 
    context.Message("checking %s type for cfortran.h ... " % F77)
    context.Result(cfortran)
    
def f90(context):
    context.Message("checking F90 compiler ... ")
    F90 = context.env.get('F90')
    if not F90:
        compilers = ['f90','f95','xlf90','pgf90','ifc','pghpf']
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
        context.Result(0)
        return
    if os.path.basename(F90) == 'ifc':
        intel(context)
        context.env.Append(F90FLAGS=' -Vaxlib')
    load_f90(context.env) 
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
    res2 = context.TryRun(main,'.f90')
    context.env['LINK'] = oldlink
    context.Result(res1 and res2[0])
    if not res1 or not res2[0]:
        sys.stderr.write("No working F90 compiler detected.\n")
        del context.env['F90']
        sys.exit(1)
    base = os.path.basename(F90)
    context.Message("checking %s type for cfortran.h ... " % base)
    cfortran = fortran.get(base,'NAGf90Fortran')
    context.env['CFORTRAN90'] = cfortran 
    context.Result(cfortran)
    context.Message("checking F90 module extension ... ")
    f90module = re.compile(r'(?:testf90|TESTF90)(\.\w+)$')
    suffix = ''
    for file in os.listdir(os.getcwd()):
        gotit = f90module.match(file)
        if gotit:
            suffix = gotit.group(1)
            os.remove(file)
            break
    context.env['F90MODSUFFIX'] = suffix
    context.Result(suffix)

def load_f90(env):
    env['F90COM']   = '$F90 $F90FLAGS $_F90INCFLAGS -c -o $TARGET $SOURCES'
    env['SHF90COM'] = '$F90 $F90FLAGS $_F90INCFLAGS -c -o $TARGET $SOURCES'
    static_obj, shared_obj = createObjBuilders(env)    
    F90Action = Action("$F90COM")
    ShF90Action = Action("$SHF90COM")
    static_obj.add_action('.f90', F90Action)
    shared_obj.add_action('.f90', ShF90Action)

def intel(context):
    '''Trying to fix wierd intel setup.'''
    libdirs = string.split(os.environ.get('LD_LIBRARY_PATH'),':')
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
    opts.Add('CC','The C compiler')
    opts.Add('CCFLAGS','General options that are passed to the C compiler',
             '-O2')
    opts.Add('CPPPATH',
             'The list of directories that the C preprocessor will search')
    opts.Add('LIBPATH',
             'The list of directories that will be searched for libraries')
    opts.Add('LIBS',
             'The list of libraries that will be linked with executables')
    opts.Add('PROGPREFIX','The prefix used for executable file names','sf')
    opts.Add('API','Support for additional languages (possible values: c++, fortran, fortran-90, python)')
    opts.Add('CXX','The C++ compiler')
    opts.Add('CXXFLAGS','General options that are passed to the C++ compiler',
             '-O2')
    opts.Add('F77','The Fortran-77 compiler')
    opts.Add('F77FLAGS','General options that are passed to the F77 compiler',
             '-O2')
    opts.Add('CFORTRAN','Type of the Fortran-77 compiler (for cfortran.h)')
    opts.Add('F90','The Fortran-90 compiler')
    opts.Add('F90FLAGS','General options that are passed to the F90 compiler',
             '-O2')
    opts.Add('CFORTRAN90','Type of the Fortran-90 compiler (for cfortran.h)')
    opts.Add('F90MODSUFFIX','Suffix of Fortran-90 module interface files')

def merge(target=None,source=None,env=None):
    sources = map(str,source)
    local_include = re.compile(r'\s*\#include\s*\"([\w\.]+)')
    includes = []
    for src in sources:
        if src in includes:
            continue
        inp = open(src,'r')
        for line in inp.readlines():
            match = local_include.match(line)            
            if match:
                other = match.group(1)
                if not other in includes:
                    includes.append(os.path.join(os.path.dirname(src),other))
        inp.close()
        includes.append(src)
    out = open(str(target[0]),'w')
    for src in includes:
        inp = open(src,'r')
        for line in inp.readlines():
            if not local_include.match(line):
                out.write(line)
        inp.close()
    out.close()
    return 0

docmerge = '''echo "import rsfdoc" > $TARGET
echo "" >> $TARGET
cat $SOURCES >> $TARGET'''

def docextra(docmerge,source,copy):
    return docmerge + '''
    echo rsfdoc.progs[\\'%s\\']=%s >> $TARGET''' % (copy,source)

#	$Id$

