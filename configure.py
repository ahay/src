import sys, os
import string, re
from SCons.Util import WhereIs
from SCons.Tool import createObjBuilders
from SCons.Action import Action
from SCons.Defaults import StaticCheckSet, SharedCheckSet

def check_all(context):
    cc(context)
    f77(context)
    f90(context)

def cc(context):
    context.Message("Checking CC compiler ... ")
    CC = context.env.get('CC')
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
    context.Message("Checking if %s works ... " % CC)
    res = context.TryRun(text,'.c')
    context.Result(res[0])
    if not res:
        sys.exit(1)
    if CC == 'gcc':
        oldflag = context.env.get('CCFLAGS')
        for flag in ('-std=gnu99 -Wall -pedantic',
                     '-std=gnu9x -Wall -pedantic',
                     '-Wall -pedantic'):
            context.Message("Checking if gcc accepts '%s' ... " % flag)
            context.env['CCFLAGS'] = oldflag + ' ' + flag
            res = context.TryCompile(text,'.c')
            context.Result(res)
            if res:
                break
        if not res:
            context.env['CCFLAGS'] = oldflag

fortran = {'g77':'f2cFortran',
           'f2c':'f2cFortran'}

def f77(context):
    context.Message("Checking F77 compiler ... ")
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
    context.Message("Checking if %s works ... " % F77)
    oldlink = context.env.get('LINK')
    context.env['LINK'] = F77
    res = context.TryRun(text,'.f')
    context.env['LINK'] = oldlink
    context.Result(res[0])
    if not res:
        sys.stderr.write("No working F77 compiler detected.\n")
        del context.env['F77']
        return
    cfortran = fortran.get(os.path.basename(F77),'NAGf90Fortran')
    context.env['CFORTRAN'] = cfortran 
    context.Message("Checking %s type for cfortran.h ... " % F77)
    context.Result(cfortran)
    
def f90(context):
    context.Message("Checking F90 compiler ... ")
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
    module = '''module test
    end module test
    '''
    context.Message("Checking if %s works ... " % F90)
    oldlink = context.env.get('LINK')
    context.env['LINK'] = F90
    res1 = context.TryCompile(module,'.f90')
    res2 = context.TryRun(main,'.f90')
    context.env['LINK'] = oldlink
    context.Result(res1 and res2[0])
    if not res1 or not res2:
        sys.stderr.write("No working F90 compiler detected.\n")
        del context.env['F90']
        return
    base = os.path.basename(F90)
    cfortran = fortran.get(base,'NAGf90Fortran')
    context.env['CFORTRAN90'] = cfortran 
    context.Message("Checking %s type for cfortran.h ... " % base)
    context.Result(cfortran)

def load_f90(env):
    env['F90COM']   = '$F90 $F90FLAGS $_F90INCFLAGS -c -o $TARGET $SOURCES'
    env['SHF90COM'] = '$F90 $F90FLAGS $_F90INCFLAGS -c -o $TARGET $SOURCES'
    static_obj, shared_obj = createObjBuilders(env)    
    F90Action = Action([StaticCheckSet,"$F90COM"])
    ShF90Action = Action([SharedCheckSet,"$SHF90COM"])
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

docmerge = '''echo "import os, sys" > $TARGET
echo "sys.path.append(os.environ.get('RSFROOT'))" >> $TARGET
echo "import rsfdoc" >> $TARGET
echo "" >> $TARGET
cat $SOURCES >> $TARGET'''

