import os, sys, string, re
import configure, rsfdoc

root = os.environ.get('RSFROOT',os.getcwd())

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')
docdir = os.path.join(root,'doc')

pydir = libdir
#for path in sys.path:
#    if os.path.isdir(os.path.join(path,'SCons')):
#        pydir = path

env = Environment()

##########################################################################
# CONFIGURATION
##########################################################################

opts = Options('config.py')
configure.options(opts)
opts.Add('RSFROOT','RSF installation root',root)
opts.Update(env)

if not os.path.isfile('config.py'):
    conf = Configure(env,custom_tests={'CheckAll':configure.check_all})
    conf.CheckAll()
    env = conf.Finish()
    
Help(opts.GenerateHelpText(env))
opts.Save('config.py',env)

config = env.Command('config.py','configure.py','')
env.Precious(config)
env.InstallAs(os.path.join(libdir,'rsfconfig.py'),'config.py')
env.InstallAs(os.path.join(libdir,'rsfconf.py'),'configure.py')
Clean(config,['#/config.log','#/.sconf_temp','configure.pyc'])
env.Alias('config',config)

##########################################################################
# SELF DOCUMENTATION
##########################################################################
Doc = Builder (action = Action(rsfdoc.selfdoc,varlist=['rsfprefix']),
               src_suffix='.c',suffix='.py')
env.Append(BUILDERS = {'Doc' : Doc})

##########################################################################
# FILT BUILD
##########################################################################
env.Prepend(CPPPATH=['../../include'],
            LIBPATH=['../../filt/lib'],
            LIBS=['rsf'])

Export('env')
dirs = ('lib','main','proc','imag')

Default('build/include')
for dir in map(lambda x: os.path.join('filt',x), dirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct')
    Default(build)

if env.has_key('F90'):
    modsuffix = env.get('F90MODSUFFIX')
    F90 = os.path.basename(env.get('F90'))
    if modsuffix:
        if 'ifort'==F90:
            mod = 'rsf'+modsuffix
        else:
            mod = 'RSF'+modsuffix
        env.Install(incdir,mod)
        Clean(build,mod)
    elif 'ifc' == F90: # old Intel compiler quirks
        env.InstallAs(os.path.join(incdir,'rsf.pc'),'work.pc')
        env.Install(incdir,'rsf.d')
        Clean(build,['rsf.d','work.pc'])

##########################################################################
# PLOT BUILD
##########################################################################
env.Prepend(LIBPATH=['../../plot/lib'],LIBS=['rsfplot'])

Export('env')
pdirs = ('lib','main','test')

Default('build/include')
for dir in map(lambda x: os.path.join('plot',x), pdirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct')
    Default(build)

##########################################################################
# PYTHON MODULES
##########################################################################

for src in ('doc','proj','prog','tex'):
    py = "rsf%s.py"% src
    pyc = py + 'c'
    env.Install(pydir,py)
    Clean(os.path.join(pydir,py),[os.path.join(pydir,pyc),pyc])
env.Install(bindir,'sfdoc')
env.Install(bindir,'sftour')

##########################################################################
# INSTALLATION
##########################################################################

env.Alias('install',[bindir,pydir,libdir,incdir,docdir])

use = os.path.join(pydir,'rsfuse.py')
env.Command(use,None,action=Action(rsfdoc.use))
Depends(use,map(lambda x: os.path.join(libdir,'sf'+x+'.py'),dirs[1:]))
Depends(use,os.path.join(libdir,'sfplot.py'))
Depends(use,os.path.join(libdir,'vpplot.py'))

index = os.path.join(docdir,'index.html')
env.Command(index,None,'PYTHONPATH=%s %s sfdoc -w %s' % 
           (libdir,WhereIs('python'),docdir))
Depends(index,use)

# 	$Id$

