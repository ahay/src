import os, sys, string, re
import configure, rsfdoc

here = os.getcwd()
root = os.environ.setdefault('RSFROOT',here)

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')

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

config = env.Command('config.py','configure.py',"")
env.Precious(config)
env.InstallAs(os.path.join(libdir,'rsfconfig.py'),'config.py')
env.InstallAs(os.path.join(libdir,'rsfconf.py'),'configure.py')
Clean(config,['#/config.log','#/.sconf_temp','configure.pyc'])
env.Alias('config',config)

##########################################################################
# PYTHON MODULES
##########################################################################
for src in ('doc','proj','prog'):
    py = "rsf%s.py"% src
    pyc = py + 'c'
    env.Install(libdir,py)
    Clean(os.path.join(libdir,py),[os.path.join(libdir,pyc),pyc])
env.Install(bindir,'sfdoc')

##########################################################################
# SELF DOCUMENTATION
##########################################################################
Doc = Builder (action = Action(rsfdoc.selfdoc),src_suffix='.c',suffix='.py')
env.Append(BUILDERS = {'Doc' : Doc})

##########################################################################
# FILT BUILD
##########################################################################
env.Append(CPPPATH=['../../include'],
           LIBPATH=['../../filt/lib'],
           LIBS=['rsf','m'])

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
    env.Install(incdir,'RSF'+modsuffix)
    if re.search(r'ifc$',env.get('F90')) and modsuffix == '.d':
        # old Intel compiler quirks
        env.InstallAs(os.path.join(incdir,'rsf.pc'),'work.pc')

##########################################################################
# PLOT BUILD
##########################################################################
env.Prepend(LIBPATH=['../../plot/lib'],LIBS=['rsfplot'])

Export('env')
dirs = ('lib','main')

Default('build/include')
for dir in map(lambda x: os.path.join('plot',x), dirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct')
    Default(build)

##########################################################################
# INSTALLATION
##########################################################################

env.Alias('install',[bindir,libdir,incdir])
