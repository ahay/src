import os, sys, string, re
import configure, rsfdoc

if os.environ.has_key('RSFROOT'):
    root = os.environ['RSFROOT']
else:
    root = os.getcwd()

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

config = env.Command('config.py','configure.py',"")
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
env.Append(CPPPATH=['../../include'],
           LIBPATH=['../../filt/lib'],
           LIBS=['rsf','m'])

if sys.platform[:5] == 'sunos':
    env.Append(LIBS=['nsl'])
    if not env['CC'].endswith('gcc'):
        # Sun's native compiler 
        env['CCFLAGS']='-xO2'

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
    if modsuffix:
        env.Install(incdir,'RSF'+modsuffix)
        Clean(build,'RSF'+modsuffix)
    elif re.search(r'ifc$',env.get('F90')): # old Intel compiler quirks
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

for src in ('doc','proj','prog'):
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

env.Command(docdir,None,'RSFROOT=%s ./sfdoc -w $TARGET' % root)
Depends(docdir,use)

# 	$Id: SConstruct,v 1.30 2004/06/23 18:29:51 fomels Exp $	
