EnsureSConsVersion(0, 96)

import os
import configure

env = Environment()

root = os.environ.get('RSFROOT',os.getcwd())

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')
docdir = os.path.join(root,'doc')
mandir = os.path.join(root,'man')

##########################################################################
# CONFIGURATION
##########################################################################
opts = configure.options('config.py')
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
env.InstallAs(os.path.join(libdir,'rsfconf.pyc'),'configure.pyc')
Clean(config,['#/config.log','#/.sconf_temp','configure.pyc'])
env.Alias('config',config)

##########################################################################
# CUSTOM BUILDERS
##########################################################################

env.Append(BUILDERS={'Include':configure.Header,
                     'Place':configure.Place,
                     'Pycompile':configure.Pycompile,
                     'Docmerge':configure.Docmerge},
           SCANNERS=[configure.Include])

##########################################################################
# FRAMEWORK BUILD
##########################################################################

system = filter(lambda x: x[0] != '.', os.listdir('system'))
user = filter(lambda x: x[0] != '.' and x != 'nobody', os.listdir('user'))

SConscript(dirs='framework',name='SConstruct',
           exports='env root bindir libdir docdir mandir system user')

##########################################################################
# API BUILD
##########################################################################
api = env.get('API',[])
api.insert(0,'c')

Default('build/include')
Default('build/lib')
for dir in map(lambda x: os.path.join('api',x), api):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',exports='env root libdir incdir')
    Default(build)


##########################################################################
# SYSTEM BUILD
##########################################################################
for dir in map(lambda x: os.path.join('system',x), system):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',exports='env root bindir libdir')
    Default(build)

##########################################################################
# USER BUILD
##########################################################################
for dir in map(lambda x: os.path.join('user',x), user):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',exports='env root bindir libdir')
    Default(build)

##########################################################################
# PLOT BUILD
##########################################################################
pdirs = ('lib','main','test','opengl')

for dir in map(lambda x: os.path.join('plot',x), pdirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',
               exports='env root libdir bindir incdir')
    Default(build)

##########################################################################
# PENS BUILD
##########################################################################
pdirs = ('fonts','include','utilities','genlib','main','docs')

for dir in map(lambda x: os.path.join('pens',x), pdirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',
               exports='env root incdir libdir bindir')
    Default(build)

##########################################################################
# SU BUILD
##########################################################################
sudirs = ('lib','main','plot')

for dir in map(lambda x: os.path.join('su',x), sudirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',
               exports='env root libdir bindir incdir')
    Default(build)

##########################################################################
# INSTALLATION
##########################################################################

rsfuser = os.path.join(libdir,'rsfuser')
env.Install(rsfuser,'__init__.py')

env.Alias('install',[incdir,bindir,libdir,rsfuser,docdir,mandir])
env.Clean('install', rsfuser)

# 	$Id$
