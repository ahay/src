EnsureSConsVersion(0, 96)

import os, sys, glob

env = Environment()

root = os.environ.get('RSFROOT',os.getcwd())

pythonpath = os.environ.get('PYTHONPATH')
framework = os.path.join(os.getcwd(),'build/framework')

if pythonpath:
    os.environ['PYTHONPATH'] = ':'.join([pythonpath,framework])
else:
    os.environ['PYTHONPATH'] = framework
sys.path = [framework,'framework'] + sys.path

import bldutil, configure, setenv

env.InstallAs('#/build/framework/rsf/conf.py','framework/configure.py') 

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')
docdir = os.path.join(root,'doc')
spcdir = os.path.join(root,'spec')
mandir = os.path.join(root,'share','man')
pkgdir = setenv.get_pkgdir(root)
etcdir = os.path.join(root, 'etc', 'madagascar')

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

Help(opts.GenerateHelpText(env,cmp))
opts.Save('config.py',env)
config = env.Command('config.py','framework/configure.py','')
env.Precious(config)

Clean(config,['#/config.log','#/.sconf_temp',
              '#/framework/bldutil.pyc',
              '#/framework/configure.pyc',
              '#/framework/setenv.pyc',
              '#/framework/rsf/doc.pyc',
              '#/framework/rsf/path.pyc'])
env.Alias('config',config)
env.Install(etcdir,'config.py') 

# ----------- Environment variable setup scripts -----------

for sh in ('sh','csh'):
    shrc  = env.Command('env.'+sh,  '', setenv.shell_script,
                        varlist=('shell','DATAPATH'),
                        shell=sh)
    env.Alias('config',shrc)
    env.Install(etcdir,shrc)

##########################################################################
# CUSTOM BUILDERS
##########################################################################

env.Append(BUILDERS={'RSF_Include':bldutil.Header,
                     'RSF_Place':bldutil.Place,
                     'RSF_Pycompile':bldutil.Pycompile,
                     'RSF_Docmerge':bldutil.Docmerge},
           SCANNERS=[bldutil.Include])

##########################################################################
# FRAMEWORK BUILD
##########################################################################

system = filter(lambda x: x[0] != '.', os.listdir('system'))
user = filter(lambda x: x[0] != '.' and x != 'nobody', os.listdir('user'))
# Avoid crashing when user places some files in RSFSRC/user
user = filter(lambda x: os.path.isdir(os.path.join('user',x)), user)
dotproj = glob.glob('book/*/*/*/.rsfproj')

frame_exports = 'env bindir libdir pkgdir docdir spcdir mandir system user dotproj'

for dir in map(lambda x: os.path.join('framework',x),Split('rsf doc ptools')):
    build = os.path.join('build',dir)
    BuildDir(build,dir)

    SConscript(dirs=build,name='SConscript',exports=frame_exports)
    Default(build)
    
##########################################################################
# API BUILD
##########################################################################
api = env.get('API',[])
if type(api) is str:
    api = [api]
api.insert(0,'c')

Default('build/include')
Default('build/lib')
for dir in map(lambda x: os.path.join('api',x), api):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    api_exports = 'env root libdir '
    if dir == 'api/python':
        api_exports += 'pkgdir'
    else:
        api_exports += 'incdir'
        
    SConscript(dirs=build,name='SConstruct',exports=api_exports)
    Default(build)

##########################################################################
# SYSTEM BUILD
##########################################################################
for dir in map(lambda x: os.path.join('system',x), system):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',exports='env root bindir pkgdir')
    Default(build)

##########################################################################
# USER BUILD
##########################################################################
for dir in map(lambda x: os.path.join('user',x), user):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    SConscript(dirs=build,name='SConstruct',
        exports='env root bindir pkgdir')
    Default(build)

##########################################################################
# PLOT BUILD
##########################################################################
pdirs = ('lib','main','test','plplot')

for dir in map(lambda x: os.path.join('plot',x), pdirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)

    if dir in ('plot/main','plot/test'):
        plot_exports = 'env root bindir pkgdir'
    elif dir == 'plot/lib':
        plot_exports = 'env root libdir incdir pkgdir'
    elif dir == 'plot/plplot':
        plot_exports = 'env root libdir incdir bindir pkgdir'

    SConscript(dirs=build,name='SConstruct', exports=plot_exports)
    Default(build)

##########################################################################
# PENS BUILD
##########################################################################
pdirs = ('fonts','include','utilities','genlib','main','docs','scripts')

for dir in map(lambda x: os.path.join('pens',x), pdirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    if dir == 'pens/main':
        pens_exports = 'env root pkgdir bindir'
        sconscript = 'SConstruct'
    elif dir == 'pens/scripts':
        pens_exports = 'env bindir pkgdir'
        sconscript = 'SConscript'
    else:
        pens_exports = 'env root incdir libdir bindir'
        sconscript = 'SConstruct'
    SConscript(dirs=build,name=sconscript,exports=pens_exports)
    Default(build)

##########################################################################
# SU BUILD
##########################################################################
sudirs = ('lib','main','plot')

for dir in map(lambda x: os.path.join('su',x), sudirs):
    build = os.path.join('build',dir)
    BuildDir(build,dir)
    if dir in ('su/main','su/plot'):
        su_exports = 'env root pkgdir bindir'
    else:
        su_exports = 'env root libdir bindir incdir'
    SConscript(dirs=build,name='SConstruct',
               exports=su_exports)
    Default(build)

##########################################################################
# INSTALLATION
##########################################################################

env.Alias('install',[incdir, bindir, pkgdir, 
                     libdir, docdir, spcdir, mandir, etcdir])

