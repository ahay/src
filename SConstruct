#EnsureSConsVersion(1,0)

import atexit, os, sys
sys.path.insert(0,'./framework')
import bldutil, configure, setenv, rsf.doc

env = Environment()

if os.path.isfile('config.py'):
    import config
    root = config.RSFROOT
else:
    root = os.environ.get('RSFROOT')

if not root:
    root = sys.prefix
    print 'Setting RSFROOT to "%s" ' % root

srcdir = os.getcwd()
bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')
shrdir = os.path.join(root,'share')
pkgdir = setenv.get_pkgdir(root)

etcdir = os.path.join(shrdir, 'madagascar', 'etc')

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
              '#/framework/rsf/__init__.pyc',
              '#/framework/rsf/doc.pyc',
              '#/framework/rsf/path.pyc'])
env.Alias('config',config)
env.Install(etcdir,'config.py') 

env.InstallAs('#/build/framework/rsf/conf.py','framework/configure.py') 

# ----------- Environment variable setup scripts -----------

etcdir2 = os.path.join(root, 'etc', 'madagascar')

for sh in ('sh','csh'):
    shrc  = env.Command('env.'+sh,  '', setenv.shell_script,
                        varlist=('shell'),shell=sh)
    env.Alias('config',shrc)
    env.Install(etcdir,shrc)

    # backward compatibility
    if os.path.isdir(etcdir2):
        env.Install(etcdir2,shrc)

env['version'] = bldutil.__read_version_file('VERSION.txt')

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
if os.path.isdir('user'):
    user = filter(lambda x: x[0] != '.' and x != 'nobody', os.listdir('user'))
    # Avoid crashing when user places some files in RSFSRC/user
    user = filter(lambda x: os.path.isdir(os.path.join('user',x)), user)
    # Avoid crashing when user places directories in RSFSRC/user
    user = filter(lambda x: os.path.isfile(os.path.join('user',x,'SConstruct')), 
        user)
else:
    user = []

trip = []
if os.path.isdir('trip'):
    tripdirs = ('iwave','rvl','iwave++')
    for dir in map(lambda x: os.path.join('trip',x), tripdirs):
        subpath = os.path.join(dir,'hsubpath')
        if os.path.exists(subpath):
            f = open(subpath,'r')
            dirs = (f.read().strip('\n')).split(':')
            f.close() 

            dirs = filter(lambda x: os.path.isdir(os.path.join(dir,x,'main')), dirs)
            trip.extend(map(lambda x: os.path.join(dir,x), dirs))

dotproj = Glob('book/[a-z]*/[a-z]*/[a-z]*/.rsfproj')

frame_exports = 'env bindir libdir pkgdir shrdir srcdir system user trip dotproj'

for dir in map(lambda x: os.path.join('framework',x),Split('rsf doc ptools')):
    build = os.path.join('build',dir)
    if configure.version[0] > 1:
        VariantDir(build,dir)
    else:
        BuildDir(build,dir)
    SConscript(dirs=build,name='SConscript',exports=frame_exports)
    Default(build)

for dir in map(lambda x: os.path.join('book',x),Split('Recipes gallery')):
    build = os.path.join('build',dir)
    if configure.version[0] > 1:
        VariantDir(build,dir)
    else:
        BuildDir(build,dir)
    SConscript(dirs=build,name='SConscript',exports='env srcdir pkgdir')
    Default(build)
    
##########################################################################
# API BUILD
##########################################################################
api = env.get('API',[])
if type(api) is str:
    api = [api]
api.insert(0,'c')
api.insert(1,'python')

Default('build/include')
Default('build/lib')
for dir in map(lambda x: os.path.join('api',x), api):
    build = os.path.join('build',dir)
    if configure.version[0] > 1:
        VariantDir(build,dir)
    else:
        BuildDir(build,dir)
    api_exports = 'env root libdir incdir'
    if dir == 'api/python':
        api_exports += ' pkgdir bindir'
        
    SConscript(dirs=build,name='SConstruct',exports=api_exports)
    Default(build)

##########################################################################
# SYSTEM BUILD
##########################################################################

for dir in map(lambda x: os.path.join('system',x), system):
    build = os.path.join('build',dir)
    if configure.version[0] > 1:
        VariantDir(build,dir)
    else:
        BuildDir(build,dir)
        
    SConscript(dirs=build,name='SConstruct',exports='env root bindir pkgdir libdir incdir')
    Default(build)

##########################################################################
# USER BUILD
##########################################################################

if os.path.isdir('user'):
    for dir in map(lambda x: os.path.join('user',x), user):
        if os.path.isfile(os.path.join(dir, 'SConstruct')):
            build = os.path.join('build',dir)
            if configure.version[0] > 1:
                VariantDir(build,dir)
            else:
                BuildDir(build,dir)
                
            SConscript(dirs=build,name='SConstruct', 
                       exports='env root bindir pkgdir libdir incdir')
            Default(build)

##########################################################################
# PLOT BUILD
##########################################################################

if os.path.isdir('plot'):
    pdirs = ('lib','main','test','plplot')
    for dir in map(lambda x: os.path.join('plot',x), pdirs):
        build = os.path.join('build',dir)
        if configure.version[0] > 1:
            VariantDir(build,dir)
        else:
            BuildDir(build,dir)
        if dir in ('plot/main','plot/test'):
            plot_exports = 'env root bindir libdir pkgdir'
        elif dir == 'plot/lib':
            plot_exports = 'env root libdir incdir pkgdir'
        elif dir == 'plot/plplot':
            plot_exports = 'env root libdir incdir bindir pkgdir'
        SConscript(dirs=build,name='SConstruct', exports=plot_exports)
        Default(build)

##########################################################################
# PENS BUILD
##########################################################################

if os.path.isdir('pens'):
    pdirs = ('fonts','include','utilities','genlib','main','docs','scripts')
    for dir in map(lambda x: os.path.join('pens',x), pdirs):
        build = os.path.join('build',dir)
        if configure.version[0] > 1:
            VariantDir(build,dir)
        else:
            BuildDir(build,dir)
        if dir == 'pens/main':
            pens_exports = 'env root pkgdir libdir bindir'
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

if os.path.isdir('su'):
    sudirs = ('lib','main','plot')
    for dir in map(lambda x: os.path.join('su',x), sudirs):
        build = os.path.join('build',dir)
        if configure.version[0] > 1:
            VariantDir(build,dir)
        else:
            BuildDir(build,dir)
        if dir in ('su/main','su/plot'):
            su_exports = 'env root pkgdir bindir'
        else:
            su_exports = 'env root libdir bindir incdir'
        SConscript(dirs=build,name='SConstruct',
                   exports=su_exports)
        Default(build)

##########################################################################
# TRIP BUILD
##########################################################################

if os.path.isdir('trip'):

    tripdirs = ('iwave','rvl','iwave++')
    for dir in map(lambda x: os.path.join('trip',x), tripdirs):
        try:
            subpathfile = os.path.join(dir,'hsubpath')
            f = open(subpathfile,'r')
            sublist = (f.read().strip('\n')).split(':')
            f.close()
        except:
            sublist = []
        for sub in sublist:
            subdir = os.path.join(dir,sub)
            build = os.path.join('build',subdir)
            if configure.version[0] > 1:
                VariantDir(build,subdir)
            else:
                BuildDir(build,subdir)
            trip_exports = 'env root libdir bindir incdir pkgdir'
            SConscript(dirs=build,exports=trip_exports)
            Default(build)

##########################################################################
# INSTALLATION
##########################################################################

docdir = os.path.join(shrdir, 'doc', 'madagascar') 	 
for docfile in Split('AUTHORS COPYING NEWS README'): 	 
    env.Install(docdir,docfile+'.txt')

env.Alias('install',[incdir, bindir, pkgdir, libdir, shrdir, etcdir])

# backward compatibility
if os.path.isdir(etcdir2):
    env.Alias('install',etcdir2)

# end-of-installation message
def msgEndInstall():
    from SCons.Script import GetBuildFailures
    if not GetBuildFailures():
        print '''
---------------------------------------------------------
To start using madagascar, source env.sh or env.csh from:
    %s/
Local documentation center at:
    %s/
Documentation wiki at http://www.ahay.org
---------------------------------------------------------
''' % (etcdir, docdir)

if 'install' in COMMAND_LINE_TARGETS:
    atexit.register(msgEndInstall)

##########################################################################
# PYTHON PACKAGE DOCUMENTATION
##########################################################################

epydoc = WhereIs('epydoc')
if epydoc:
    epydir = os.path.join(docdir,'epydoc')
    envcmd = 'PYTHONPATH=%s %s' % (setenv.get_local_site_pkgs(root),epydoc)   
    epyargs = '--exclude rsf.use --exclude rsf.vplot --exclude rsf.sf* '
    epyargs += '--html -qqq --no-private --graph classtree rsf'
    env.Command(epydir, pkgdir, envcmd + ' -o $TARGET ' + epyargs)
