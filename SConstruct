import os, sys
import configure

here = os.getcwd()
root = os.environ.setdefault('RSFROOT',here)

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')

env = Environment()

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
env.Alias('config',config)
env.InstallAs(os.path.join(libdir,'rsfconfig.py'),'config.py')
Clean(config,['#/config.log','#/.sconf_temp'])

##########################################################################
# DOCUMENTATION
##########################################################################
for src in ('doc','proj','prog'):
    py = "rsf%s.py"% src
    env.Install(libdir,py)
    Clean(os.path.join(libdir,py),os.path.join(libdir,py+'c'))
env.Install(bindir,'sfdoc')
##########################################################################

Export('env')
libdirs = ['seis/rsf','vplot/lib']
libdirs = ['seis/rsf']
progdirs = ['seis/main','seis/proc','seis/imag','vplot/main']

BuildDir('seisrsf','seis/rsf')
SConscript(dirs='seisrsf',name='SConstruct')
#SConscript(dirs=libdirs,name='SConstruct')
#SConscript(dirs=progdirs,name='SConstruct')

env.Alias('install',[bindir,libdir,incdir])
Default(libdirs+progdirs)

