import os

root = os.environ.setdefault('RSFROOT',os.getcwd())

bindir = os.path.join(root,'bin')
libdir = os.path.join(root,'lib')
incdir = os.path.join(root,'include')

env = Environment(CCFLAGS='-std=gnu9x -Wall -pedantic -g',
                  CPPPATH=incdir,
                  LIBPATH=libdir,
                  LIBS=['rsfplot','rsf','m'],
                  PROGPREFIX='sf')

for src in ('doc','proj','prog'):
    py = "rsf%s.py"% src
    env.Install(libdir,py)
    Clean(os.path.join(libdir,py),os.path.join(libdir,py+'c'))
env.Install(bindir,'sfdoc')

Export('env')
SConscript(dirs=['seis/rsf','vplot/lib'],
           name='SConstruct')
SConscript(dirs=['seis/main','seis/proc','seis/imag','vplot/main'],
           name='SConstruct')

env.Alias('all',[bindir,libdir,incdir])


