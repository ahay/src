from rsfproj import *

def modmig(par):
    # Layered model
    Flow('layermodel',None,'vam nz=%(nz)d nx=%(nx)d dz=%(dz)g dx=%(dx)g' % par)
    Result('layermodel',
           'grey color=j allpos=y bias=1.5 scalebar=y wanttitle=n barreverse=y')

    wave = '''
    wave nt=%(nt)d dt=%(dt)g isz=%(isz)d 
    isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d
    igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d 
    fpeak=%(fpeak)g source=${TARGETS[1]} nm=%(nm)d
    ''' % par

    # Variable velocity
    Flow('trace wave','layermodel',wave)
    Plot('wave','grey gainpanel=all',view=1)
    Result('trace','grey title=Data')

    # Constant velocity
    Flow('hom','layermodel','math output=1.5')
    Flow('tracehom wavehom','hom',wave)
    Plot('wavehom','grey gainpanel=all',view=1)
    Result('tracehom','grey title="Data (Homogeneous Model)"')

    # Difference
    Flow('tracediff','trace tracehom','add scale=1,-1 ${SOURCES[1]}')
    Result('tracediff','grey title="Data (Reflection)"')

    
