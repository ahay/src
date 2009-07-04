from rsfproj import *

def modmig(par):
    Flow('layermodel',None,'vam nz=%(nz)d nx=%(nx)d dz=%(dz)g dx=%(dx)g' % par)
    Result('layermodel','grey color=j allpos=y bias=1.5 scalebar=y wanttitle=n')

    Flow('trace wave','layermodel',
         '''
         wave nt=%(nt)d dt=%(dt)g isz=%(isz)d isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d
         igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d fpeak=%(fpeak)g source=${TARGETS[1]} nm=%(nm)d
         ''' % par)
    Plot('wave','grey')
    Result('trace','grey title=Data')
