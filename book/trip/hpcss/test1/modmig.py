from rsf.proj import *

def Movie(title):
    Plot(title,'grey gainpanel=all wanttitle=n',view=1)
    Alias('movies',title+'.vpl')

def modmiggeom(par):
    # Layered model
    Flow('layermodel',None,'vam nz=%(nz)d nx=%(nx)d dz=%(dz)g dx=%(dx)g' % par)
    Result('layermodel',
           '''
           grey color=j allpos=y bias=1.5 
           scalebar=y wanttitle=n barreverse=y
           ''')

    wave = '''
    nt=%(nt)d dt=%(dt)g isz=%(isz)d 
    isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d
    igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d 
    fpeak=%(fpeak)g source=${TARGETS[1]} nm=%(nm)d
    ''' % par

    # Variable velocity
    Flow('trace wave','layermodel','wavegeom ' + wave)
    if par['nm'] != 0:
        Movie('wave')
        Result('wave','grey gainpanel=all wanttitle=n')
        
    Result('trace','grey title=Data')

    # Constant velocity
    Flow('hom','layermodel','math output=1.5')
    Flow('tracehom wavehom','hom','wavegeom ' + wave)
    if par['nm'] != 0:
        Movie('wavehom')
    Result('tracehom','grey title="Data (Homogeneous Model)"')

    # Difference
    Flow('tracediff','trace tracehom','add scale=1,-1 ${SOURCES[1]}')
    Result('tracediff','grey title="Data (Reflection)"')
    Plot('trace','grey title="trace"')
    Plot('tracehom','grey title="trace hom"')
    Plot('tracediff','grey title="diff"')
    Result('traces','trace tracehom tracediff','SideBySideAniso')

    # Reverse-time migration
    rtm = '''
    rtmgeom trace=${SOURCES[1]} 
    isz=%(isz)d isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d 
    igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d --readwrite=y
    ''' % par

    # Flow('source','hom tracediff',rtm + ' source=$TARGET',stdout=0)
    Flow('source','hom tracediff',rtm + ' source=$TARGET',stdout=0)
    
    #Flow('rtm rtmmov','hom tracediff source',
    #     rtm + '  source=${SOURCES[2]} receiver=${TARGETS[1]}')
    Flow('rtm rtmmov','hom tracediff source',
         rtm + '  source=${SOURCES[2]} receiver=${TARGETS[1]}')
    
    Movie('rtmmov')
    Result('rtm','grey title=Image')

    # higher accuracy
    Flow('tracehom24 wavehom24','hom','wave24 ' + wave)
    if par['nm'] != 0:
        Movie('wavehom24')


def modmig(par):
    # Layered model
    Flow('layermodel',None,'vam nz=%(nz)d nx=%(nx)d dz=%(dz)g dx=%(dx)g' % par)
    Result('layermodel',
           '''
           grey color=j allpos=y bias=1.5 
           scalebar=y wanttitle=n barreverse=y
           ''')

    wave = '''
    nt=%(nt)d dt=%(dt)g isz=%(isz)d 
    isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d
    igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d 
    fpeak=%(fpeak)g source=${TARGETS[1]} nm=%(nm)d
    ''' % par

    # Variable velocity
    Flow('trace wave','layermodel','wave ' + wave)
    if par['nm'] != 0:
        Movie('wave')
        Result('wave','grey gainpanel=all wanttitle=n')
        
    Result('trace','grey title=Data')

    # Constant velocity
    Flow('hom','layermodel','math output=1.5')
    Flow('tracehom wavehom','hom','wave ' + wave)
    if par['nm'] != 0:
        Movie('wavehom')
    Result('tracehom','grey title="Data (Homogeneous Model)"')

    # Difference
    Flow('tracediff','trace tracehom','add scale=1,-1 ${SOURCES[1]}')
    Result('tracediff','grey title="Data (Reflection)"')
    Plot('trace','grey title="trace"')
    Plot('tracehom','grey title="trace hom"')
    Plot('tracediff','grey title="diff"')
    Result('traces','trace tracehom tracediff','SideBySideAniso')

    # Reverse-time migration
    rtm = '''
    rtm trace=${SOURCES[1]} 
    isz=%(isz)d isxbeg=%(isxb)d isxend=%(isxe)d iskip=%(isxd)d 
    igz=%(igz)d igxbeg=%(igxb)d igxend=%(igxe)d --readwrite=y
    ''' % par

    # Flow('source','hom tracediff',rtm + ' source=$TARGET',stdout=0)
    Flow('source','hom tracediff',rtm + ' source=$TARGET',stdout=0)
    
    #Flow('rtm rtmmov','hom tracediff source',
    #     rtm + '  source=${SOURCES[2]} receiver=${TARGETS[1]}')
    Flow('rtm rtmmov','hom tracediff source',
         rtm + '  source=${SOURCES[2]} receiver=${TARGETS[1]}')
    
    Movie('rtmmov')
    Result('rtm','grey title=Image')

    # higher accuracy
    Flow('tracehom24 wavehom24','hom','wave24 ' + wave)
    if par['nm'] != 0:
        Movie('wavehom24')

