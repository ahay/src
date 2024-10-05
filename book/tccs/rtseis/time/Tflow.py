'extra flows'

import string, re, sys, os
from rsf.proj import *


def Tflow(target, source, command, prefix='sf',):

    if sys.platform == 'darwin':
        time_nm='gtime'
    else:
        time_nm='time'

    timer=WhereIs(time_nm)
    if timer==None:
        sys.stderr.write('Tflow need %s.'%time_nm)
        sys.exit(1)

    if isinstance(target, list):
        tfiles = target
    else:
        tfiles = target.split(' ')
#	tfiles.insert(0, tfiles[0]+'_runtime')
    pars=command.split()
    p0=pars.pop(0)

    p1=WhereIs(p0) # sfdip
    if p1==None :
        p1=WhereIs(prefix+p0) # dip
    if re.match(r'[^/]+\.exe$',p0) :
        p1=os.path.join('.',p0) # Mdip.exe

    pars.insert(0,p1)
    cmd=' '.join(pars)
    Flow(tfiles, source,
            '''
            ( %s -f "%%S %%U" %s <${SOURCES[0]} >/dev/null 2> time.out ) &&
            (tail -1 time.out;
            echo in=time0.rsf n1=2 data_format=ascii_float)
                    > time0.rsf &&
            dd form=native < time0.rsf | stack axis=1 norm=n
            > time1.rsf &&
            cmplx ${SOURCES[2]} time1.rsf > ${TARGETS[0]} &&
            %s -f time.out time0.rsf && rm time1.rsf
            '''%(timer, cmd, WhereIs('rm')),
            stdin=None, stdout=None)
            
    #Flow(tfiles, source,
    #        '''
    #        ( %s -f "%%S %%U" %s <${SOURCES[0]} >/dev/null )
    #                    >& time.out &&
    #        (tail -1 time.out;
    #        echo in=time0.rsf n1=2 data_format=ascii_float)
    #                > time0.rsf &&
    #        dd form=native < time0.rsf | stack axis=1 norm=n
    #        > time1.rsf &&
    #        cmplx ${SOURCES[2]} time1.rsf > ${TARGETS[0]} &&
    #        %s -f time.out time0.rsf && rm time1.rsf
    #        '''%(timer, cmd, WhereIs('rm')),
    #        stdin=None, stdout=None)
