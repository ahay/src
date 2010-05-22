import string

from rsf.proj import *

def Helderiv(name,eps=0.001,na=16):
    '''Creates helix derivative filter'''

    def tag(name,rsf):
        return '%s-%s' % (name,rsf)

    # ------------------------------------------------------------
    Flow(tag(name,'slag0.asc'),None,
         'echo 1 1000 n=1000,1000 n1=2 in=$TARGET data_format=ascii_int')
    Flow(tag(name,'slag'),tag(name,'slag0.asc'),'dd form=native')

    # ------------------------------------------------------------
    Flow(tag(name,'ss0.asc'),tag(name,'slag'),
         '''echo -1 -1 a0=%g n1=2
         lag=$SOURCE in=$TARGET data_format=ascii_float''' % (2.0+0.5*eps))
    Flow(tag(name,'ss'),tag(name,'ss0.asc'),'dd form=native')

    # ------------------------------------------------------------    
    Flow(tag(name,'alag0.asc'),None,
         'echo %s n=1000,1000 n1=%d in=$TARGET data_format=ascii_int' %
         (string.join(map(str, range(1,na+1) + range(1001-na,1001)),' '),2*na))
    Flow(tag(name,'alag'),tag(name,'alag0.asc'),'dd form=native')

    # ------------------------------------------------------------
    Flow([name,tag(name,'lag')],[tag(name,'ss'),tag(name,'alag0.asc')],
         'wilson lagin=${SOURCES[1]} lagout=${TARGETS[1]}')

