from rsfproj import *

def Helderiv(name,eps=0.001,na=16):
    '''Creates helix derivative filter'''

    def tag(rsf):
        return '-'.join([name,rsf])
    
    Flow(tag('slag0'),None,
         'echo 1 1000 n=1000,1000 n1=2 in=$TARGET data_format=ascii_int')
    Flow(tag('slag'),tag('slag0'),'dd data_format=native_int')
    Flow(tag('ss0'),tag('slag'),
         '''echo -1 -1 a0=%g n1=2
         lag=$SOURCE in=$TARGET data_format=ascii_float''' % (2.0+0.5*eps))
    Flow(tag('ss'),tag('ss0'),'dd data_format=native_float')
    
    Flow(tag('alag0'),None,
         'echo %s n=1000,1000 n1=%d in=$TARGET data_format=ascii_int' %
         (' '.join(map(str, range(1,na+1) + range(1001-na,1001))),2*na))
    Flow(tag('alag'),tag('alag0'),'dd data_format=native_int')

    Flow([name,tag('lag')],[tag('ss'),tag('alag0')],
         'wilson lagin=${SOURCES[1]} lagout=${TARGETS[1]}')

