from string import split

def Helderiv(out,inp,eps=0.001,na=16,cut=0):
    '''Applies helix derivative to create out from inp'''

    def uniq(name):
        return out+name

    def F(dst,src,flow,**kw):
        if src:
            return apply(Flow,(map(uniq,split(dst)),
                               map(uniq,split(src)),flow),kw)
        else:
            return apply(Flow,(map(uniq,split(dst)),None,flow),kw)
    
    F('slag0',None,
      'echo 1 1000 n=1000,1000 n1=2 in=$TARGET data_format=ascii_int')
    F('slag','slag0','dd data_format=native_int')
    F('ss0','slag',
      '''echo -1 -1 a0=%g n1=2
      lag=$SOURCE in=$TARGET data_format=ascii_float''' % (2.0+0.5*eps))
    F('ss','ss0','dd data_format=native_float')
    
    F('alag0',None,
      'echo %s n=1000,1000 n1=%d in=$TARGET data_format=ascii_int' %
      (' '.join(map(str, range(1,na+1) + range(1001-na,1001))),2*na))
    F('alag','alag0','dd data_format=native_int')

    F('aa lag','ss alag0',
      'wilson lagin=${SOURCES[1]} lagout=${TARGETS[1]}')

    if cut:
        F('tt0','slag',
          '''echo -1 -1 a0=2.0005 n1=2
          lag=$SOURCE in=$TARGET data_format=ascii_float''')
        F('tt','tt0','dd data_format=native_float')
        
        F('bb blag','tt alag0',
          'wilson lagin=${SOURCES[1]} lagout=${TARGETS[1]}')

        Flow(out,[inp,uniq('aa'),uniq('bb')],
             '''
             helicon filt=${SOURCES[1]} div=1 |
             helicon filt=${SOURCES[2]} div=0
             ''')
    else:
        Flow(out,[inp,uniq('aa')],'helicon filt=${SOURCES[1]}')
