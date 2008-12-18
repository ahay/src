
def generate( env ):
    
    clone = env.Clone()
    clone.Tool('rsfroot')
    
    env['RSFROOT'] = clone['INSTALL_PREFIX']
    env.Append( CPPPATH=clone['include_prefix'],
                LIBPATH=clone['lib_prefix'],
                LIBS=['rsf'] )
    
    return


def exists( env ):
    return 1
