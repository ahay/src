
from os.path import join
def generate( env ):
    clone = env.Clone()
    clone.Tool('rsfroot')
    
    root = clone['RSFSRC']
    env.Tool('rsfcc_debug')
    env.Tool('rsfc')
    
    env['lib_prefix'] = join( root, 'build/debug/lib' )
    env['include_prefix'] = join( root, 'build/debug/include' )
    
    
    if clone['include_prefix'] in env['CPPPATH']:
        env['CPPPATH'].remove( clone['include_prefix'] )
        
    if clone['lib_prefix'] in env['LIBPATH']:
        env['LIBPATH'].remove( clone['lib_prefix'] )
        
    env.Append( CPPPATH=env['include_prefix'] )
    env.Append( LIBPATH=env['lib_prefix'] )
    
def exists( env ):
    return 1
