
from os.path import join
def generate( env ):
    clone = env.Clone()
    clone.Tool('rsfroot')
    
    root = clone['RSFSRC']
    env.Tool('rsfcc_debug')
    
    env['lib_prefix'] = join( root, 'build/debug/lib' )
    env['inlude_prefix'] = join( root, 'build/debug/include' )
    
    
    if clone['inlude_prefix'] in env['CPPPATH']:
        env['CPPPATH'].remove( clone['inlude_prefix'] )
        
    if clone['lib_prefix'] in env['LIBPATH']:
        env['LIBPATH'].remove( clone['lib_prefix'] )
        
    env.Append( CPPPATH=env['inlude_prefix'] )
    env.Append( LIBPATH=env['lib_prefix'] )
    
def exists( env ):
    return 1
