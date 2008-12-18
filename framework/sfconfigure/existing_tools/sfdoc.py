

from SCons.Script import Builder,Action

import rsfdoc


def generate( env ):
    getprogs_act = Action(rsfdoc.getprogs,varlist=['dirs'])
    
    getprogs_bld = Builder( action=getprogs_act,
                            suffix='.py'
                            )
    
    Doc = Builder (action = Action(rsfdoc.selfdoc,varlist=['rsfprefix','lang']),
               src_suffix='.c',suffix='.py')
    
    
    use_action=Action(rsfdoc.use,varlist=['book'])
    
    Use = Builder( action=use_action)
    
    env.Append(BUILDERS = {'Doc' : Doc,
                           'GetProgs':  getprogs_bld,
                           'Use' : Use}
    )
    
    



def exists( env ):
    return 1
