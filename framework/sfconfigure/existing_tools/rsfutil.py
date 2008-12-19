
from sfconfigure.custom_builders.utils import pycompile
from sfconfigure.custom_builders.utils import header
from sfconfigure.custom_builders.utils import included
from sfconfigure.custom_builders.utils import placeholder
from sfconfigure.custom_builders.utils import docmerge
from sfconfigure.custom_builders.utils import merge
from sfconfigure.custom_builders.utils import pycompile_emit

from SCons.Script import Builder, Action, Scanner
def generate( env ):    
    """
    Doc of generate
    """
    pass
    Pycompile = Builder(action=pycompile)
    
    Header = Builder (action = Action(header,varlist=['prefix']),
                  src_suffix='.c',suffix='.h')

    Include = Scanner(name='Include',function=included,skeys=['.c'])

    Place = Builder (action = Action(placeholder,varlist=['var','package']))
    Docmerge = Builder(action=Action(docmerge,varlist=['alias']),
                   emitter=pycompile_emit)
    
    Merge = Builder(action=Action(merge) )
    env.Append(BUILDERS={'Include':   Header,
                         'Place':     Place,
                         'Pycompile': Pycompile,
                         'Docmerge':  Docmerge,
                         'Merge' :Merge},
               SCANNERS=[Include])


def exists( env ):
    return 1



