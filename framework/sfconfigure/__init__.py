#  Copyright (C) 2008 The University of British Columbia
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  Written By Sean Ross-Ross
#

from SCons.Script import EnsureSConsVersion as _EnsureSConsVersion 
from SCons.Script import Mkdir as _MkDir, Tool as _scons_tool
from os.path import dirname as _dname, join as _join, isdir as _isdir
from sfconfigure.project import add_to_default_toolpath
from sfconfigure.project import add_to_default_tools
from os import environ
import check_requirements

from autotoolgen import ToolCreator

def CreateOptionsFileName( env , name ):
    libpath =  env.Dir('#/config-etc/variables').abspath
    if not _isdir(libpath): env.Execute( _MkDir( libpath ) )
    
    vpath = _join( libpath, name ) 
    env.Clean( '.', vpath )
    env.Clean( '.', "#/config-etc" )
    env.Clean( '.', "#/.sconsign.dblite" )
    return vpath

def DoConfig( env ):
    """
    Only do config stuff if not cleaning or help or dryrun
    """
    return not (env.GetOption('clean') or env.GetOption('help') or env.GetOption('no_exec') )

    
_EnsureSConsVersion( 0, 97 )


_this_path = _dname(__file__)
sconfig_toolpath = _join( _this_path, 'tools' )


add_to_default_toolpath( sconfig_toolpath )

SCONFIG_TOOLS = environ.get( "SCONFIG_TOOLS" ) 
if SCONFIG_TOOLS:
    tooldirs = SCONFIG_TOOLS.split( ":" )
    add_to_default_toolpath( tooldirs )
    pass

add_to_default_tools( 'sconf' )


from stdtools import standard_tools_config
from stdtools import LoadConfigFile
from stdtools import StdInstall

from config_checks import options
from config_checks import check_all



