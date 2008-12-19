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
"""
"""

from os.path import split,join,isfile
from os.path import isdir,abspath

import os
from SCons.Script import AddOption,SConscript,Environment,Dir
from glob import glob
import distutils.sysconfig

from sfconfigure import DoConfig
from sfconfigure import add_to_default_toolpath
from sfconfigure.autotoolgen import ToolCreator
#from sfconfigure.sort_components import ComponentsCmp

HAVE_SCONFIG_TOOL_DIR_OPTION = False
HAVE_SCONFIG_FATAL_ERROR_OPTION = False

def load_options( ):
    global HAVE_SCONFIG_TOOL_DIR_OPTION
    global HAVE_SCONFIG_FATAL_ERROR_OPTION

    if not HAVE_SCONFIG_TOOL_DIR_OPTION:
        AddOption('--tool-prefix',
          dest='tool_prefix',
          nargs=1, type='string',
          action='store',
          metavar='DIR',
          help='tool installation prefix')
        HAVE_SCONFIG_TOOL_DIR_OPTION = True
    
    if not HAVE_SCONFIG_FATAL_ERROR_OPTION:
        HAVE_SCONFIG_FATAL_ERROR_OPTION = True
        AddOption( '--fatal-error',
           dest='fatal',
           action='store_true',
           default=False,
           help='do not ignore errors' )

def do_missing_tool( env, ex_tool ):
    tool_name = split(ex_tool)[-1]
    if not (env.GetOption('clean') or env.GetOption('help') ):
        print  "    | Tool: '%s' (disabled)" %(tool_name)
    default_py = join( ex_tool,"default")
    default_py_d = {}
    tools_to_create = []
    if isfile(default_py):
        try:
            exec open(default_py) in default_py_d
            tools_to_create.extend( default_py_d.get('tools') )
        except:
            pass 
    if not tools_to_create:
        tools_to_create.append( tool_name )
        
    for tool_to_create in tools_to_create:
        print "Writing empty tool template for '%s', exists ... no" %(tool_to_create)
        tc = ToolCreator( tool_to_create, ex_tool )
        tc.CreateTool(env)



HAVE_CONFIGFILE_OPT = False
CONFIGFILE_OPT = {}

def StdInstall( env, toolname, dest ):
    
    PREFIX = distutils.sysconfig.PREFIX
    root = os.environ.get('RSFROOT',PREFIX)
    
    
    AddOption( '--prefix',
           dest='prefix',
           nargs=1, type='string',
           action='store',
           metavar='DIR',
           default=root,
           help='installation prefix' )

    AddOption( '--lib-prefix',
           dest='lib_prefix',
           nargs=1, type='string',
           action='store',
           metavar='DIR',
           help='installation prefix' )

    AddOption( '--bin-prefix',
           dest='bin_prefix',
           nargs=1, type='string',
           action='store',
           metavar='DIR',
           help='installation prefix' )

    AddOption( '--python-prefix',
           dest='python_prefix',
           nargs=1, type='string',
           action='store',
           metavar='DIR',
           help='installation prefix' )

    AddOption( '--include-prefix',
           dest='include_prefix',
           nargs=1, type='string',
           action='store',
           metavar='DIR',
           help='installation prefix' )

    AddOption( '--no-std-python-install',
           dest='std_python',
           action='store_false',
           metavar='DIR',
           default=True,
           help='install python modules to /prefix/lib instead of standard location' )
    
    tc = ToolCreator( toolname, dest )
    
    
    INSTALL_PREFIX = abspath(env.GetOption('prefix') or root )
    lib_prefix = env.GetOption('lib_prefix') or join(INSTALL_PREFIX ,'lib' )
    bin_prefix = env.GetOption('bin_prefix') or join(INSTALL_PREFIX ,'bin')

    include_prefix = env.GetOption('include_prefix') or join(INSTALL_PREFIX, 'include' )
    
    
    std_python = env.GetOption('std_python')

    if not std_python:
        python_prefix = join(INSTALL_PREFIX, 'lib' )
    else:
        python_prefix = distutils.sysconfig.get_python_lib( prefix=INSTALL_PREFIX )
        
    
    python_prefix = env.GetOption('python_prefix') or python_prefix
    
    tc.Replace( INSTALL_PREFIX=INSTALL_PREFIX,
                lib_prefix=lib_prefix,
                bin_prefix=bin_prefix,
                python_prefix=python_prefix,
                include_prefix=include_prefix,
                tool_dest = dest,
                RSFSRC=Dir('#').abspath
                
                )
    tc.Exists(True)
    
    tc.CreateTool( env )
#    env.Install( dest, "%s.py" %toolname )

HAVE_TOOLPATH_OPT = False

def CommandLineToolPath( env ):
    global HAVE_TOOLPATH_OPT
    
    if not HAVE_TOOLPATH_OPT:
        AddOption( '--local-toolpath',
               dest='local_toolpath',
               nargs=1, type='string',
               action='store',
               metavar='DIR',
               help='configuration defaults' )
        
        HAVE_TOOLPATH_OPT = True
   
    
    local_toolpath = env.GetOption( 'local_toolpath' )
    
    if local_toolpath:
        local_toolpath = local_toolpath.split(":")
    else:
        local_toolpath = []
        
    for toolpath in local_toolpath:
        add_to_default_toolpath(toolpath)
    
    return
    
    

def LoadConfigFile( nenv ):
    
    global HAVE_CONFIGFILE_OPT
    if not HAVE_CONFIGFILE_OPT:
        
        AddOption( '--config-file',
               dest='config_file',
               nargs=1, type='string',
               action='store',
               metavar='DIR',
               default=None,
               help='configuration defaults' )
        
        HAVE_CONFIGFILE_OPT = True
    
        config_file = nenv.GetOption( 'config_file' )
        
        if config_file:
            config_file = abspath(config_file)
            exec open(config_file) in CONFIGFILE_OPT
    
    nenv.Replace( **CONFIGFILE_OPT )
    
    return


#===============================================================================
# Local variables and functions
#===============================================================================

_sconscript = lambda path: isfile(path) and SConscript( path )

def add_list(env,name):
    
    disable_ = env.GetOption(name)
    if not disable_:
        disable_ = []
    else:
        disable_ = disable_.split(",")
        
    return disable_

#def AddSLABOptions( env=None ):
#    
#    AddOption( '--fatal-error',
#           dest='fatal',
#           action='store_true',
#           default=False,
#           help='do not ignore errors' )
#
#    AddOption( '--disable-users',
#           dest='disable_users',
#           action='store',
#           default="",
#           help='do not configure users folders' )
#    
#    AddOption( '--disable-components',
#           dest='disable_components',
#           action='store',
#           default="",
#           help='do not configure or run components' )
#    
#    AddOption( '--disable-tools',
#           dest='disable_tools',
#           action='store',
#           default="",
#           help='do not configure specified external tools' )
#    
#    if env:
#        
#        env['fatal_error'] = env.GetOption('fatal') or False
#                    
#        env['disable_users'] = add_list(env,'disable_users')
#        env['disable_components'] = add_list(env,'disable_components')
#        env['disable_tools'] = add_list(env,'disable_tools')
#
#    return 


def standard_tools_config( env=None, toolsdir='external_tools' ):
    """
    """
    
    if env is None:
        env = Environment( )

    unused_external_tools = [ cnt for cnt in glob( join( toolsdir,"_*" ) ) if isdir(cnt) ]
    external_tools = [ cnt for cnt in glob( join( toolsdir,"*" ) ) if isdir(cnt) and cnt not in unused_external_tools]
    
    if DoConfig( env ):
        print "    +-------------------------------+"
        print "    |        Configure Tools        |"
        print "    +-------------------------------+"
    
    ccmp = ComponentsCmp( external_tools, None )
    external_tools.sort( ccmp )
    external_tools.sort( ccmp )
    for ex_tool in external_tools:
        tool_name = split(ex_tool)[-1]
        if tool_name in env.get('disable_tools',[]):
            do_missing_tool( env, ex_tool )
        else:
            if not (env.GetOption('clean') or env.GetOption('help') ):
#                print  "+ Configuring external tool:", ex_tool
                print  "    | Tool: '%s'" %(tool_name)
            try:
                _sconscript( join(ex_tool,"SConfig") )
                add_to_default_toolpath( abspath( ex_tool ) ) 
            except Exception, e:
                if env.GetOption('fatal'):
                    raise
                print "***********************************************************"
                print "SConfigure ERROR: Not including component '%(ex_tool)s': got exception in %(ex_tool)s" %vars()
                print "Exception: %(e)s" %vars()
                print "This may make the slab package unstable"
                print "     To be safe use the option 'disable_component=%(ex_tool)s'" %vars()
                print "     also you can use the '--fatal-error' command line option to "
                print "     see the full error message"
                print "***********************************************************"
                do_missing_tool( env, ex_tool )
    return 
 
 