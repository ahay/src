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

from os import environ
from SCons.Script import AddOption
from sfconfigure.stdtools import load_options
"""
Doc of RSF tool
"""

from sfconfigure.custom_builders.AliasIfExists import AliasIfExists
from sfconfigure.custom_builders.Install import InstallPythonPackage,\
    InstallTool, InstallPythonExecutable
from sfconfigure.custom_builders.Install import InstallPythonModule
from sfconfigure.custom_builders.Install import InstallLibrary
from sfconfigure.custom_builders.Install import InstallExecutable
from sfconfigure.custom_builders.Install import InstallInclude

#from sfconfigure.custom_builders.MainDoc import SelfDocBuilder
#from sfconfigure.custom_builders.MainDoc import SelfDocEmitter
#from sfconfigure.custom_builders.MainDoc import BuildMainsDoc
#from sfconfigure.custom_builders.MainDoc import BuildMainsDocEmitter

HAVE_SCONFIG_TOOL_DIR_OPTION = False
HAVE_SCONFIG_FATAL_ERROR_OPTION = False

def generate( env ):    
    """
    Doc of generate
    """
    
    global HAVE_SCONFIG_TOOL_DIR_OPTION
    global HAVE_SCONFIG_FATAL_ERROR_OPTION
    
#    DocMain = env.Builder( action=SelfDocBuilder,
#                          emitter=SelfDocEmitter,
#                          suffix={'.cpp':'.hpp',
#                                  '.c':'.h' } 
#                          )
#    
#    DocMains = env.Builder( action=BuildMainsDoc,
#                          emitter=BuildMainsDocEmitter,
#                          )
    
    
    env.AddMethod(  AliasIfExists , "AliasIfExists" )
    env.AddMethod(  InstallPythonPackage , "InstallPythonPackage" )
    env.AddMethod(  InstallPythonModule , "InstallPythonModule" )
    env.AddMethod(  InstallLibrary , "InstallLibrary" )
    env.AddMethod(  InstallInclude , "InstallInclude" )
    env.AddMethod(  InstallExecutable , "InstallExecutable" )
    env.AddMethod(  InstallTool , "InstallTool" )
    env.AddMethod(  InstallPythonExecutable , "InstallPythonExecutable" )
    
#    env["BUILDERS"]["DocMain"] = DocMain
#    env["BUILDERS"]["BuildMainsDoc"] = DocMains
    
    
    load_options( )
    
    SCONFIG_TOOL_DIR = env.GetOption('tool_prefix')
    
    if not SCONFIG_TOOL_DIR:
        SCONFIG_TOOL_DIR = environ.get( 'SCONFIG_TOOLS' , None )
    
    if SCONFIG_TOOL_DIR:
        env['tool_prefix'] = SCONFIG_TOOL_DIR

    
    
    return


    
def exists( env ):
    """
    Doc of exist
    """
    return True
    

