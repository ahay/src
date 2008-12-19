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
from SCons.Script import Tool as _Tool

def check( mains ):
    """
    Check if SCons nodes have any dependencies on tools that 
    do not exists.
    """
    failed_tools = set( )
    
    if not mains:
        return failed_tools
    
    for source in mains:
        if not (hasattr(source,'env') and source.env):
            continue
        
        failed_tools = check_env( source.env )
        
        if failed_tools:
            pass
        else:
            sources = source.sources
            if source in sources:
                return failed_tools
            failed_tools = check( sources )
            if failed_tools:
                break
    
    return failed_tools

def check_env( env ):
    """
    check if all of the tools in an envoronment exist
    """
    usedtools = set( env['TOOLS'] )
    failed_tools = set()
    for toolname in usedtools:
        tool = _Tool( toolname )
        if not tool.exists( env ):
            failed_tools.add( toolname )
        
    return failed_tools


