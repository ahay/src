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

from os.path import join
"""
 \page sconsbuilders 
 \arg \subpage sf_install_main
     install a slab main program
     
 \arg \subpage sf_install_libs
     install a slab c,c++ or python library
     
 \arg \subpage sf_install_inc
     install slab c/c++ include files or directory 
 
 \page sf_install_main sf_install_main
     Nothing in here yet

 \page sf_install_libs sf_install_libs
     Nothing in here yet

 \page sf_install_inc sf_install_inc
     Nothing in here yet
     

"""
from os.path import join as _join
    
    
def InstallLibrary( nenv, source, **kw):
    env = nenv.Clone()
    env.Replace( **kw )
    
    prefix = env['lib_prefix']
    
    inst = env.Install( target=prefix, source=source )
    
    env.AliasIfExists( ["build","buildlib"], source )
    env.AliasIfExists( ["lib","install"], inst )
    
    return inst 

def InstallExecutable( nenv, source, **kw):
    
    env = nenv.Clone()
    env.Replace( **kw )
    
    prefix = env['bin_prefix']
    
    inst = env.Install( target=prefix, source=source )
    
#    env.Depends(source,"install")
    env.AliasIfExists( "buildbin", source )
    env.AliasIfExists( ["bin","main"], inst )
    
    return inst 

def InstallInclude( nenv, source, **kw):
    
    env = nenv.Clone()
    env.Replace( **kw )
    
    prefix = env['include_prefix']
    
    inst = env.Install( target=prefix, source=source )
    
#    env.Depends(source,"install")
    env.AliasIfExists( "buildlib", source )
    env.AliasIfExists( ["lib","install"], inst )
    
    return inst 

def InstallPythonExecutable( nenv, source, **kw):
    
    env = nenv.Clone()
    env.Replace( **kw )
    
    insts =[]
    for src in nenv.Flatten( source ):
        
        Src = env.File(src)
        target = env.subst( "${bin_prefix}/${PROGPREFIX}${SOURCE.filebase}" , source=Src )
        
        inst = env.InstallAs( target=target, source=Src )
        insts.append(inst)
    
#    env.Depends(source,"install")
    env.AliasIfExists( "buildbin", source )
    env.AliasIfExists( ["bin","main"], insts )
    
    return insts


def InstallPythonModule( env, source, **kw):
#    qq = env.subst(source)
#    
#    if qq == 'rsf.py':
#        import pdb;pdb.set_trace()
    nenv = env.Clone()
    nenv.Replace( **kw )
    
    python_prefix = nenv['python_prefix']
    
    source_c = nenv.Pycompile( source )
    inst = nenv.Install( target=python_prefix, source=source_c )
    

    nenv.AliasIfExists( ["buildlib","build","buildpython"], source )
    nenv.AliasIfExists( ["lib","python","install"], inst )
    
    return inst 

def InstallPythonPackage( nenv, source, **kw):

    env = nenv.Clone()
    env.Replace( **kw )
    
    if isinstance(source, list):
        result = []
        for item in source:
            result.append( InstallPythonPackage(env, item, **kw) )
        return result
    
    python_prefix = env['python_prefix']
    
    package_name = env.get('package_name', None )
    
    if package_name is None:
        package_name = env.Dir(source).path 
        
#    source_pack = source_root.replace(package_root, "" ).strip("/")
    
    target = _join( python_prefix, package_name )
    
    py_sources = env.Glob( _join(env.Dir(source).path,"*.py") )
    env.Clean( py_sources, env.Glob( _join(env.Dir(source).path,"*.pyc") ) )
    env.Clean( target, env.Glob( _join( python_prefix, package_name,  "*.pyc" ) ) )
           
    inst = env.Install( target=target, source=py_sources )
    
    env.AliasIfExists( ["build","buildpython"], source )
    env.AliasIfExists( ["python","install"], inst )
    
    return inst 

def InstallTool( env, source, **kw ):
    
    nenv = env.Clone()
    nenv.Replace( **kw )
    
    tool_prefix = nenv['tool_prefix']
    
    inst = nenv.Install( target=tool_prefix, source=source )
    
    nenv.AliasIfExists( ["buildtool","build"], source )
    nenv.AliasIfExists( ["tool","install"], inst )
    
    return inst 
