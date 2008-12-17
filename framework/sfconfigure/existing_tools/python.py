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
"""
Doc of RSF tool
"""

#import rsfconfig
from os.path import join
import distutils.sysconfig

def generate( env ):    
    """
    Doc of generate
    """
    vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS', 'CCSHARED', 'LDSHARED', 'SO')
    for i in range(len(vars)):
        if vars[i] is None:
            vars[i] = ""
            
    (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext) = vars
    env['LDMODULE'] = cc
    env['CC'] =cc
    env['CXX'] =cxx
#    print ldshared
    ldshared  = ldshared.split()
     
#    print ldshared
#    env['SHLINK']=ldshared[0],
#    env['SHLINKFLAGS']=ldshared[1:],
    env['LDMODULEPREFIX']=""
    env['LDMODULESUFFIX']=so_ext

    env.MergeFlags( basecflags + " " + opt )
#    import pdb; pdb.set_trace()
    env.Append(  CPPPATH=[distutils.sysconfig.get_python_inc()],
                 LDMODULEFLAGS=ldshared[1:] 
            )


    
def exists( env ):
    """
    Doc of exist
    """
    return True
    

