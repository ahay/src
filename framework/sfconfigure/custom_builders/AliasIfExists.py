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

import sfconfigure

def AliasIfExists( env, alias_name, source=None, **kw ):
    
    env.Replace( **kw )
    failed = sfconfigure.check_requirements.check( source )
    
    failed.update( sfconfigure.check_requirements.check_env(env) )
    
    failed = list(failed)

    if not failed:
        return env.Alias( alias_name, source )
    else:
        if env.get('WARNING',True):
            str_sources = []
            for s in source:
                try:
                    name = s.name
                except:
                    name = str(s)
                str_sources.append( name )
            print "WARNING: not adding target(s) %r to alias(es) %s because tool(s) %r do(es) not exist" %(str_sources,alias_name,failed)
        return None
        
#        env.Alias( what,       target )
#        env.Alias( "intstall", target )
