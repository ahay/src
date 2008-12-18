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

from SCons.Script import Flatten
from os.path import join as _join


#import sconfig
ToolTemplate = '''
"""
This Tool was automatically generated, do not edit!
"""

my_appends = %(appends)r
my_prepends = %(prepends)r  
my_sets =    %(sets)r

def generate( env ):
    env.Append( **my_appends )
    env.Prepend( **my_prepends )
    env.Replace( **my_sets )

def exists( env ):
    return %(exists)r

'''

class ToolCreator(object):
    
    def __init__(self,name,dest='.'):
        self.name = name
        self.dest = dest
        self.appends = {}
        self.prepends = {}
        self.sets = {}
        self.exists = 0
    
    def Append(self,**kw):
        for key,value in kw.iteritems():
            fvalue = Flatten( value)
            lst = self.appends.setdefault(key,[])
            lst.extend(fvalue)
#        self.appends.update(kw)

    def Prepend(self,**kw):
        for key,value in kw.iteritems():
            fvalue = Flatten( value)
            lst = self.prepends.setdefault(key,[])
            lst.extend(fvalue)

    def Replace(self,**kw):
        self.sets.update(kw)
        
    def Exists(self,val=True):
        self.exists = val
        
    def UpdateExists(self, val ):
        if not self.exists:
            return
        if not val:
            self.exists = False
        return
        
    def CreateTool( self, env ):
        from sfconfigure import DoConfig
        tool_p = _join( self.dest,"%(name)s.py" %self.__dict__)
        
        env.Clean(".", tool_p)
        env.Clean(".", env.Glob("*.pyc") )

        Va = env.Value(self.appends)
        Vp = env.Value(self.prepends)
        Vs = env.Value(self.sets)
        Ve = env.Value(self.exists)

#        return env.Command( tool_p, [Va,Vp,Vs,Ve] , self.CreateToolAction)
        if DoConfig(env):
            toolfile = open( tool_p, 'w')
            toolfile.write( ToolTemplate %self.__dict__ )
            
    def CreateToolAction(self,target,source,env):
        
        toolfile = open( str(target[0]), 'w' )
        toolfile.write( ToolTemplate %self.__dict__ )
    
