##   Copyright (C) 2004 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import rsfproj, os, re, string

susuffix = '.su'
pssuffix = '.eps'

topdir = os.environ.get('CWPROOT')
bindir = os.path.join(topdir,'bin')

suprogs = os.listdir(bindir)
suplots = []
for prog in suprogs:
    if prog[0] == 'x' and 'ps'+prog[1:] in suprogs:
        suplots.append(prog[1:])
re_plots = re.compile(r'\b(%s)\b' % string.join(suplots,'|'))

class SUProject(rsfproj.Project):
    def __init__(self,**kw):
        apply(rsfproj.Project.__init__,(self,),kw)
        self['ENV']['PATH'] = self['ENV']['PATH'] + ':' + bindir
    def Flow(self,target,source,flow,suffix=susuffix,src_suffix=susuffix,**kw):
        kw.update({'rsf':0,'suffix': suffix,'src_suffix':src_suffix})
        return apply(rsfproj.Project.Flow,(self,target,source,flow),kw)
    def Plot(self,target,source,flow=None,**kw):
        if not flow: # two arguments
            flow = source
            source = target

        xflow  = re_plots.sub('x\\1',flow)
        kw.update({'suffix':'.view','stdout':-1})
        apply(self.Flow,(target,source,xflow),kw)
        
        psflow = re_plots.sub('ps\\1',flow)
        kw.update({'suffix': pssuffix,'stdout':1})
        apply(self.Flow,(target,source,psflow),kw)

# Default project
project = SUProject()
def Flow(target,source,flow,**kw):
    return apply(project.Flow,(target,source,flow),kw)
def Plot(target,source,flow=None,**kw):
    return apply(project.Plot,(target,source,flow),kw)
def End(**kw):
    return apply(project.End,[],kw)
