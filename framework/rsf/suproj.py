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

from __future__ import print_function, division, absolute_import
import rsf.proj, os, re, string, sys, array, types

susuffix = '.su'
pssuffix = '.eps'

topdir = os.environ.get('CWPROOT','')
bindir = os.path.join(topdir,'bin')

try:
    suprogs = os.listdir(bindir)
except:
    print()
    sys.stderr.write("It looks like SU is not installed.\n")
    sys.stderr.write("exception from suprogs=listdir(bindir) (bindir=%s)\n"%
                     bindir)
    sys.stderr.write("bindir is topdir(=%s) with /bin appended.\n"%topdir)
    sys.stderr.write("bindir is os.environ.get('CWPROOT','').\n")
    sys.stderr.write("Maybe CWPROOT not set in .bash_profile.\n\n")
    sys.exit(1)

suplots = ['plot']
for prog in suprogs:
    if prog[0] == 'x' and 'ps'+prog[1:] in suprogs:
        suplots.append(prog[1:])
re_plots = re.compile(r'\b(su|s|)(?:x|ps)?(%s)\b' % '|'.join(suplots))

class SUProject(rsf.proj.Project):
    def __init__(self,**kw):
        rsf.proj.Project.__init__(*(self,), **kw)
        self['ENV']['PATH'] = self['ENV']['PATH'] + ':' + bindir
        self['ENV']['CWPROOT'] = topdir
        self.plots = []
        self.views = []
        self.bindir = bindir
    def Flow(self,target,source,flow,suffix=susuffix,src_suffix=susuffix,**kw):
        kw.update({'rsfflow':0,'suffix': suffix,'src_suffix':src_suffix})
        return rsf.proj.Project.Flow(*(self,target,source,flow), **kw)
    def Plot(self,target,source,flow=None,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        if flow == 'Merge':
            if not type(source) is list:
                sources = source.split()
                source = [x+pssuffix for x in sources]
            self.Command(target+pssuffix,source,
                         '%s %s > $TARGET' % \
                         (os.path.join(bindir,'psmerge'),
                          ' '.join([' in=${SOURCES[%d]}' % x for x in range(len(source))])))
        else:
            # X output
            xflow  = re_plots.sub('\\1x\\2',flow)
            kw.update({'suffix':'.x11','stdout':-1})
            self.Flow(*(target,source,xflow), **kw)
            # Postscript output
            psflow = re_plots.sub('\\1ps\\2',flow)
            kw.update({'suffix': pssuffix,'stdout':1})
            self.Flow(*(target,source,psflow), **kw)
    def Result(self,target,source,flow=None,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        target2 = os.path.join(self.resdir,target)
        self.Plot(*(target2,source,flow), **kw)
        self.Default (target2+pssuffix)
        if flow != 'Merge':
            self.Alias(target+'.view',target2+'.x11')
            self.views.append(target)
        self.plots.append(target)

        locked = os.path.join(self.figdir,target+pssuffix)
        self.InstallAs(locked,target2+pssuffix)
        self.Alias(target + '.lock',locked)
        self.lock.append(locked)

        test = self.Test('.test_'+target,target2+pssuffix,
                         figdir=self.figs,bindir=self.bindir)
        self.test.append(test)
        self.Alias(target + '.test',test)
    def End(self):
        self.Command('.suproj',self.lock,
                     action=self.Action(self.Info,var='su,suffix'),
                     su=1,suffix='.su')
        if self.plots: # if any results
            self.Alias('lock',self.lock+['.suproj'])
            self.Alias('test',self.test)
        if self.views:
            self.Alias('view',[x+'.view' for x in self.views])
        else:
            self.Command('view',None,'echo "There is nothing to view" ')

def little_endian():
    "check for endianness"
    return ord(array.array("i",[1]).tostring()[0])

# Default project
project = SUProject()
def Flow(target,source,flow,**kw):
    return project.Flow(*(target,source,flow), **kw)
def Plot(target,source,flow=None,**kw):
    return project.Plot(*(target,source,flow), **kw)
def Result(target,source,flow=None,**kw):
    return project.Result(*(target,source,flow), **kw)
def Fetch(file,dir,private=0,server='ftp://ftp.cwp.mines.edu/pub',**kw):
    return project.Fetch(*(file,dir,private,server), **kw)
def End():
    return project.End()
