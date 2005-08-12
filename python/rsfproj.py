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

import os, stat, sys, types 
import re, string, urllib, ftplib

try:
    import filecmp
except:
    import cmp

import rsfdoc
import rsfprog
import rsfconf

import SCons

# The following adds all SCons SConscript API to the globals of this module.
if SCons.__version__ == '0.96.90':
    from SCons.Script import *
else:
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

##############################################################################
# BEGIN CONFIGURATION VARIABLES
##############################################################################

#prefix for rsf commands
sfprefix = 'sf'
#prefix for vpl commands
plprefix = 'vp'
#suffix for rsf files
sfsuffix = '.rsf'
# suffix for vplot files
vpsuffix = '.vpl'

# path for binary files
datapath = os.environ.get('DATAPATH')
if not datapath:
    try:
        file = open('.datapath','r')
    except:
        try:
            file = open(os.path.join(os.environ.get('HOME'),'.datapath'),'r')
        except:
            file = None
    if file:
        for line in file.readlines():
            check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],line)
            if check:
                datapath = check.group(1)
        file.close()
    if not datapath:
        datapath = './' # the ultimate fallback
dataserver = os.environ.get('RSF_DATASERVER','ftp://egl.beg.utexas.edu/')
tmpdatapath = os.environ.get('TMPDATAPATH',datapath)

# directory tree for executable files
top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')
libdir = os.path.join(top,'lib')
incdir = os.path.join(top,'include')
figdir = os.environ.get('RSFFIGS',os.path.join(top,'figs'))

resdir = None

def set_dir(ref='.',dir='Fig'):
     global resdir
     resdir = os.path.join(ref,dir)

set_dir()

# temporary (I hope)
sep = os.path.join(os.environ.get('SEP',''),'bin/')

#############################################################################
# CUSTOM BUILDERS
#############################################################################

def test(target=None,source=None,env=None):
    src = str(source[0])
    locked = re.sub('.*\/([^\/]+)\/([^\/]+)\/([^\/]+)\/Fig\/',
                    figdir+'/\\1/\\2/\\3/',os.path.abspath(src))
    print "Comparing %s and %s" % (locked,src)
    if os.path.isfile(locked):
        try:
            if not filecmp.cmp(locked,src,shallow=0):
                return 1
        except:
            if not cmp.cmp(locked,src):
                return 1
    else:
        print 'No locked file "%s" ' % locked
    return 0

def retrieve(target=None,source=None,env=None):
    "Fetch data from the web"
    global dataserver
    top = env.get('top')
    folder = top + os.sep +env['dir']
    private = env.get('private')
    if private:
        login = private['login']
        password = private['password']
        server = private['server']
        try:
            session = ftplib.FTP(server,login,password)
            session.cwd(folder)
        except:
            print 'Could not establish connection with "%s/%s" ' % (server,
                                                                    folder)
            return 3
        for file in map(str,target):
            remote = os.path.basename(file)
            try:
                 download = open(file,'wb')
                 session.retrbinary('RETR '+remote,
                                    lambda x: download.write(x))
                 download.close()
            except:
                 print 'Could not download file "%s" ' % file
                 return 1
            if not os.stat(file)[6]:
                print 'Could not download file "%s" ' % file
                os.unlink(file)
                return 4
        session.quit()
    else:
        for file in map(str,target):
            remote = os.path.basename(file)
            rdir =  string.join([dataserver,folder,remote],'/')
            try:
                urllib.urlretrieve(rdir,file)

                if not os.stat(file)[6]:
                    print 'Could not download file "%s" ' % file
                    os.unlink(file)
                    return 2
            except:
                print 'Could not download "%s" from "%s" ' % (file,rdir)
                return 5
    return 0

View = Builder(action = sep + "xtpen $SOURCES",src_suffix=vpsuffix)
Retrieve = Builder(action = Action(retrieve,varlist=['dir','private','top']))
Test = Builder(action=Action(test))

#############################################################################
# PLOTTING COMMANDS
#############################################################################

combine ={
    'SideBySideAniso': lambda n:
    sep + "vppen yscale=%d vpstyle=n gridnum=%d,1 $SOURCES" % (n,n),
    'OverUnderAniso': lambda n:
    sep + "vppen xscale=%d vpstyle=n gridnum=1,%d $SOURCES" % (n,n),
    'SideBySideIso': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=%d,1 $SOURCES" % n,
    'OverUnderIso': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=1,%d $SOURCES" % n,
    'TwoRows': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=%d,2 $SOURCES" % (n/2),
    'Overlay': lambda n:
    sep + "vppen erase=o vpstyle=n $SOURCES",
    'Movie': lambda n:
    sep + "vppen vpstyle=n $SOURCES"
    }

#############################################################################
class Project(Environment):
    def __init__(self,**kw):
        apply(Environment.__init__,(self,),kw)
        self.EnsureSConsVersion(0,96)
        opts = Options(os.path.join(libdir,'rsfconfig.py'))
        rsfconf.options(opts)
        opts.Add('TIMER','Whether to time execution')
        opts.Update(self)
        cwd = os.getcwd()
        dir = os.path.basename(cwd)
        if datapath[:2] == './':
            self.path = datapath
        else:
            self.path = datapath + dir + os.sep
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self.SConsignFile(self.path+'.sconsign')
        self.resdir = resdir
        self.figdir = re.sub('.*\/((?:[^\/]+)\/(?:[^\/]+)\/(?:[^\/]+))$',figdir+'/\\1',cwd)
        self.progsuffix = self['PROGSUFFIX']
        self.Append(ENV={'DATAPATH':self.path,
                         'TMPDATAPATH': tmpdatapath,
                         'PYTHONPATH': os.environ.get('PYTHONPATH',libdir), 
                         'XAUTHORITY':
                         os.path.join(os.environ.get('HOME'),'.Xauthority'),
                         'DISPLAY': os.environ.get('DISPLAY'),
                         'VPLOTFONTDIR': os.environ.get('VPLOTFONTDIR'),
                         'HOME': os.environ.get('HOME'),
                         'RSFROOT':top},
                    BUILDERS={'View':View,
                              'Retrieve':Retrieve,
                              'Test':Test},
                    LIBPATH=[libdir],
                    CPPPATH=[incdir],
                    LIBS=['rsf','m'],
                    PROGSUFFIX='.exe')
        if sys.platform[:6] == 'cygwin':
            self['ENV']['PATH'] = self['ENV']['PATH'] + ':/usr/X11R6/bin'
            self['ENV']['SYSTEMROOT'] = os.environ.get('SYSTEMROOT')
        self['PROGPREFIX']=''
        self.view = []
        self.lock = []
        self.test = []
        self.coms = []
        sys.path.append('../../../packages')
    def Exe(self,source,**kw):
        target = source.replace('.c','.x')
        return apply(self.Program,(target,source),kw)
    def Flow(self,target,source,flow,stdout=1,stdin=1,rsf=1,
             suffix=sfsuffix,prefix=sfprefix,src_suffix=sfsuffix):
        if not flow:
            return None        
        sources = []
        if source:
            if type(source) is types.ListType:
                files = source
            else:
                files = string.split(source)
            for file in files:
                if ('.' not in file):
                    file = file + src_suffix
                sources.append(file)
        else:
            stdin=0
        lines = string.split(flow,'&&')
        steps = []
        for line in lines:
            substeps = []
            sublines = string.split(line,'|')
            for subline in sublines:           
                pars = string.split(subline)
                # command is assumed to be always first in line
                command = pars.pop(0)
                # check if this command is in our list
                if rsf:
                    rsfprog = prefix + command            
                    if rsfdoc.progs.has_key(rsfprog):
                        command = os.path.join(bindir,rsfprog+self.progsuffix) 
                        sources.append(command)
                        if rsfprog not in self.coms:
                            self.coms.append(rsfprog)
                elif re.match(r'[^/]+\.exe$',command): # local program
                    command = os.path.join('.',command)                    
                #<- assemble the command line
                pars.insert(0,command)
                substeps.append(string.join(pars,' '))
            #<-
            steps.append(string.join(substeps," | "))
        #<- assemble the pipeline
        command = string.join(steps," && ")
        if stdout==1:
            command = command + " > $TARGET"
        elif stdout==0:
            command = command + " >/dev/null"
        if stdin:
            command = "< $SOURCE " + command
        timer = self.get('TIMER')
        if timer and timer[0] != 'n' and timer[0] != '0':
            command = WhereIs('time') + ' ' + command
        targets = []
        if type(target) is types.ListType:
            files = target
        else:
            files = string.split(target)
        for file in files:
            if not re.search(suffix + '$',file):
                file = file + suffix
            targets.append(file)
        if suffix == sfsuffix:            
            datafiles = [] 
            for target in targets:
                if os.sep not in target:
                    datafile = self.path + target + '@'
                    datafiles.append(datafile)
            targets = targets + datafiles
        return self.Command(targets,sources,command)
    def Plot (self,target,source,flow=None,suffix=vpsuffix,vppen=None,
              view=None,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        if combine.has_key(flow):
            if not type(source) is types.ListType:
                source = string.split(source)
            flow = apply(combine[flow],[len(source)])
            if vppen:
                flow = flow + ' ' + vppen
            kw.update({'src_suffix':vpsuffix,'stdin':0})
        if view:
            flow = flow + ' | %s pixmaps=y' % (sep+'xtpen')
            kw.update({'stdout':-1})
        kw.update({'suffix':suffix})        
        return apply(self.Flow,(target,source,flow),kw)
    def Result(self,target,source,flow=None,suffix=vpsuffix,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        target2 = os.path.join(self.resdir,target)
        plot = apply(self.Plot,(target2,source,flow),kw)
        self.Default (plot)
        self.view.append(self.View(target + '.view',plot))
        lock = self.InstallAs(os.path.join(self.figdir,target+suffix),
                              target2+suffix)
        
        self.lock.append(lock)
        self.Alias(target + '.lock',lock)
        test = self.Test('.test_'+target,target2+suffix)
        self.test.append(test)
        self.Alias(target + '.test',test)
        return plot
    def End(self):
        if self.view: # if any results
            self.Alias('view',self.view)
            self.Alias('lock',self.lock)
            self.Alias('test',self.test)
        self.Command('.sf_uses',None,'echo %s' % string.join(self.coms,' '))
    def Fetch(self,file,dir,private=None,top='data'):
        return self.Retrieve(file,None,dir=dir,private=private,top=top)

# Default project
project = Project()
def Flow(target,source,flow,**kw):
    return apply(project.Flow,(target,source,flow),kw)
def Plot (target,source,flow=None,**kw):
    return apply(project.Plot,(target,source,flow),kw)
def Result(target,source,flow=None,**kw):
    return apply(project.Result,(target,source,flow),kw)
def Fetch(file,dir,private=0):
    return project.Fetch(file,dir,private)
def Exe(source,**kw):
    return apply(project.Exe,[source],kw)
def End(**kw):
    return apply(project.End,[],kw)
def Program(*arg,**kw):
    return apply(project.Program,arg,kw)
def Get(name):
    return project['ENV'].get(name)
def Program90(prog):
    sources = [prog]
    rsfconf.depends90(project,sources,prog)
    return project.Program(prog,
                           map(lambda x: x + '.f90',sources),
                           F90PATH=[incdir],
                           LIBS=['rsff90','rsf','m'],
                           LINK=project.get('F90'),
                           LINKFLAGS=project.get('F90FLAGS'))

if __name__ == "__main__":
     import pydoc
     pydoc.help(Project)
     
#   $Id$
