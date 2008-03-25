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

import os, stat, sys, types, copy
import re, string, urllib, ftplib

# database for .sconsign
try:
    import dbhash
    db='dbhash'
except:
    try:
        import gdbm
        db='gdbm'
    except:
        db=None
        
try:
    import filecmp
except:
    import cmp

import rsfdoc
import rsfprog
import rsfconf

import SCons

# The following adds all SCons SConscript API to the globals of this module.
version = map(int,string.split(SCons.__version__,'.'))
if version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
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
        pathfile = open('.datapath','r')
    except:
        try:
            pathfile = open(os.path.join(os.environ.get('HOME'),'.datapath'),
                            'r')
        except:
            pathfile = None
    if pathfile:
        for line in pathfile.readlines():
            check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],line)
            if check:
                datapath = check.group(1)
        pathfile.close()
    if not datapath:
        datapath = './' # the ultimate fallback
dataserver = os.environ.get('RSF_DATASERVER','ftp://egl.beg.utexas.edu/')
tmpdatapath = os.environ.get('TMPDATAPATH',datapath)
# split into directory
pathdir = os.path.dirname(datapath)
pathbase = os.path.basename(datapath)

# directory tree for executable files
top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')
libdir = os.path.join(top,'lib')
incdir = os.path.join(top,'include')
figdir = os.environ.get('RSFFIGS',os.path.join(top,'figs'))

libs = os.environ.get('LIBS',"")

resdir = None

def set_dir(dir='Fig'):
     global resdir
     resdir = dir

set_dir()

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

def symlink(target, source, env):
    "Create symbolic link"
    src = str(source[0])
    tar = str(target[0])
    os.symlink(src,tar)
    return 0

def retrieve(target=None,source=None,env=None):
    "Fetch data from the web"
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
        server = env.get('server')
        for file in map(str,target):
            remote = os.path.basename(file)
            rdir =  string.join([server,folder,remote],'/')
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

vppen = os.path.join(bindir,'vppen')
xtpen = os.path.join(bindir,'xtpen')
pspen = os.path.join(bindir,'pspen')

View  = Builder(action = xtpen + " $SOURCES",src_suffix=vpsuffix)

printer = os.environ.get('PSPRINTER',os.environ.get('PRINTER','postscript'))

Print = Builder(action = pspen + " printer=%s $SOURCES" % printer,src_suffix=vpsuffix)
Retrieve = Builder(action = Action(retrieve,
                                   varlist=['dir','private','top','server']))
Test = Builder(action=Action(test))

#############################################################################
# PLOTTING COMMANDS
#############################################################################

combine ={
    'SideBySideAniso': lambda n:
        vppen + " yscale=%d vpstyle=n gridnum=%d,1 $SOURCES" % (n,n),
    'OverUnderAniso': lambda n:
        vppen + " xscale=%d vpstyle=n gridnum=1,%d $SOURCES" % (n,n),
    'SideBySideIso': lambda n:
        vppen + " size=r vpstyle=n gridnum=%d,1 $SOURCES" % n,
    'OverUnderIso': lambda n:
        vppen + " size=r vpstyle=n gridnum=1,%d $SOURCES" % n,
    'TwoRows': lambda n:
        vppen + " size=r vpstyle=n gridnum=%d,2 $SOURCES" % (n/2),
    'TwoColumns': lambda n:
        vppen + " size=r vpstyle=n gridnum=2,%d $SOURCES" % (n/2),
    'Overlay': lambda n:
        vppen + " erase=o vpstyle=n $SOURCES",
    'Movie': lambda n:
        vppen + " vpstyle=n $SOURCES"
    }

# Environmental variables to pass to SCons
keepenv = ('DISPLAY','VPLOTFONTDIR','HOME','LD_LIBRARY_PATH','RSFMEMSIZE')

#############################################################################
class Project(Environment):
    def __init__(self,**kw):
        apply(Environment.__init__,(self,),kw)
        self.EnsureSConsVersion(0,96)
        opts = Options(os.path.join(libdir,'rsfconfig.py'))
        rsfconf.options(opts)
        opts.Add('TIMER','Whether to time execution')
        opts.Add('CHECKPAR','Whether to check parameters')
        opts.Add('ENVIRON','Additional environment settings')
        opts.Add('CLUSTER','Nodes available on a cluster')
        opts.Update(self)
        cwd = os.getcwd()
        self.cwd = cwd
        self.path = pathdir
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        if datapath[:2] != './':
            # create a hierarcical structure
            (book,chap,proj) = (os.path.basename(
                os.path.dirname(os.path.dirname(cwd))),
                                os.path.basename(os.path.dirname(cwd)),
                                os.path.basename(cwd))
            for level in (book,chap,proj):
                if level:
                    self.path = os.path.join(self.path,level)
                    if not os.path.exists(self.path):
                        os.mkdir(self.path)
        self.path = os.path.join(self.path,pathbase)
        if db=='gdbm':
            self.SConsignFile(self.path+'.sconsign.'+db,gdbm)
        elif db=='dbhash':
            self.SConsignFile(self.path+'.sconsign.'+db,dbhash)
        else:
            self.SConsignFile(self.path+'.sconsign')
        self.resdir = resdir
        self.figdir = re.sub('.*\/((?:[^\/]+)\/(?:[^\/]+)\/(?:[^\/]+))$',
                             figdir+'/\\1',cwd)
        self.progsuffix = self['PROGSUFFIX']
        for env in keepenv:
            getenv = os.environ.get(env)
            if getenv:
                self.Append(ENV={env:getenv})
        if sys.platform[:6] == 'cygwin':
            exe = ''
        else:
            exe = '.exe'
        self.Append(ENV={'DATAPATH':self.path,
                         'TMPDATAPATH': tmpdatapath,
                         'PYTHONPATH': os.environ.get('PYTHONPATH',libdir), 
                         'XAUTHORITY':
                         os.environ.get('XAUTHORITY',
                                        os.path.join(os.environ.get('HOME'),
                                                     '.Xauthority')),
                         'RSFROOT':top},
                    BUILDERS={'View':View,
                              'Print':Print,
                              'Retrieve':Retrieve,
                              'Test':Test},
                    LIBPATH=[libdir],
                    CPPPATH=[incdir],
                    LIBS=[libs],
                    PROGSUFFIX=exe)
        self.Prepend(LIBS=['rsf'])
        if sys.platform[:6] == 'cygwin':
            self['ENV']['PATH'] = self['ENV']['PATH'] + ':/usr/X11R6/bin'
            self['ENV']['SYSTEMROOT'] = os.environ.get('SYSTEMROOT')
        self['PROGPREFIX']=''
        self.view = []
        self.prnt = []
        self.lock = []
        self.test = []
        self.coms = []
        self.data = []
        sys.path.append('../../../packages')

        timer = self.get('TIMER')
        if timer and timer[0] != 'n' and timer[0] != '0':
            self.timer = WhereIs('time') + ' '
        else:
            self.timer = ''

        checkpar = self.get('CHECKPAR')
        self.checkpar = checkpar and checkpar[0] != 'n' and checkpar[0] != '0'

        self.environ = self.get('ENVIRON','')

        
        self.jobs = GetOption('num_jobs')
        cluster = self.get('CLUSTER','localhost 1')
        hosts = string.split(cluster)
        self.nodes = []
        for i in range(1,len(hosts),2):
            nh = int(hosts[i])
            self.nodes.extend([hosts[i-1]]*nh)
        self.ip = 0

        # self.nodes is a list of CPUs
        # self.jobs is the number of jobs
        # self.ip is the current CPU

        for key in self['ENV'].keys():
            self.environ = self.environ + ' %s=%s' % (key,self['ENV'][key]) 

    def Flow(self,target,source,flow,stdout=1,stdin=1,rsf=1,
             suffix=sfsuffix,prefix=sfprefix,src_suffix=sfsuffix,
             split=[],axis=[3,1],reduce='cat',local=0):

        if not flow:
            return None     

        if type(target) is types.ListType:
            tfiles = target
        else:
            tfiles = string.split(target)

        if source:
            if type(source) is types.ListType:
                sfiles = source
            else:
                sfiles = string.split(source)
        else:
            sfiles = []

        if split and self.jobs > 1 and rsf and sfiles:
            # Split the flow into parallel flows

            if self.jobs < axis[1]:
                jobs = self.jobs            
                w = axis[1]/self.jobs # length of one chunk
            else:
                jobs = axis[1]
                w = 1
                
            par_sfiles = copy.copy(sfiles)
            par_targets = {}
            for tfile in tfiles:
                par_targets[tfile] = []
                
            for i in range(jobs):
                for j in split:
                    source = sfiles[j] + '__' + str(i)
                    par_sfiles[j] = source

                    if i==jobs-1: # last chunk
                        self.Flow(source,sfiles[j],
                                  'window f%d=%d squeeze=n' % (axis[0],i*w))
                    else:
                        self.Flow(source,sfiles[j],
                                  'window n%d=%d f%d=%d squeeze=n' % 
                                  (axis[0],w,axis[0],i*w))
                

                par_tfiles = []
                for j in range(len(tfiles)):
                    tfile = tfiles[j]
                    par_tfile = tfile + '__' + str(i)
                    
                    par_tfiles.append(par_tfile)
                    par_targets[tfile].append(par_tfile)
 
                # operation on one chunk    
                self.Flow(par_tfiles,par_sfiles,flow,
                          stdout,stdin,1,
                          suffix,prefix,src_suffix)

            # Reduce parallel TARGETS down to original TARGETS:
            for tfile in tfiles:
                self.Flow(tfile,par_targets[tfile],
                          '%s axis=%d ${SOURCES[1:%d]}' % \
                          (reduce,axis[0],jobs),local=1)
            return

        sources = []
        if sfiles:
            for file in sfiles:
                if ('.' not in file):
                    file = file + src_suffix
                sources.append(file)
        else:
            stdin=0

        # May need to do it remotely
        if local:
            node = 'localhost'
        else: # get it from the rotating list
            node = self.nodes[self.ip]
            self.ip = self.ip + 1
            if self.ip == len(self.nodes):
                self.ip = 0
                
        if node != 'localhost':
            remote = '%s %s ' % (WhereIs('env'),self.environ)
        else:
            remote = ''
            
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
                    if command[:2]==prefix:
                        # assuming prefix is two chars (sf)
                        rsfprog = command
                    else:
                        rsfprog = prefix + command            
                    if rsfdoc.progs.has_key(rsfprog):
                        if self.checkpar:
                            for par in pars:
                                if rsfdoc.progs[rsfprog].check(par):
                                    sys.stderr.write('Failed on "%s"\n' % 
                                                     subline)
                                    sys.exit(1)
                        command = os.path.join(bindir,rsfprog+self.progsuffix) 
                        sources.append(command)
                        if rsfprog not in self.coms:
                            self.coms.append(rsfprog)
                else:
                    rsfprog = None
                    if re.match(r'[^/]+\.exe$',command): # local program
                        command = os.path.join('.',command)
                pars.insert(0,command)
                # special rule for solvers
                if rsfprog == prefix+'conjgrad' or \
                       rsfprog == prefix+'cconjgrad':
                    command = pars.pop(1)
                    # check if this command is in our list
                    if rsf:
                        if command[:2]==prefix:
                            # assuming prefix is two chars
                            rsfprog = command
                        else:
                            rsfprog = prefix + command            
                        if rsfdoc.progs.has_key(rsfprog):
                            command = os.path.join(bindir,
                                                   rsfprog+self.progsuffix) 
                            sources.append(command)
                            if rsfprog not in self.coms:
                                self.coms.append(rsfprog)
                    else:
                        rsfprog = None
                        if re.match(r'[^/]+\.exe$',command): # local program
                            command = os.path.join('.',command)                 
                    pars.insert(1,command)
                #<- assemble the command line
                substep = remote + string.join(pars,' ')
                substeps.append(substep)
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
        command = self.timer + command

        # May need to do it remotely
        if remote:
            command = re.sub('"','\\"',command)
            command = string.join([WhereIs('ssh'),node,'\"cd ',
                                   self.cwd,';',command,'\"'])
                        
        targets = []
        for file in tfiles:
            if (not re.search(suffix + '$',file)) and ('.' not in file):
                file = file + suffix
            targets.append(file)
        if suffix == sfsuffix and \
               sys.platform[:6] != 'cygwin' and \
               sys.platform[:7] != 'interix':            
            targets = targets + \
                      map(lambda x, self=self: self.path + x + '@',
                          filter(lambda x, suffix=suffix:
                                 x[-4:] == suffix,targets))
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
            flow = flow + ' | %s pixmaps=y' % xtpen
            kw.update({'stdout':-1})
        kw.update({'suffix':suffix})        
        return apply(self.Flow,(target,source,flow),kw)
    def Result(self,target,source,flow=None,suffix=vpsuffix,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        target2 = os.path.join(self.resdir,target)
        plot = apply(self.Plot,(target2,source,flow),kw)
        target2 = target2 + suffix
        self.Default (plot),
        self.view.append(self.View(target + '.view',plot))
        self.prnt.append(self.Print(target + '.print',plot))
        locked = os.path.join(self.figdir,target+suffix)
        self.InstallAs(locked,target2)
        lock2 = self.Command(target2+'@',locked,symlink)
        self.Alias(target + '.lock',lock2)
        self.lock.append(lock2)
        self.Command(target + '.flip',target2,
                     '%s $SOURCE %s' % (xtpen,locked))
        test = self.Test('.test_'+target,target2)
        self.test.append(test)
        self.Alias(target + '.test',test)
        return plot
    def End(self):
        if self.view: # if any results
            self.Alias('view',self.view)
            self.Alias('print',self.prnt)
            self.Alias('lock',self.lock)
            self.Alias('test',self.test)
        else:
            self.Command('test',None,'echo "Nothing to test"')
        self.Command('.sf_uses',None,'echo %s' % string.join(self.coms,' '))
        self.Command('.sf_data',None,'echo %s' % string.join(self.data,' '))
    def Fetch(self,files,dir,private=None,server=dataserver,top='data'):
        if private:
            self.data.append('PRIVATE')
        else:
            if not type(files) is types.ListType:
                files = string.split(files)
            for fil in files:
                self.data.append(os.path.join(top,dir,fil))
        return self.Retrieve(files,None,
                             dir=dir,private=private,
                             top=top,server=server)

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
