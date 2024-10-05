# Copyright (C) 2004 University of Texas at Austin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from __future__ import division, absolute_import, print_function
import os, stat, sys, types, copy, re, string, ftplib, socket, json
import rsf.conf, rsf.path, rsf.flow, rsf.prog, rsf.node
import SCons

from urllib.request import urlopen
from shutil import copyfileobj
import ssl

context = ssl.create_default_context()
context.check_hostname = False
context.verify_mode = ssl.CERT_NONE

# The following adds all SCons SConscript API to the globals of this module.
version = list(map(int,SCons.__version__.split('.')[:3]))
from SCons.Script import *
"""
if version[0] >= 1 or version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
else:
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())
"""

SCons.Defaults.DefaultEnvironment(tools = [])

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

# Convert Chinese to two-byte-Unicode string
def ch2uni(str1):
    str2=""
    for i in str1:
        if not i == "\n":
            i = i.encode("unicode_escape")
            i = str(i)
            i = i[5:-1]
            i1 = i[:2]
            i2 = i[2:]
            str2 += "\\v%d \\v%d "%(int(i1,16),int(i2,16))
        else: str2+="\\n "
    return(str2)

def get_geolocation(address=""):
    """Get geolocation from http://ip-api.com."""
    if address == "":
        url = "http://ip-api.com/json"
    else:
        url = "http://ip-api.com/json/" + address

    try:
        response = urllib_request.urlopen(url,timeout=5)
        data = json.load(response)
        return data['countryCode']
    except:
#        print("Fail to get country code")
        return None

def get_dataserver():
    'Set the default data server'
    country = get_geolocation()
    if country == "CN":
        #    dataserver = os.environ.get('RSF_DATASERVER','http://49.235.136.252')
        # return os.environ.get('RSF_DATASERVER','https://reproducibility.org')
        return os.environ.get('RSF_DATASERVER','https://ahay.org')
    else:
        # return os.environ.get('RSF_DATASERVER','https://reproducibility.org')
        return os.environ.get('RSF_DATASERVER','https://ahay.org')

dataserver = None
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
    figdir = env.get('figdir')
    bindir = env.get('bindir')

    locked = re.sub('.*\/([^\/]+)\/([^\/]+)\/([^\/]+)\/Fig\/',
                    figdir+'/\\1/\\2/\\3/',os.path.abspath(src))
    print("Comparing %s and %s" % (locked,src))
    if os.path.isfile(locked):
        if locked.endswith(vpsuffix):
            diff = os.system(' '.join([os.path.join(bindir,sfprefix+'vplotdiff'),
                                   locked,src]))
        else:
            diff = 0
        return diff
    else:
        print('No locked file "%s" ' % locked)
        return 0

def echo(target,source,env):
    obj = env.get('out','')
    if obj:
        trg = open(str(target[0]),'w')
        if type(obj) is list:
            obj = ' '.join(obj)
        trg.write(obj+'\n')
        trg.close()
    err = env.get('err','')
    if err:
        if type(err) is list:
            err = ' '.join(err)
        sys.stderr.write(err+'\n')
    return 0

def retrieve_emit(target=None, source=None, env=None):
    if sys.platform[:6] != 'cygwin':
        usedatapath = env.get('usedatapath')
    else:
        usedatapath = False
    server = env.get('server')
    if usedatapath and server != 'local':
        for file in list(map(str,target)):
            localfile=env.path+os.path.basename(file)
            target.append(localfile)
    return target, source

def retrieve(target=None,source=None,env=None):
    "Fetch data from the web"
    top = env.get('top')
    if top:
        folder = top + os.sep +env['dir']
    else:
        folder = env['dir']
    private = env.get('private')
    if sys.platform[:6] != 'cygwin':
        usedatapath = env.get('usedatapath')
    else:
        usedatapath = False
    if private:
        login = private['login']
        password = private['password']
        server = private['server']
        if not server:
            print('Cannot access proprietary data server')
            return 7
        try:
            session = ftplib.FTP(server,login,password)
            session.cwd(folder)
        except:
            print('Could not establish connection with "%s/%s" ' % (server,
                                                                    folder))
            return 3
        for file in [x for x in list(map(str,target)) if not os.path.abspath(x).startswith(env.path)]:
            remote = os.path.basename(file)
            if usedatapath:
                localfile=env.path+remote
            else:
                localfile=file
            try:
                download = open(localfile,'wb')
                session.retrbinary('RETR '+remote,
                                   lambda x: download.write(x))
                download.close()
            except:
                print('Could not download file "%s" ' % file)
                return 1
            if not os.stat(localfile)[6]: # if zero size file
                print('Could not download file "%s" ' % file)
                os.unlink(localfile)
                return 4
            if usedatapath:
                if os.path.isfile(file):
                    os.unlink(file)
                os.symlink(localfile,file)
        session.quit()
    else:
        server = env.get('server')
        if server == 'local':
            for file in list(map(str,target)):
                remote = os.path.basename(file)
                remote = os.path.join(folder,remote)
                try:
                    os.symlink(remote,file)
                except:
                    print('Could not link file "%s" ' % remote)
                    os.unlink(file)
                    return 6
        else:
            for file in [x for x in list(map(str,target)) if not os.path.abspath(x).startswith(env.path)]:
                remote = os.path.basename(file)
                rdir =  '/'.join([server,folder,remote])
                if usedatapath:
                    localfile=env.path+remote
                else:
                    localfile=file
                try:
                    with urlopen(rdir, context=context) as in_stream, open(localfile, 'wb') as out_file:
                        copyfileobj(in_stream, out_file)
                        
                    if not os.stat(localfile)[6]:
                        print('Could not download file "%s" ' % localfile)
                        os.unlink(localfile)
                        return 2
                except:
                    print('Could not download "%s" from "%s" ' % (localfile,rdir))
                    return 5
                if usedatapath:
                    if os.path.isfile(file):
                        os.unlink(file)
                    os.symlink(localfile,file)
    return 0

printer = os.environ.get('PSPRINTER',os.environ.get('PRINTER','postscript'))

Retrieve = Builder(action = Action(retrieve,
                                   varlist=['dir','private','top','server','usedatapath']),
                   emitter=retrieve_emit)
Test = Builder(action=Action(test),varlist=['figdir','bindir'])
Echo = Builder(action=Action(echo),varlist=['out','err'])

#############################################################################
# PLOTTING COMMANDS
#############################################################################

combine ={
    'SideBySideAniso': lambda v, n: v + " yscale=%d vpstyle=n gridnum=%d,1 ${SOURCES[:%d]}" % (n,n,n),
    'OverUnderAniso':  lambda v, n: v + " xscale=%d vpstyle=n gridnum=1,%d ${SOURCES[:%d]}" % (n,n,n),
    'SideBySideIso':   lambda v, n: v + " size=r vpstyle=n gridnum=%d,1 ${SOURCES[:%d]}" % (n,n),
    'OverUnderIso':    lambda v, n: v + " size=r vpstyle=n gridnum=1,%d ${SOURCES[:%d]}" % (n,n),
    'TwoRows':         lambda v, n: v + " size=r vpstyle=n gridnum=%d,2 ${SOURCES[:%d]}" % ((n+1)/2,n),
    'TwoColumns':      lambda v, n: v + " size=r vpstyle=n gridnum=2,%d ${SOURCES[:%d]}" % ((n+1)/2,n),
    'Overlay':         lambda v, n: v + " erase=o vpstyle=n ${SOURCES[:%d]}" % n,
    'Movie':           lambda v, n: v + " vpstyle=n ${SOURCES[:%d]}" % n
    }

# Environmental variables to pass to SCons
keepenv = ('DISPLAY','VPLOTFONTDIR','HOME','USER','WORK', 'SCRATCH',
           'LD_LIBRARY_PATH','DYLD_LIBRARY_PATH','MIC_LD_LIBRARY_PATH',
           'RSFMEMSIZE','PYTHONPATH','SFPENOPTS','MPICH_HOME',
           'MODULEPATH','OMPI_MCA_btl','I_MPI_ROOT')

#############################################################################
class Project(Environment):
    def __init__(self,**kw):
        Environment.__init__(*(self,), **kw)
        #        self.EnsureSConsVersion(0,96)
        opts = {
            'TIMER':'Whether to time execution',
            'CHECKPAR':'Whether to check parameters',
            'ENVIRON':'Additional environment settings',
            'CLUSTER':'Nodes available on a cluster',
            'BATCH' : 'Parameter file for batch jobs on a cluster'
            }
        rsf.conf.set_options(self,opts)

        root = self.get('RSFROOT',rsf.prog.RSFROOT)
        self.bindir = os.path.join(root,'bin')

        self.sfpen = os.path.join(self.bindir,'sfpen')
        self.pspen = os.path.join(self.bindir,'pspen')
        self.vppen = os.path.join(self.bindir,'vppen')
#        self.runonnode = os.path.join(self.bindir,'runonnode')

        self.figs = os.environ.get('RSFFIGS',os.path.join(root,'share','madagascar','figs'))

        try:
            sys.path.append(os.path.abspath(os.path.join(root,'share','madagascar','etc')))
            import config
            self.OMP = config.OMP
            self.MPI = config.MPICC
        except:
            self.OMP=False
            self.MPI=False

        cwd = os.getcwd()
        self.cwd = cwd

        # path for binary files
        self.path = rsf.path.getpath(cwd)
        tmpdatapath = os.environ.get('TMPDATAPATH',self.path)
        rsf.path.sconsign(self)

        self.resdir = resdir
        self.figdir = re.sub('.*\/((?:[^\/]+)\/(?:[^\/]+)\/(?:[^\/]+))$',
                             self.figs+'/\\1',cwd)
        self.progsuffix = self['PROGSUFFIX']

        # Keep certain environmental variables in the environment
        for env in keepenv:
            getenv = os.environ.get(env)
            if getenv:
                self.Append(ENV={env:getenv})

        self.hostname = socket.gethostname()

        # Keep environmental variables needed for SLURM
        for env in list(os.environ.keys()):
            if 'SLURM_' == env[:6] or 'TACC_' == env[:5] or '_ModuleTable' == env[:12]:
                self.Append(ENV={env:os.environ[env]})

        if sys.platform[:6] == 'cygwin':
            exe = ''
        else:
            exe = '.exe'

        libdir = os.path.join(root,'lib')
        incdir = os.path.join(root,'include')

        self.Append(ENV={'RSFROOT':root,
                         'DATAPATH':self.path,
                         'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS',rsf.node.cpus()),
                         'TMPDATAPATH': tmpdatapath,
                         'XAUTHORITY':
                         os.environ.get('XAUTHORITY',
                                        os.path.join(os.environ.get('HOME'),
                                                     '.Xauthority'))},
                    BUILDERS={'Retrieve':Retrieve,
                              'Test':Test,
                              'Echo':Echo},
                    LIBPATH=[libdir],
                    CPPPATH=[incdir],
                    F90PATH=[incdir],
                    LIBS=[libs],
                    PROGSUFFIX=exe)
        self.Prepend(LIBS=[self.get('DYNLIB','')+'rsf'])

        minesjtk = self.get('MINESJTK',None)
        usejava = self.get('JAVA_HOME',None)
        if not usejava: usejava = self.get('JAVA_SDK',None)
        if minesjtk or usejava:
            classpath = []
            classpath.append(os.path.join(libdir,'rsf.jar'))
            userclasspath = os.environ.get('CLASSPATH',None)
            if userclasspath: classpath.append(userclasspath)
            if minesjtk:
                classpath.append(minesjtk)
                self.Append(JAVACLASSPATH=':'.join(classpath))
            classpath.append('.')
            self.Append(ENV={'CLASSPATH':':'.join(classpath)})

        path = {'darwin': '/opt/local/bin',
                'irix': '/usr/freeware/bin',
                'cygwin':
                '/usr/X11R6/bin:/usr/lib/lapack:' + libdir}
        for plat in list(path.keys()):
            if sys.platform[:len(plat)] == plat:
                self['ENV']['PATH'] = ':'.join([path[plat],
                                                self['ENV']['PATH']])
        pythonpath = os.path.join(sys.prefix,'bin')
        if os.path.isdir(pythonpath):
            self['ENV']['PATH'] = ':'.join([pythonpath,
                                            self['ENV']['PATH']])

        if sys.platform[:6] == 'cygwin':
            self['ENV']['SYSTEMROOT'] = os.environ.get('SYSTEMROOT')

        self['PROGPREFIX']=''
        self.view = []
        self.prnt = []
        self.lock = []
        self.test = []
        self.coms = []
        self.data = []
        self.rest = []
        sys.path.append('../../../Recipes')

        timer = self.get('TIMER')
        if timer and timer[0] != 'n' and timer[0] != '0':
            self.timer = WhereIs('time') + ' '
        else:
            self.timer = ''

        self.mpirun = self.get('MPIRUN',WhereIs('ibrun') or WhereIs('mpirun'))

        checkpar = self.get('CHECKPAR')
        self.checkpar = checkpar and checkpar[0] != 'n' and checkpar[0] != '0'

        self.environ = self.get('ENVIRON','')

        self.batch = self.get('BATCH')

        self.jobs = GetOption('num_jobs') # getting information from scons -j

        cluster = self.get('CLUSTER',os.environ.get('RSF_CLUSTER','localhost 1'))
        hosts = cluster.split()
        self.nodes = []
        for i in range(1,len(hosts),2):
            nh = int(hosts[i])
            self.nodes.extend([hosts[i-1]]*nh)
        self.ip = 0

        # self.nodes is a list of CPUs
        # self.jobs is the number of jobs
        # self.ip is the current CPU

        for key in list(self['ENV'].keys()):
            self.environ = self.environ + " %s='%s'" %(key,self['ENV'][key])

    def __Split(self,split,reduction,
                sfiles,tfiles,flow,stdout,stdin,jobmult,suffix,prefix,src_suffix):
        '''Split jobs for pscons'''

        if self.jobs < split[1]:
            jobs = self.jobs*jobmult
            w = int(1+float(split[1])/jobs) # length of one chunk
        else:
            jobs = split[1]
            w = 1
            jobmult = 1

        par_sfiles = copy.copy(sfiles)
        par_targets = {}
        for tfile in tfiles:
            par_targets[tfile] = []

        prev_par_tfiles = []
        bigjobs = split[1] - jobs*(w-1)
        cflow = flow
        for i in range(jobs):
            if i < bigjobs:
                chunk=w
                skip=i*w
            else:
                chunk=w-1
                skip=bigjobs*w+(i-bigjobs)*chunk
            for j in split[2]:
                if 0 == j and stdin and \
                        flow.rfind('$SOURCE') < 0 and \
                        flow.rfind('${SOURCES[0]}') < 0:
                    # For the first input (if stdin), use windowing
                    # to avoid creation of a chunk file
                    par_sfiles[j] = sfiles[j]
                    cflow = '''
                    window n%d=%d f%d=%d squeeze=n icpu=%d ncpu=%d |
                    ''' % (split[0],chunk,split[0],skip,
                           i%self.jobs,self.jobs) + flow
                else:
                    source = sfiles[j] + '_' + str(i)
                    par_sfiles[j] = source

                    self.Flow(source,sfiles[j],
                              '''
                              window n%d=%d f%d=%d squeeze=n icpu=%d ncpu=%d
                              ''' % (split[0],chunk,split[0],skip,
                                     i%self.jobs,self.jobs),
                              noderotate=stdin)

            par_tfiles = []
            for j in range(len(tfiles)):
                tfile = tfiles[j]
                par_tfile = tfile + '__' + str(i)
                par_tfiles.append(par_tfile)

                par_targets[tfile].append(par_tfile)
            par_sfiles0 = copy.copy(par_sfiles)
            if i >= self.jobs:
                par_sfiles0.append (tfile + '__' + str(i % self.jobs))

            # operation on one chunk
            self.Flow(par_tfiles,par_sfiles0,cflow,
                      stdout,stdin,1,
                      suffix,prefix,src_suffix)

        # Reduce parallel TARGETS down to original TARGETS:
        for tfile in tfiles:
            self.Flow(tfile,par_targets[tfile],
                      '%s ${SOURCES[1:%d]}' % (reduction,jobs),
                      local=1)

    def Flow(self,target,source,flow,stdout=1,stdin=1,rsfflow=1,
             suffix=sfsuffix,prefix=sfprefix,src_suffix=sfsuffix,
             split=[],np=1,reduce='cat',jobmult=1,local=0,noderotate=1,
             workdir=None,wall=''):

        if not flow:
            return None

        if type(target) is list:
            tfiles = target
        else:
            tfiles = target.split()

        if source:
            if type(source) is list:
                sfiles = source
            else:
                sfiles = source.split()
        else:
            sfiles = []

        if self.hostname[-15:] == 'tacc.utexas.edu':
            mpirun = '%s tacc_affinity' % self.mpirun
        else:
            mpirun = '%s -np %s' % (self.mpirun,np)

        if split:
            if len(split) < 2:
                split.append(1)
            if len(split) < 3:
                split.append(list(range(len(sfiles))))

            if reduce.find('axis=') < 0:
                reduction = '%s axis=%d' % (reduce,split[0])
            else:
                reduction = reduce

            if split[1] == 'omp' or split[1] == 'mpi':
                if (split[1] == 'omp' and self.OMP) or \
                   (split[1] == 'mpi' and self.MPI):
                    splitpar = 'split=%d ' % split[0]
                    if reduce == 'add':
                        splitpar += ' join=0'
                    else:
                        join = re.search('cat\s+axis=(\d)',reduce)
                        if join:
                            splitpar += ' join=%s' % join.group(1)
                    flow = '|'.join([' '.join([split[1],splitpar,x]) for x in flow.split('|')])
                    for k in split[2]:
                        # par=${SOURCES[k]} -> _par=${SOURCES[k]}
                        flow = re.sub(r'(\S+=\${SOURCES\[%d\]})' % k,'_\\1',flow)
            elif self.jobs > 1 and rsfflow and sfiles:
                # Split the flow into parallel flows
                self.__Split(split,reduction,
                             sfiles,tfiles,flow,stdout,stdin,
                             jobmult,suffix,prefix,src_suffix)
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
            if noderotate:
                self.ip = self.ip + 1
                if self.ip == len(self.nodes):
                    self.ip = 0

        if node != 'localhost':
            remote = '%s $( %s $)' % (WhereIs('env'),self.environ)
        else:
            remote = ''

        command = rsf.flow.Flow(sources,flow,self.bindir,rsfflow,
                                self.checkpar,self.coms,prefix,self.progsuffix,
                                remote,stdout,stdin,self.timer,mpirun,workdir,
                                self.batch,np,wall)

        # May need to do it remotely
        if remote:
            command = re.sub('"','\\"',command)
            command = ' '.join(['$( ssh',node,'$) \"cd ',self.cwd,';',command,'\"'])

        targets = []
        for file in tfiles:
            if (not re.search(suffix + '$',file)) and ('.' not in file):
                file = file + suffix
            targets.append(file)
        if suffix == sfsuffix and re.search('/',targets[0]) and 1 == stdout:
            subdir = os.path.dirname(os.path.join(self.path,targets[0]))
            rsf.path.mkdir(subdir)
            command = command + ' --out=%s' % os.path.join(self.path,'${TARGET}@')

        if workdir:
            command = re.sub(r'\$(SOURCE|TARGET)',r'${\1.abspath}',command)
            command = re.sub(r'\$\{(SOURCES|TARGETS)(\[[^\]]+\])\}',r'${\1\2.abspath}',command)

        flow = self.Command(targets,sources,command)

        if workdir:
            Clean(flow,workdir)

        if suffix == sfsuffix:
            binaries = list(map(lambda x, self=self: self.path + x + '@',
                           list(filter(lambda x, suffix=suffix:
                                      x[-len(suffix):] == suffix,targets))))
            if binaries:
                Clean(flow,binaries)

        self.Default(flow)
        return flow

    def Plot (self,target,source,flow=None,suffix=vpsuffix,vppen=None,
              view=None,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        if 'Annotate'==flow:
            if not type(source) is list:
                source = source.split()
            flow = os.path.join(self.bindir,'vpannotate') + \
              ' text=${SOURCES[1]} batch=y ${SOURCES[0]} $TARGET'
            kw.update({'src_suffix':vpsuffix,'stdin':0,'stdout':-1})
        elif flow in combine:
            if not type(source) is list:
                source = source.split()
            flow = combine[flow](*[self.vppen,len(source)])
            if vppen:
                flow = flow + ' ' + vppen
            kw.update({'src_suffix':vpsuffix,'stdin':0})
        if view:
            if suffix==vpsuffix and not 'matplotlib' in flow:
                flow = flow + ' | %s pixmaps=y' % self.sfpen
            kw.update({'stdout':-1})
        kw.update({'suffix':suffix})
        return self.Flow(*(target,source,flow), **kw)
    def Result(self,target,source,flow=None,suffix=vpsuffix,**kw):
        if not flow: # two arguments
            flow = source
            source = target
        target2 = os.path.join(self.resdir,target)
        if 'matplotlib' in flow:
            pngflow = flow + ' format=png'
            pngsuffix = '.png'
            kw.update({'suffix':pngsuffix})
            pngplot = self.Plot(*(target,source,pngflow), **kw)

            flow += ' format=pdf'
            suffix = '.pdf'

        kw.update({'suffix':suffix})
        plot = self.Plot(*(target2,source,flow), **kw)
        target2 = target2 + suffix
        if suffix == vpsuffix:
            viewer = self.sfpen
        elif suffix == '.pdf':
            viewer = WhereIs('acroread') or WhereIs('kpdf') \
              or WhereIs('evince') or WhereIs('xpdf') or WhereIs('gv') \
              or WhereIs('open')
        elif suffix == '.eps':
            viewer = WhereIs('evince') or WhereIs('gv') or WhereIs('open')
        else:
            viewer = None

        if viewer:
            view = self.Command(target + '.view',plot,viewer + " $SOURCES",
                                src_suffix=suffix)
            self.view.append(view)

        prnt = self.Command(target + '.print',plot,
                            self.pspen + " printer=%s $SOURCES" % printer,
                            src_suffix=vpsuffix)
        self.prnt.append(prnt)

        locked = os.path.join(self.figdir,target+suffix)
        self.InstallAs(locked,target2)
        self.Alias(target + '.lock',locked)
        self.lock.append(locked)

        if suffix == vpsuffix:
            self.Command(target + '.flip',target2,
                        '%s $SOURCE %s' % (self.sfpen,locked))
        test = self.Test('.test_'+target,target2,
                         figdir=self.figs,bindir=self.bindir)
        self.test.append(test)
        self.Alias(target + '.test',test)
        self.rest.append(target)
        return plot
    def End(self):
        self.Command('.rsfproj',self.lock,action=Action(self.Info))
        self.Echo('results',None,err=self.rest)
        if self.view: # if any results
            self.Alias('view',self.view)
            self.Alias('print',self.prnt)
            self.Alias('lock',self.lock+['.rsfproj'])
            self.Alias('test',self.test)
        else:
            self.Echo('test',None,err='Nothing to test')
    def Info(self,target=None,source=None,env=None):
        'Store project information in a file'
        infofile = str(target[0])
        info = open(infofile,'w')
        info.write('uses=' + str(self.coms) + '\n')
        info.write('data=' + str(self.data) + '\n')
        info.close()
        sizes = os.path.join(self.get('RSFROOT',rsf.prog.RSFROOT),'bin','sfsizes')
        su = env.get('su',0)
        suffix = env.get('suffix',sfsuffix)
        os.system('%s files=n su=%d *%s >> %s' % (sizes,su,suffix,infofile))
        return 0
    def Fetch(self,files,dir,private=None,server=dataserver,top='data',usedatapath=True):
        global dataserver
        if private:
            self.data.append('PRIVATE')
        elif server=='local':
            self.data.append('LOCAL')
        else:
            if not type(files) is list:
                files = files.split()
            if server == None:
                dataserver = get_dataserver()
                server = dataserver
            for fil in files:
                if server != dataserver:
                    if top:
                        self.data.append(os.path.join(server,top,dir,fil))
                    else:
                        self.data.append(os.path.join(server,dir,fil))
                elif top:
                    self.data.append(os.path.join(top,dir,fil))
                else:
                    self.data.append(os.path.join(dir,fil))
        return self.Retrieve(files,None,
                             dir=dir,private=private,
                             top=top,server=server,usedatapath=usedatapath)

# Default project
project = Project()

def Flow(target,source,flow,**kw):
    return project.Flow(*(target,source,flow), **kw)
def Plot (target,source,flow=None,**kw):
    return project.Plot(*(target,source,flow), **kw)
def Result(target,source,flow=None,**kw):
    return project.Result(*(target,source,flow), **kw)
def Fetch(file,dir,private=0,**kw):
    return project.Fetch(*(file,dir,private), **kw)
def End(**kw):
    return project.End(*[], **kw)
def Program(*arg,**kw):
    return project.Program(*arg, **kw)
def Get(name):
    return project['ENV'].get(name)

if __name__ == "__main__":
    import pydoc
    pydoc.help(Project)
