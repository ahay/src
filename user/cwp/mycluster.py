'''
This python module overrides the standard Madagascar
python SCons methods.  We do this, carefully so that 
the SCons methods still work unless we set the Global 
variable CSCONS.  If CSCONS is set, then SCons does 
not work, and rather we generate a series of PBS
job wrappers for our commands.

This file is customized in a number of areas.  To quickly
find all changes, search for the tag:  CWPONLY

The python executable cscons automatically turns on 
the PBS capabilities of this module.  See the documentation
online for Mio (on the iTeam website) to learn how to 
configure your SConstructs to use this module.

author:  Jeff Godwin, Jan 12 2011, CWP

Some new functions are added to make cscons working in an 
iterative inversion scheme. These functions include reading
user input from command line, setup different pbs directory
name at each iteration, executing an command line input, etc.

author: Tonging Yang, Jun 16 2011, CWP

'''
import rsf.proj, os, string, subprocess, sys
import SCons
import SCons.Script
import SCons.Script.SConsOptions as SConsOptions

# Default project
project = rsf.proj.Project()


###############
######################
# CWPONLY - From here down, the rest of this file is basically custom. 
###############
JOB_NAME = None
EMAIL    = None
DEFAULT_TIME= None
DEFAULT_PPN = None
NODETYPE = None
CLUSTER = None
AUTO_SUBMIT = False
SAVED = []

JOB_QUEUE    = []
CURRENT_JOB  = None

pbs_dirt = ''

CSCONS  = False
CLUSTER_CONFIGURED = False

SCONSIGNS = []

def NamePbsDirectory(name):
    global pbs_dirt

    pbs_dirt = name

class FakeOptionParser(object):
    """
    A do-nothing option parser, used for the initial OptionsParser variable.

    During normal SCons operation, the OptionsParser is created right
    away by the main() function.  Certain tests scripts however, can
    introspect on different Tool modules, the initialization of which
    can try to add a new, local option to an otherwise uninitialized
    OptionsParser object.  This allows that introspection to happen
    without blowing up.

    """
    class FakeOptionValues(object):
        def __getattr__(self, attr):
            return None
    values = FakeOptionValues()
    def add_local_option(self, *args, **kw):
        pass

OptionsParser = FakeOptionParser()

def OptionInitial():

    global OptionsParser

    parser = SConsOptions.Parser('')
    values = SConsOptions.SConsValues(parser.get_default_values())

    OptionsParser = parser

    sconsflags = os.environ.get('SCONSFLAGS', '')
    all_args = sconsflags.split() + sys.argv[1:]

    options, args = parser.parse_args(all_args, values)
    options = parser.values

def AddOption(*args, **kw):

    if 'default' not in kw:
        kw['default'] = None
    result = OptionsParser.add_local_option(*args, **kw)

    return result

def GetOption(name):

    return getattr(OptionsParser.values, name)

def createPBSfile(name,email,nodes,ppn,time,last,next,tasks,relaunch,run,content,nodetype=None,parallel=False):
    '''
    Where we actually create our pbs job files. Modify this to change this to fit
    your cluster, or changes to Mio and Ra.
    '''

    global pbs_dirt

    lines = []
    lines.append('#!/bin/bash')

    if nodetype:
        lines.append('#PBS -l nodes=%d:ppn=%d:%s' % (nodes,ppn,nodetype))
    else:
        lines.append('#PBS -l nodes=%d:ppn=%d' % (nodes,ppn))
    lines.append('#PBS -l naccesspolicy=singlejob')
    lines.append('#PBS -l walltime=%d:00:00' % time)
    lines.append('#PBS -N %s' % name)
    lines.append('#PBS -o %s/%s.out' % (pbs_dirt,name))
    lines.append('#PBS -e %s/%s.err' % (pbs_dirt,name))
    lines.append('#PBS -V')
    if email:
        lines.append('#PBS -m a')
        lines.append('#PBS -M %s' % email)
    lines.append('#-----------')
#    lines.append('export DATAPATH=%s' % rsf.path.getpath(os.getcwd()))
    lines.append('export SCONS_OVERRIDE=1')
    lines.append('cd %s' % os.getcwd())
    lines.append('set -e # Catch all errors immediately and exit')
#    lines.append('echo "JOB: %s running" >> pbs/jobs.txt' % name)

    if last:
        depends = string.join(['%s' % item for item in last])
#        lines.append("""for i in %s; do  if [ -n "$(grep "$i done" pbs/jobs.txt)" ]; then echo "found $i"; elif [ -n "$(grep "$i running" pbs/jobs.txt)" ]; then echo "running $i"; if [ `grep "$i check" pbs/jobs.txt | wc -l` -ge 5 ]; then echo "already checked $i 5 times"; exit; else echo "JOB: $i check" >> pbs/jobs.txt; sleep 300; echo "Resubmitting myself"; ssh $PBS_O_HOST "source ~/.bash_profile; cd $PBS_O_WORKDIR; qsub pbs/%s" ; exit; fi; elif [ -n "$(grep "$i error" pbs/jobs.txt)" ]; then echo "error $i"; exit; else echo "did not find $i"; exit; fi; done""" % (depends, name))

    if parallel:
        lines.append('sort -u $PBS_NODEFILE > %s/%s-nodes' % (pbs_dirt,name))
        lines.append('match %s/%s-nodes %s/%s-shell > %s/%s-app ' % (pbs_dirt,name,pbs_dirt,name,pbs_dirt,name))
        lines.append('mpiexec -app %s/%s-app' % (pbs_dirt,name))
    else:
        if relaunch or run:
            lines.append("""
            ssh $PBS_O_HOST "source ~/.bash_profile; cd $PBS_O_WORKDIR; %s"
            """ %content)
        else:
            lines.append('scons -f %s/SConstruct-%s' % (pbs_dirt,name))

        #lines.append(tasks+'\n')

#    lines.append('echo "JOB: %s done" >> pbs/jobs.txt' % name)

    if next == None:
        if not relaunch:
            global SCONSIGNS
            lines.append('sfdbmerge outdb=%s %s ' % (project.path+'.sconsign.dbhash', ' '.join(SCONSIGNS)))

    file = open('%s/%s' % (pbs_dirt,name),'w')
    text = string.join(lines,'\n')
    if CLUSTER=='ra':
        text.replace('qsub','msub')
    file.write(text)
    file.close()


class Job:
    def __init__(self,time=0,nodes=1,ppn=0,notes=''):
        if not time or not nodes or not ppn:
            raise Exception('Must specify, time, nodes and ppn, got %f %d %d' % (time,nodes,ppn))
        self.name  = ''
        self.time  = time
        if not time: self.time = DEFAULT_TIME
        self.nodes = nodes
        self.ppn   = ppn
        if not ppn: self.ppn = DEFAULT_PPN
        self.tasks = []
        self.notes = notes
        self.all = []
        self.last = None
        self.next = None
        self.sconsign = []
        self.content = ''
        self.relaunch = False
        self.run = False

    def make(self):
        createPBSfile(self.name,EMAIL,self.nodes,self.ppn,
            self.time,self.last,self.next,'\n'.join(self.tasks),self.relaunch,self.run,self.content,nodetype=NODETYPE)

    def keep(self):
        if len(self.tasks) > 0:
            return True
        else:
            return False

    def setLast(self,last):
        self.last = last
    def setNext(self,next):
        self.next = next

    def prepare(self):
        global pbs_dirt
        if self.name:
            sconsign = '%s.sconsign.dbhash' % self.name
            self.sconsign.append(sconsign)
            self.all = [self.name]
            file = open('%s/SConstruct-%s' % (pbs_dirt,self.name),'w')
            file.write('from rsf.proj import *\n')
            file.write( \
'''
import dbhash 	 
proj = Project()
proj.SConsignFile("%s",dbhash)	 
''' % (sconsign))
            file.write(string.join(self.tasks,'\n'))
            file.write('\n')
            file.write('End()\n')
            file.close()

    def add(self, command):
        self.tasks.append(command)

    def store(self,content):
        self.content = content

    def restart(self):
        self.all = [self.name]

    def __str__(self):
        return 'Serial:%d %d %d %s' % (self.time,self.nodes,self.ppn,self.notes)

class Parallel(Job):
    def __init__(self,time=0, nodes=0, ppn=0,ipn=1):
        Job.__init__(self,time, nodes, ppn)
        self.ipn = ipn
        self.scripts = []
        self.nnl = self.nodes  # Number of nodes for the last job

    def make(self):
        for i in range(len(self.all)):
            next = None
            nodes = self.nodes
            if i == len(self.all)-1:
                next = self.next
                nodes = self.nnl
            createPBSfile(self.all[i],EMAIL,nodes,self.ppn,
                self.time,self.last,next,'',False,False,'',nodetype=NODETYPE,parallel=True)

    def keep(self):
        if len(self.scripts) > 0:
            return True
        else:
            return False

    def flush(self):
        if len(self.tasks) > 0: 
            text = '\n'.join(self.tasks)
            self.scripts.append(text)
            self.tasks = []

    def prepare(self):
        global pbs_dirt
        if self.name:
            nscripts = len(self.scripts)
           
            inode = 0
            ijob  = 0
            j = 1

            names = []
            nodescripts = []
            for i in range(len(self.scripts)):
                nodescripts.append(self.scripts[i])

                if j == self.ipn or i == len(self.scripts) -1:
                    scriptname = self.name+'-%02d-%02d' % (ijob,inode)
                    names.append(scriptname)

                    file = open('%s/SConstruct-%s' % (pbs_dirt,scriptname),'w')
                    file.write('from rsf.proj import *\n')
                    sconsign = '%s.sconsign.dbhash' % scriptname
                    file.write( \
'''
import dbhash 	 
proj = Project()
proj.SConsignFile("%s",dbhash)	 
''' % (sconsign))

                    self.sconsign.append(sconsign)

                    for tscript in nodescripts:
                        file.write(tscript)
                        file.write('\n')
                    
                    file.write('End()\n')
                    file.close()
                    
                    inode += 1
                    j = 1
                    nodescripts = []
                else: 
                    j += 1

                if inode == self.nodes or i == len(self.scripts)-1:
                    jobname = self.name+'-%02d' % ijob
                    file = open('%s/%s-shell'%(pbs_dirt,jobname),'w')
                    file.write(string.join(['scons -f %s/SConstruct-%s'% (pbs_dirt,name) for name in names],'\n'))
                    file.close()
                    ijob += 1
                    if i == len(self.scripts)-1 and inode != self.nodes:
                        self.nnl = inode
                    inode = 0
                    names = []
                    self.all.append(jobname)

            
    def __str__(self):
        return 'Parallel: %d %d %d %d'  % (self.time,self.nodes,self.ppn,self.ipn)

def Cluster(name, email=None, cluster='mio', time=1,ppn=8,nodetype=None,submit=False):
    global DEFAULT_TIME, DEFAULT_PPN, JOB_NAME,EMAIL,CLUSTER,NODETYPE,CLUSTER_CONFIGURED
    DEFAULT_TIME = time
    DEFAULT_PPN  = ppn
    JOB_NAME = name
    EMAIL = email
    CLUSTER = cluster
    NODETYPE = nodetype
    AUTO_SUBMIT =submit 

    CLUSTER_CONFIGURED = True

def Serial(time=None,ppn=None,tasks=None):
    global CURRENT_JOB, JOB_QUEUE, CLUSTER_CONFIGURED
    if not CLUSTER_CONFIGURED:
        print "FATAL ERROR: You must initialize default parameters first, by calling the Cluster function"
        sys.exit(155)

    if not time: time = DEFAULT_TIME
    if not ppn: ppn = DEFAULT_PPN
    CURRENT_JOB = Job(time,1,ppn)
    JOB_QUEUE.append(CURRENT_JOB)
    if tasks:
        for task in tasks:
            CURRENT_JOB.add(task)

        Serial()

def Fork(time=0,nodes=0,ipn=0,ppn=0):
    global CURRENT_JOB, JOB_QUEUE, CLUSTER_CONFIGURED
    if not CLUSTER_CONFIGURED:
        print "FATAL ERROR: You must initialize default parameters first, by calling the Cluster function"
        sys.exit(155)
    if not time or not nodes or not ipn:
        raise Exception('Must enter a value for nodes, time, and ipn')
    if not ppn: ppn = DEFAULT_PPN
    CURRENT_JOB = Parallel(time=time,nodes=nodes,ipn=ipn,ppn=ppn)
    JOB_QUEUE.append(CURRENT_JOB)

def Iterate():
    global CURRENT_JOB
    try:
        CURRENT_JOB.flush()
    except:
        print 'Could not find a valid Job type.  Did you forget to Fork?'

def Join():
    global CURRENT_JOB
    CURRENT_JOB.flush()
    CURRENT_JOB = None

def __createstring(rsfcommand,target,source,command,kw):

    outsou = ''
    if type(source) == str:
        outsou = '"%s"' % source 
    else:
        outsou = str(source)

    if type(target) == str:
        outtar = '"%s"' % target
    else:
        outtar = str(target)

    if source == None and rsfcommand != 'Flow':
        command = """%s(%s,\n\t'''%s'''"""% (rsfcommand,outtar,command)
    else:
        command = """%s(%s,%s,\n\t'''%s'''""" % (rsfcommand,outtar,outsou,command)

    for key in kw.keys():
        if type(kw[key]) == str:
            command+=',%s="%s"' % (str(key),str(kw[key]))
        else:
            command+=',%s=%s' % (str(key),str(kw[key]))

    command += ')'
    return command


def Flow(target,source,flow,**kw):
    global CURRENT_JOB
    if CSCONS:

        if not CLUSTER_CONFIGURED:
            raise Exception("You must call Cluster() to initialize parameters")

        if not CURRENT_JOB: # No job
            Serial()

        mpi = None
        nodes = -1
        if kw.has_key('mpi'):
            mpi = kw['mpi']
            kw.pop('mpi')
            if not kw.has_key('nodes'):
                raise Exception('Must specify # of nodes for mpi jobs')
            else:
                nodes = kw.pop('nodes')

            if kw.has_key('ppn'):
                ppn = kw.pop('ppn')
                if not ppn: ppn= DEFAULT_PPN
            else:
                ppn = DEFAULT_PPN

            if not kw.has_key('time'):
                raise Exception('you must specify time for mpi jobs')
            else:
                time = kw.pop('time')
                if not time: raise Exception("must specify time for mpi jobs")

            if not kw.has_key('np'):
                if nodes <= 0:
                #try:
                #    kw['np'] = nodes
                    raise Exception('You must specify the number of nodes for mpi jobs')
            else:
                kw.pop('np')

                #if not kw['np']: 
                #    kw['np'] = nodes


#        command = project.Flow(target,source,flow,**kw)

        if mpi and nodes > 0:
            CURRENT_JOB = Job(nodes=nodes,time=time,ppn=ppn,notes='MPI')
            JOB_QUEUE.append(CURRENT_JOB)
       
        command = __createstring('Flow',target,source,flow,kw)
        CURRENT_JOB.add(command)

    else:
        if kw.has_key('time'): 
            kw.pop('time')
        if kw.has_key('np'):
            kw['np'] = -1 # We are running locally
            kw.pop('np')
        if kw.has_key('mpi'):
            kw.pop('mpi')
        if kw.has_key('nodes'):
            kw.pop('nodes')
        if kw.has_key('ppn'):
            kw.pop('ppn')
        return apply(project.Flow,(target,source,flow),kw)

def Plot (target,source,flow=None,**kw):
    if CSCONS:
        global CURRENT_JOB, JOB_QUEUE
        if not CURRENT_JOB: # No job
            Serial()

#        command = project.Plot(target,source,flow,**kw)
        if flow == None:
            flow = source
            source = None
        command = __createstring('Plot',target,source,flow,kw)
        CURRENT_JOB.add(command)
    else:
        return apply(project.Plot,(target,source,flow),kw)

def Result(target,source,flow=None,**kw):
    if CSCONS:
        global CURRENT_JOB, JOB_QUEUE
        if not CURRENT_JOB: Serial()
        if flow == None:
            flow = source
            source = None
        command = __createstring('Result',target,source,flow,kw)
        CURRENT_JOB.add(command)
    else:
        return apply(project.Result,(target,source,flow),kw)

def Force(content):
    global CURRENT_JOB

    CURRENT_JOB = None
    Serial()

    CURRENT_JOB.add(content)
    CURRENT_JOB.store(content)
    CURRENT_JOB.relaunch = True

def Run(content):
    global CURRENT_JOB

    CURRENT_JOB = None
    Serial()

    CURRENT_JOB.add(content)
    CURRENT_JOB.store(content)
    CURRENT_JOB.relaunch = False 
    CURRENT_JOB.run = True

    CURRENT_JOB = None

####
# No changes here
def Fetch(file,dir,private=0,**kw):
    return apply(project.Fetch,(file,dir,private),kw)
#####

def Save(file):
    Flow(file+'.hh',file,'window out=stdout')

def End(**kw):
    global pbs_dirt

    if CSCONS and CLUSTER_CONFIGURED:
        if len(JOB_QUEUE) > 0:
            print "Removing extra jobs..."
            final_queue = []
            # We could possibly have jobs with no tasks
            # Remove these from the queue
            for job in JOB_QUEUE:
                if job.keep():
                    final_queue.append(job)

            print "Found: %d jobs" % len(final_queue)
            print 'Creating job directory...'
            if not os.path.exists(pbs_dirt):
                os.mkdir(pbs_dirt)
            else:
                if os.path.isdir(pbs_dirt):
                    print '......pbs directory already exists'
                else:
                    raise Exception('pbs directory exists, but is not suitable for script files?')

            if not os.path.exists('Fig'):
                os.mkdir('Fig')
            else:
                if os.path.isdir('Fig'):
                    print '......Fig directory already exists'
                else:
                    raise Exception('Fig directory exists but is not suitable for vpl files?')

            print 'Preparing jobs...'
            
            global SCONSIGNS

            i = 0; n = len(final_queue);
            for job in final_queue:
                if isinstance(job,Parallel):
                    job.name = JOB_NAME+'-l-%02d' % i
                else:
                    job.name = JOB_NAME+'-s-%02d' % i
                print 'Job %d/%d - %s - %s' % (i,n,str(job),job.name)
                if not job.relaunch:
                  job.prepare()
                  SCONSIGNS.extend(job.sconsign)
                else:
                  job.restart()

                i += 1

            SCONSIGNS = map(lambda x: pbs_dirt + '/' + x, SCONSIGNS)

            print 'Making job files...'
            i = 0
            for job in final_queue:
                if i > 0:
                    last = final_queue[i-1].all
                    job.setLast(last)
                if i < len(final_queue)-1:
                    next = final_queue[i+1].all
                    job.setNext(next)

                job.make()
                i += 1

            print 'Submitting jobs...\n'

            pbs_names = []

            j = 0
            for job in final_queue:
                subjobs = []
                for pbs in job.all:
                    command = 'qsub %s/' %pbs_dirt
                    if CLUSTER == 'ra':
                        command = command.replace('qsub','msub')
                    command += pbs
                    if job.last != None:
                        command +=' -W depend=afterok:%s' % (':'.join(pbs_names[j-1]))
                    print 'Executing...',command
                    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
                    stdout,stderr  = process.communicate()
                    if len(stdout) < 3:
                        print 'WARNING: submission was not successful?',stdout
                    print 'Job submitted to: %s' % stdout
                    subjobs.append(stdout.strip('\n'))

                pbs_names.append(subjobs)
                j += 1

        else:
            raise Exception("Did not find any jobs?")
            
    return apply(project.End,[],kw)
