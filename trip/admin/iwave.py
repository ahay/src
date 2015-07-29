####### IWAVE/IWAVE++ numerical experiment SConstruct template ########
#
# Stores results (Flow targets) in working directory, from two 
# kinds of flows:
# 
# I. flows invoking non-IWAVE commands: these should be written 
# in normal Madagascar fashion, and placed in the 'local commands'
# section of this SConstruct.
# 

# II. IWAVE flows: defined in the flow loop at the bottom of this
# file, which the user should not alter. This utility assumes that
# IWAVE/IWAVE++ commands will read their input from parameter files,
# which should reside in this working directory. A parameter file
# is identified by the suffix ".par". 

# Note that other commands which read all input data from a single
# parameter file with suffix ".par" should work the same way - will
# refer to all such commands as "IWAVE commands".

# For parfile 'foo.par', 'foo' is a key in the jobs dictionary. an
# entry in jobs looks like

# 'foo' : { 'cmd' : <command name>, 
#           'dep' : <dependency list>, 
#           'tgt' : <target dictionary>,
#           'exe' : (optional) <execution environment dictionary> }
# 
# That is, the root parfile name serves as the keyword, and the
# associated value is another dictionary. The four keywords of 
# the value dictionary and their associated values are
#
# 'cmd' : path to the IWAVE/IWAVE++ command 
# 'dep' : results on which this job depends - each result is
#         a file, which is a target of a Flow defined in this 
#         SConstruct, and resides in this directory
# 'tgt' : another dictionary, in which each keyword is a 
#         keyword in the parfile specifying an output file, and 
#         the corresponding value is the filename under 
#         which the result should be archived in this 
#         working directory. Each such file is a Flow
#         target. So for example if foo[tgt] contains the pair 
#         'datafile : mydata.su', then the flow for target 'mydata.su'
#         moves the file associated to keyword 'datafile' in the 
#         job parfile to the file mydata.su from whereever it currently
#         resides (usually a scratch directory for storing various process
#         outputs from another flow) to "mydata.su" in this working 
#         directory.

# 'exe' : (optional) yet another dictionary, specifying execution
#         environment.  If given, the keyword 'platf', specifies
#         execution platform, which falls into one of three major
#         categories:

#         SERIAL - If either 'exe' is not a key, or if the value of
#         'platf' is '', execution takes place in the foreground, by
#         serial (non-mpi) command line invocation of the command.

#         COMMAND-LINE MPI: A second predefined value of 'platf' is
#         'mpi', for command-line mpi execution: for this parallel
#         mode, you must specify integer values for the keyword 'ppn'
#         (processes per node - by convention, this form of
#         parallelism takes place on a single "node").

#         Thus foo:exe defined as
#                     'exe'     : { 'platf' : 'mpi',
#                                   'ppn'   : '16',
#                                 }
#         with foo:cmd = <cmd> results in the foreground execution of 
#             mpirun -np 16 <cmd> par=foo.par

#         BATCH: Other 'platf' values should be descriptive names for
#         environments in which the user intends to run iwave jobs
#         under batch control. Additional parameters for batch
#         submissions are:

#             wall  = Wallclock time limits (hr:mn:se)
#             nodes = number of nodes
#             ppn   = processes/cores/threads per node

#         Each batch environment name must be a keyword in the 'penv'
#         dictionary appearing just before the jobs dictionary. 'penv'
#         defines several standard inputs needed by all batch
#         systems. Here is an example:

#             penv = { (...)
#                     'stampede' : { 'batch' : 'slurm',
#                                    'queue' : 'normal',
#                                    'acct'  : 'FDTD3D-Cont',
#        	                     'mail'  : 'symes@caam.rice.edu',
#                                  }
#                      (...)
#                    }
#         The value of penv:platf:batch must be a keyword in the
#         bsys dictionary, defined in the common code section below
#         bsys defines various implementation-indedpendent attributes
#         of batch systems - see the code below for the currently 
#         accommodated systems, and feel free to add another and 
#         share!
#
#         So foo:exe defined as follows:
#                     'exe'     : { 'platf' : 'stampede',
#                                   'wall'  : '00:59:00',
#                                   'nodes' : '4',
#                                   'ppn'   : '16',
#                                  }

#         would cause the a batch file foo.bat to be written, and
#         submitted for execution on 4 nodes, 16 cores per node, of
#         environment 'stampede' via invocation of 'sbatch'

#         Thus 'bsys:batch' specifies attributes common to all
#         implementations of the batch system 'batch', 'penv:platform'
#         specifies batch parameters common to all uses of 'platform'
#         (including a choice of batch system), and 'jobs:jobname:exe'
#         specifies parameters particular to the job 'jobname',
#         including the platform on which it is to be executed.
 
# IWAVE commands execute in scratch (sub)directories generated by the flows that 
# invoke them; these flows generate the work directory names from the parfiles. 
# For example, the parfile foo.par generates a commond which executes in 
# foo.work. The subdirectories appear in WORKPATH, specified in the 
# common definitions section immediately following


import os

THISPATH        = os.getcwd()

def getFlowSignature(job, jobdict, envdict, workpath):

    # test for command definition
    if 'cmd' not in jobdict[job].keys():
        print 'iwave.signature: job dictionary ' + job + ' does not include key = cmd,'
        print 'so command cannot be determined'
        return
    
    # path to build directory
    workdir = os.path.join(workpath, job + '.work')

    # create command to put abs pathnames in par file
    parcmd = getParamCommand(job,jobdict)
    if parcmd == None:
        print 'iwave.signature: job dictionary ' + job + ' does not include key=dep'
        print 'and/or key=tgt: cannot determine dependencies'
        return
    
    # batch execution
    if isBatch(job,jobdict,envdict):
        # detect batch system attributes
        bsys = getBatchAttributes(envdict[jobdict[job]['exe']['platf']]['batch'])
        if bsys == None:
            print 'iwave.signature: failed to parse batch submission system'
            return
        
        # create batch file names
        thisbat = job + '.' + bsys['suff']

        # write batch files
        writeScript(envdict[jobdict[job]['exe']['platf']]['batch'],
                    job,
                    workdir,	
                    envdict[jobdict[job]['exe']['platf']]['queue'],
                    envdict[jobdict[job]['exe']['platf']]['mail'],
                    envdict[jobdict[job]['exe']['platf']]['acct'],
                    jobdict[job]['exe']['wall'],
                    jobdict[job]['exe']['nodes'],
                    jobdict[job]['exe']['ppn'],
                    jobdict[job]['cmd'])
        
        workcmd = 'cd ' + workdir + \
            '; ' + parcmd + \
            '; /bin/mv ' + os.path.join(THISPATH,thisbat) + \
            ' ./' + thisbat  + '; /bin/rm -f jobid; ' + \
            bsys['bshl'] + ' -c "' + \
            bsys['bcmd'] + ' ' + thisbat + '"' \
            '; (while [ ! -f ./jobid ]; do (sleep 1); done)'

    elif isMPI(job,jobdict):

        MPIROOT = os.getenv('MPIROOT')
        mpirun = os.path.join(MPIROOT,'bin/mpirun')
        # this should be checked for existance
        workcmd = 'cd '          + workdir + \
            '; ' + parcmd + \
            '; ' + mpirun + ' -np ' + jobdict[job]['exe']['ppn'] + \
            ' ' + jobdict[job]['cmd'] + ' par=./parfile'
        
    else:      

        workcmd = 'cd '          + workdir + \
            '; ' + parcmd + \
            '; ' + jobdict[job]['cmd'] + ' par=./parfile'

    # flow target list = target list from dictionary 
    worktgt = []
    for k in jobdict[job]['tgt'].keys():
        worktgt = worktgt + [os.path.join(THISPATH,jobdict[job]['tgt'][k])]

    # flow source list = dependency list from dictionary plus par file plus workdir
    worksrc = jobdict[job]['dep'] + [os.path.join(THISPATH, job + '.par'), workdir]

    if len(worktgt) == 0 or len(worksrc) == 0:
        return
    
    workdict = { 'tgt' : worktgt,
                 'src' : worksrc,
                 'cmd' : workcmd,
                 'dir' : workdir,
                 'wcd' : 'if [ ! -d ' + workdir + ' ] ; then ( /bin/mkdir ' + workdir + ' ); fi'
               }
    return workdict

def getBatchImmediateAncestor(job,jobdict,envdict):
    # extract immediate ancestor from amongst batch jobs on which this one depends
    bdep=[]
    # first develop list of batch targets which are dependencies of this:
    for j in jobdict[job]['dep']:
        for k in jobdict.keys():
            if isBatch(k,jobdict,envdict):
                for kk in jobdict[k]['tgt'].keys():
                    if j==jobdict[k]['tgt'][kk]:
                        bdep = bdep + [k]
                        
    # then eliminate any which precede others in dependence:                
    ldep=[]
    if (len(bdep) > 0):
        for k in bdep:
            f=True
            if (len(bdep) > 1):
                bdepk = bdep.remove(k)
                if bdepk == None:
                    f = True
                else:
                    for kk in jobdict[k]['tgt'].keys():
                        for l in bdepk:
                            for ll in jobdict[l]['dep']:
                                if ll==jobdict[k]['tgt'][kk]:
                                    f=False
            if f:
                ldep = ldep + [k]
    return ldep

def writeScript(batch,job,path,queue,mail,acct,wall,nodes,ppn,exe):
    if batch == 'pbs':
        f = open(job +'.pbs','w')
        f.write('#PBS -N ' + job + '\n')
	f.write('#PBS -A ' + acct + '\n')
        f.write('#PBS -V\n')
        f.write('#PBS -l walltime=' + wall + '\n')
        f.write('#PBS -l nodes=' + nodes + ':ppn=' + ppn + '\n')
        f.write('#PBS -q '+queue + '\n')
        f.write('#PBS -m abe\n')
        f.write('#PBS -M '+ mail + '\n')
        f.write('#PBS -o '+path+'/cout.txt\n')
        f.write('#PBS -e '+path+'/cerr.txt\n')
        f.write('cd '+path+'\n')
        f.write('mpiexec ' + exe + ' par=./parfile\n')
        f.write('echo $PBS_JOBID > jobid\n')
        f.close()

    if batch == 'slurm':

        f = open(job +'.bat','w')
	f.write('#!/bin/bash\n')
        f.write('#SBATCH -J ' + job + '\n')
        f.write('#SBATCH -A ' + acct + '\n')
        f.write('#SBATCH -t '+wall + '\n')
        f.write('#SBATCH -n ' + str(int(ppn) * int(nodes)) + '\n')
        f.write('#SBATCH -p '+queue + '\n')
        f.write('#SBATCH --mail-type=begin\n')
        f.write('#SBATCH --mail-type=end\n')
        f.write('#SBATCH --mail-user='+mail + '\n')
        f.write('#SBATCH --uid=' + os.getenv('USER') + '\n')
        f.write('#SBATCH --get-user-env\n')
        f.write('#SBATCH -o '+path+'/cout.txt\n')
        f.write('#SBATCH -e '+path+'/cerr.txt\n')
        f.write('cd '+path+'\n')
        f.write('ibrun ' + exe + ' par=./parfile\n')
        f.write('echo $SLURM_JOBID > jobid\n')
        f.close()

def isBatch(job,jobdict,envdict):
    if ('exe' in jobdict[job].keys()):
        # check spec of platform
        if ('platf' in jobdict[job]['exe'].keys() and \
            'wall' in jobdict[job]['exe'].keys() and \
            'nodes' in jobdict[job]['exe'].keys() and \
            'ppn' in jobdict[job]['exe'].keys()):
            # batch case - spec'd in envdict 
            if (jobdict[job]['exe']['platf'] in envdict.keys()):
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def isMPI(job,jobdict):
# test for presence of MPI info, for command line MPI
    MPIROOT = os.getenv('MPIROOT')
    if ('exe' in jobdict[job].keys()):
        # check spec of platform
        if ('platf' in jobdict[job]['exe'].keys() and \
            'ppn' in jobdict[job]['exe'].keys()):
            if (jobdict[job]['exe']['platf'] == 'mpi'):
                if MPIROOT == None:
                    print 'Note: MPIROOT not defined' 
                    print '  -> command line MPI not available'
                    return False
                else:
                    return True    
            else:
                return False
        else:
            return False
    else:
        return False
    
def getParamCommand(job,jobdict):
    # assigns full pathnames to dependencies, targets - overrides
    # target filenames in parfile by appending target pairs 
    if ('dep' in jobdict[job].keys()) and ('tgt' in jobdict[job].keys()):
        parlist = ['/bin/cat ' + os.path.join(THISPATH, job + '.par |')]
        for j in range(len(jobdict[job]['dep'])-1):
            parlist = parlist + ['sed s#' + jobdict[job]['dep'][j] + '#' + os.path.join(THISPATH,jobdict[job]['dep'][j]) + '# |']
        parlist = parlist + ['sed s#' + jobdict[job]['dep'][len(jobdict[job]['dep'])-1] + '#' + os.path.join(THISPATH,jobdict[job]['dep'][len(jobdict[job]['dep'])-1]) + '# > ./parfile']
        for k in jobdict[job]['tgt'].keys():
            parlist = parlist + \
                ['; echo ' + k + ' = ' + os.path.join(THISPATH,jobdict[job]['tgt'][k]) + ' >> ./parfile']
        parcmd =  ' '.join(parlist)
        return parcmd
    else:
        return

def getBatchAttributes(batch):
    # presumes that host shell is bash or other bourne-shell derivative
    # also that standard cluster directory structure is defined in
    # environment

    exp = []
    for i in ['USER', 'HOME', 'WORK', 'SCRATCH']:
        var = os.getenv(i)
        if var != None:
            exp = exp + ['export ' + i + '=' + var + '; ']
            
    syslist = ['slurm', 'pbs']
    bcmdlist = ['/usr/bin/sbatch', '/usr/bin/qsub']
    wcmdlist = ['--dependency=afterok:', '-W depend=afterok:']
    sufflist = ['bat', 'pbs']
    postlist = ["| grep 'batch job'|sed s/'Submitted batch job '//", '']

    found = False
    for i in range(len(syslist)):
        if batch == syslist[i]:
            bcmd = bcmdlist[i]
            wcmd = wcmdlist[i]
            suff = sufflist[i]
            post = postlist[i]
            found = True

    if not found:
        print 'Note: iwave.py supports only the following batch submission systems:'
        print ' '.join(syslist)
        return

    bshl = ' '.join(exp) + ' bash'
    
    bsys = {'bcmd' : bcmd, 
            'bshl' : bshl,
            'suff' : suff,
            'wcmd' : wcmd, 
            'post' : post
           }
    return bsys
        
