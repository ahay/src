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

def getFlowSignature(workdir, jobdict, envdict):

    # batch execution
    if isBatch(jobdict,envdict):
        # detect batch system attributes
        bsys = getBatchAttributes(envdict[jobdict['exe']['platf']]['batch'])
        if bsys == None:
            print 'iwave.signature: failed to parse batch submission system'
            return

        # create batch file names
        thisbat = jobdict['job'] + '.' + bsys['suff']
        
        # write batch files
        writeScript(envdict[jobdict['exe']['platf']]['batch'],
                    jobdict['job'],
                    workdir,	
                    envdict[jobdict['exe']['platf']]['queue'],
                    envdict[jobdict['exe']['platf']]['mail'],
                    envdict[jobdict['exe']['platf']]['acct'],
                    jobdict['exe']['wall'],
                    jobdict['exe']['nodes'],
                    jobdict['exe']['ppn'],
                    jobdict['cmd'])

        workcmd = jobdict['pre'] + \
            '; /bin/mv ' + os.path.join(THISPATH,thisbat) + \
            ' ./' + thisbat  + \
            '; /bin/rm -f jobid; ' + \
            bsys['bshl'] + envdict[jobdict['exe']['platf']]['local'] + \
            ' bash -c "' + bsys['bcmd'] + ' ' + thisbat + '"' + \
            '; (while [ ! -f ./jobid ]; do (sleep 1); done)'

#            bsys['bcmd'] + ' ' + thisbat + '"' \

    elif isMPI(jobdict):

        MPIROOT = os.getenv('MPIROOT')
        mpirun = os.path.join(MPIROOT,'bin/mpirun')
        # this should be checked for existence
        workcmd = jobdict['pre'] + '; cd ' + workdir + \
            '; ' + mpirun + ' -np ' + jobdict['exe']['ppn'] + \
            ' ' + jobdict['cmd']
        
    else:      

        workcmd = jobdict['pre'] + '; cd ' + workdir + \
            '; ' + jobdict['cmd']

    workdict = { 'tgt' : jobdict['tgt'],
                 'src' : jobdict['src'],
                 'cmd' : workcmd,
                 'dir' : workdir,
                 'wcd' : '/bin/rm -rf ' + workdir + '; /bin/mkdir ' + workdir
               }
    return workdict

def writeScript(batch,job,path,queue,mail,acct,wall,nodes,ppn,exe):
    if batch == 'pbs':
        f = open(job + '.pbs','w')
        f.write('#PBS -N ' + job + '\n')
	f.write('#PBS -A ' + acct + '\n')
        f.write('#PBS -V\n')
        f.write('#PBS -l walltime=' + wall + '\n')
        f.write('#PBS -l nodes=' + nodes + ':ppn=' + ppn + '\n')
        f.write('#PBS -q '+queue + '\n')
        f.write('#PBS -m abe\n')
        f.write('#PBS -M '+ mail + '\n')
        f.write('#PBS -o '+os.path.join(THISPATH,path)+'/cout.txt\n')
        f.write('#PBS -e '+os.path.join(THISPATH,path)+'/cerr.txt\n')
        f.write('cd '+ os.path.join(THISPATH,path)+'\n')
        f.write('mpiexec ' + exe + '\n')
        f.write('echo $PBS_JOBID > jobid\n')
        f.close()

    if batch == 'slurm':
        f = open(job + '.bat','w')
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
        f.write('#SBATCH -o '+os.path.join(THISPATH,path)+'/cout.txt\n')
        f.write('#SBATCH -e '+os.path.join(THISPATH,path)+'/cerr.txt\n')
        f.write('cd '+os.path.join(THISPATH,path)+'\n')
        f.write('ibrun ' + exe + '\n')
        f.write('echo $SLURM_JOBID > jobid\n')
        f.close()

def isBatch(jobdict,envdict):
    if ('exe' in jobdict.keys()):
        # check spec of platform
        if ('platf' in jobdict['exe'].keys() and \
            'wall' in jobdict['exe'].keys() and \
            'nodes' in jobdict['exe'].keys() and \
            'ppn' in jobdict['exe'].keys()):
            # batch case - spec'd in envdict 
            if (jobdict['exe']['platf'] in envdict.keys()):
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def isMPI(jobdict):
# test for presence of MPI info, for command line MPI
    MPIROOT = os.getenv('MPIROOT')
    if ('exe' in jobdict.keys()):
        # check spec of platform
        if ('platf' in jobdict['exe'].keys() and \
            'ppn' in jobdict['exe'].keys()):
            if (jobdict['exe']['platf'] == 'mpi'):
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
    
def getBatchAttributes(batch):
    # presumes that host shell is bash or other bourne-shell derivative
    # also that standard cluster directory structure is defined in
    # environment

    exp = []
    #    for i in ['USER', 'HOME', 'WORK', 'SCRATCH']:
    #        var = os.getenv(i)
    #        if var != None:
    #            exp = exp + ['export ' + i + '=' + var + '; ']

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

    bshl = ' '.join(exp)
    
    bsys = {'bcmd' : bcmd, 
            'bshl' : bshl,
            'suff' : suff,
            'wcmd' : wcmd, 
            'post' : post
           }
    return bsys

def getThreads(dict):
    nthread=1
    if ('ppn' in dict.keys()):
        nthread = nthread * int(dict['ppn'])
    if ('nodes' in dict.keys()):
        nthread = nthread * int(dict['nodes'])
    return str(nthread)
                                                       
