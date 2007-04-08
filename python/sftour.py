#!/usr/bin/env python
import sys, string, os, signal, types

def handler(signum, frame):
    'signal handler for abortion [Ctrl-C]'
    sys.stderr.write('\n%s: aborting...\n' % comm)
    if child:
        os.kill (signal.SIGINT,child)
    sys.exit(-1)

signal.signal(signal.SIGINT,handler) # handle interrupt

child = None
def syswait(comm):
    'Interruptable system command'
    global child
    child = os.fork()
    if child:
        (pid,exit) = os.waitpid(child,0)
        child = 0
        return exit
    else:
        os.system(comm)
        os._exit(0)

def tour(dirs=[],comm='',verbose=1):
    'Visit every directory in dirs running a command comm'
    if not verbose: # no output to stdout
        sys.stdout = open("/dev/null","w")
    sys.stderr.write('Executing "%s"...\n' % comm)
    sys.stderr.write(string.join(dirs,'::') + '\n')

    cwd = os.getcwd()
    for subdir in dirs:
        if type(comm) is types.ListType:
            mycomm = comm.pop(0)
        else:
            mycomm = comm
            
        os.chdir (cwd)
        try:
            os.chdir (subdir)
        except:
            sys.stderr.write('\n%s: wrong directory %s...\n' % (mycomm,subdir))
            sys.exit(1)
        os.environ['PWD'] = os.path.join(cwd,subdir)
        sys.stderr.write(string.join(['+' * 44,subdir,'\n'],' '))
        if mycomm:
            syswait(mycomm)
        sys.stderr.write(string.join(['-' * 44,subdir,'\n'],' '))
    sys.stderr.write('Done.\n')
    os.chdir (cwd)

if __name__ == "__main__":
    import glob
    
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    if len(sys.argv) < 2:
        print '''
        Usage: %s [-q] command
        visits lower-case subdirectories and executes command
        -q     quiet (suppress stdout)
        ''' % sys.argv[0]
        sys.exit(0)

    ########
    comm = sys.argv.pop(0)

    if comm == "-q":
        verbose = 0
        comm = sys.argv.pop(0)
    else:
        verbose = 1

    comm = string.join(sys.argv,' ')
    dirs = filter(lambda x: x[-5:] != '_html',
                  filter(os.path.isdir,glob.glob('[a-z]*')))

    tour(dirs,comm,verbose)
    sys.exit(0)
