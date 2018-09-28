#!/usr/bin/env python

# Copyright (C) 2007 University of Texas at Austin
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

import sys, string, os, signal, types

def handler(signum, frame):
    'signal handler for abortion [Ctrl-C]'
    sys.stderr.write('\n[Ctrl-C] Aborting...\n')
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
        if type(comm) is list:
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
            mycomm = mycomm.replace('%',subdir,1)
            syswait(mycomm)
        sys.stderr.write(string.join(['-' * 44,subdir,'\n'],' '))
    sys.stderr.write('Done.\n')
    os.chdir (cwd)

if __name__ == "__main__":
    import glob
    
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    if len(sys.argv) < 2:
        print('''
        Usage: %s [-q] command
        visits lower-case subdirectories and executes command
        -q     quiet (suppress stdout)

        The '%%' character is replaced with the current directory
        ''' % sys.argv[0])
        sys.exit(0)

    ########
    comm = sys.argv.pop(0)

    if comm == "-q":
        verbose = 0
        comm = sys.argv.pop(0)
    else:
        verbose = 1

    comm = string.join(sys.argv,' ')
    dirs = [x for x in list(filter(os.path.isdir,glob.glob('[a-z]*'))) if x[-5:] != '_html']

    tour(sorted(dirs),comm,verbose)
    sys.exit(0)
