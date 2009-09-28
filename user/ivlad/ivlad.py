#! /usr/bin/env python
"""
NAME
	ivlad
DESCRIPTION
	Utilities for python metaprograms in RSFSRC/user/ivlad
SOURCE
	user/ivlad/ivlad.py
"""
# Copyright (C) 2007, 2009 Ioan Vlad
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

import os, sys, commands, math

# Operating system return codes
unix_success = 0
unix_error   = 1

# Other constants
pipe = ' | '
ext  = '.rsf'

# Try to avoid old platform-dependent modules that will be deprecated.
# See http://docs.python.org/dev/lib/module-subprocess.html
try:
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False
    from os import system
    from commands import getoutput

###############################################################################

def send_to_os(prog, arg=None, stdin=None, stdout=None, want=None, verb=False):
    '''Sends command to the operating system. Arguments:
    - prog. Executable to be run. STRING. The only non-optional argument.
    - arg. List of strings with program arguments. LIST (can be STRING for only 1 arg)
    - stdin. Filename STRING.
    - stdout. Filename STRING. Must not be specified if want=stdout
    - want: What to return. STRING: 'stdout' or 'stderr'. [Stderr not implemented yet]
    - verb: whether to print command before executing it. BOOL
    If stdout= is given, then the function writes to that file and returns None
    If want=stdout, then a STRING with stripped newlines is returned.'''

    if want:
        if want not in ('stdout','stderr'):
            sys.stderr.write('The "want" argument to send_to_os must be "stdout" or "stderr"')
            sys.exit(error)
        if stdout and want == 'stdout':
            sys.stderr.write('No stdout= should be given to send_to_os when want="stdout"')
            sys.exit(error)

    # Build the [prog, args] list
    if arg:
        if str(arg) == arg: # it's a string
            arg = [arg]     # make it a list
        arg.insert(0,prog)
        cmdlist = arg
    else:
        cmdlist = [prog]

    # Build command string for printing or Python < 2.4
    if verb or not have_subprocess:
        cmd4print = cmdlist[:]
        if stdin:
            cmd4print.append('<')
            cmd4print.append(stdin)
        if stdout:
            cmd4print.append('>')
            cmd4print.append(stdout)
        command = ' '.join(cmd4print)
        if verb:
            print command

    if have_subprocess:
        if stdin:
            finp = open(stdin,'r')
        else:
            finp = None
        if stdout:
            fout = open(stdout,'w')
        else:
            fout = None
        if want == 'stdout':
            s = subprocess.Popen(cmdlist,stdin=finp,stdout=subprocess.PIPE)
            output = s.communicate()[0]
            return output.rstrip('\n') # Remove the newline character
        else:
            subprocess.Popen(cmdlist,stdin=finp,stdout=fout).wait()
            return None
    else: # no subprocess module present
        if want == 'stdout':
            s = getoutput( command )
            return s
        elif stdout:
            system(command)
            return None

###############################################################################

def readaxis( inp, axisnr, verb=False ):
    '''Reads n,o,d for one axis of a Madagascar hypercube
    - inp: filename. STRING.
    - axisnr: INT
    - verb: BOOL'''

    ax = str(axisnr)

    n = send_to_os('sfget',
                   arg = ['parform=n','n'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    o = send_to_os('sfget',
                   arg = ['parform=n','o'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    d = send_to_os('sfget',
                   arg = ['parform=n','d'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    return( int(n), float(o), float(d))

###############################################################################

def ndims(filename):
    'Returns nr of dims in a file, ignoring trailing length 1 dimensions'

    max_dims = 9
    nlist = [1] * max_dims

    for dim in range(1,max_dims+1):
        command = '<'+filename+' sfget parform=n n' + str(dim)
        sfget_out = commands.getoutput(command)
        if sfget_out[:15] != 'sfget: No key n':
            curr_n = int(sfget_out)
            if curr_n > 1:
                nlist[dim-1] = curr_n

    for dim in range(max_dims,0,-1):
        if nlist[dim-1] > 1:
            break

    return dim

####################################################################

def add_zeros(i, n):
    '''Generates string of zeros + str(i)'''

    ndigits_n = int(math.floor(math.log10(n-1)))
    if i == 0:
        nzeros = ndigits_n
    else:
        nzeros = ndigits_n - int(math.floor(math.log10(i)))

    return nzeros*'0'+str(i)

####################################################################

def execute(command, verb=False):
    '''Echoes a command to screen, then executes it'''

    if verb:
        print command
    os.system(command)
