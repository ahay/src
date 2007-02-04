#! /usr/bin/env python
"""
NAME
	ivlad
DESCRIPTION
	Utilities for python metaprograms in RSFSRC/user/ivlad
SOURCE
	user/ivlad/ivlad.py
"""
# Copyright (C) 2007 Ioan Vlad
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

import sys
error = 1
success = 0

# Try to avoid old platform-dependent modules that will be deprecated.
# See http://docs.python.org/dev/lib/module-subprocess.html
try:
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False
    from os import system
    from commands import getoutput


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
    else: # no want= given
        if not stdout:
            sys.stderr.write('send_to_os must receive either stdout= or want=')
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
            subprocess.Popen(cmdlist,stdin=finp,stdout=fout).wait()
            return None
        elif want == 'stdout':
            s = subprocess.Popen(cmdlist,stdin=finp,stdout=subprocess.PIPE)
            output = s.communicate()[0]
            return output.rstrip('\n') # Remove the newline character
    else: # no subprocess module present
        if want == 'stdout':
            s = getoutput( command )
            return s
        elif stdout:
            system(command)
            return None