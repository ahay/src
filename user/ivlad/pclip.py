#! /usr/bin/env python
"""NAME
	sfpclip
DESCRIPTION
	Percentile clip. Shell for sfclip(sfquantile(input)).
SYNOPSIS
	sfpclip inp= out= verb=n pclip=99.0
PARAMETERS
	string inp=		Input file
	string out=		Output file
	bool   verb=n [y/n] 	If y, print system commands, outputs
	float  pclip=99.0	Percentile clip
SOURCE
	user/ivlad/pclip.py
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

try:
    import rsf
except:
    import rsfbak as rsf

try:
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False
    from os import system
    from commands import getoutput

def assemble_cmd(cmd, arg, stdin, stdout):
    command = cmd+' '+arg
    if stdin:
        command = '< ' + stdin + ' ' + command
    if stdout:
        command = command+'>'+stdout
    return command

def print_cmd(verb, cmd, arg, stdin, stdout):
    if verb:
        command=assemble_cmd(cmd, arg, stdin, stdout)
        print(command)

def main(argv=sys.argv):

    success = 0
    error   = 1
    par = rsf.Par(argv)

    inp = par.string('inp')
    out = par.string('out')
    if inp == None or out == None:
        sys.stderr.write(__doc__)
        return error

    verb  = par.bool('verb', False)
    pclip = par.float('pclip',99)

    if pclip <0 or pclip>100:
        sys.stderr.write('pclip must be between 0 and 100\n')
        return error

    cmd = 'sfquantile'
    arg = 'pclip=' + str(pclip)
    finp= open(inp,'r')

    if have_subprocess:
        if verb:
            print_cmd(verb, cmd, arg, inp, None)
        clip=subprocess.Popen([cmd,arg],stdin=finp,stdout=subprocess.PIPE).communicate()[0]
    else:
        command = assemble_cmd(cmd, arg, inp, None)
        if verb:
            print command
        clip = getoutput( command )

    if not clip:
        sys.stderr.write('sfquantile did not return anything!\n')
        return error

    if verb:
        print(clip)

    cmd = 'sfclip'
    arg = 'clip=' + str(clip)
    finp= open(inp,'r')
    fout= open(out,'w')

    if have_subprocess:
        if verb:
            print_cmd(verb, cmd, arg, inp, out)
        clip=subprocess.Popen([cmd,arg],stdin=finp,stdout=fout).wait()
    else:
        command = assemble_cmd(cmd, arg, inp, out)
        if verb:
            print command
        system(command)

    return success

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main