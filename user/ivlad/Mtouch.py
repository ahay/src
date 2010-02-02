#! /usr/bin/env python
'''Applies the Unix command touch to binaries of RSF datasets in a directory.
Will go down recursively in subdirectories. Current date and time is used.
Useful for determining disk leaks: Orphan binaries (those without headers) will
not be touched. You can remove them with commands such as:
find $DATAPATH -type f -mmin -15 -exec rm -f {} \;'''

# Copyright (C) 2010 Ioan Vlad
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

# This program is dependent on the output of sfin info=n

import sys, os, glob, commands
import rsfprog

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

def visit(verb, dirname, names):
    'Function to execute once in each dir'
    for f in names:
        if os.path.splitext(f)[1] == ivlad.ext: # File has rsf extension 
            fname = os.path.abspath(os.path.join(dirname,f))
            cmd = 'sfin info=n ' + fname
            bfile = commands.getoutput(cmd)
            if(verb):
                print fname + ': \t' + bfile
            cmd = 'touch -c ' + bfile # Do not create if it does not exist
            os.system(cmd)

###############################################################################

def main(argv=sys.argv):

    par = rsf.Par(argv)
    verb = par.bool('verb', False) # display header and corresponding binary
    mydir = par.string('dir') # directory with files

    if mydir == None:
        rsfprog.selfdoc()
        return ivlad.unix_error
    if not os.path.isdir(mydir):
        print mydir + ' is not a valid directory'
        return ivlad.unix_error
    if not os.access(mydir,os.X_OK):
        print mydir + ' lacks +x permissions for ' + os.getlogin()
        return ivlad.unix_error
    if not os.access(mydir,os.R_OK):
        print mydir + ' lacks read permissions for ' + os.getlogin()
        return ivlad.unix_error

    os.path.walk(mydir, visit, verb)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main

