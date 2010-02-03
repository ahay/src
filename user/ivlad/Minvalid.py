#! /usr/bin/env python
'''Finds RSF files with missing or incomplete binaries or headers.
Delete them all with shell constructs like: rm -f `sfinvalid dir=.`'''

# Copyright (C) 2009 Ioan Vlad
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

# This program is dependent on the output of sfin

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

###############################################################################

def visit(verb, dirname, names):
    'Function to execute once in each dir'

    sfin = 'sfin info='
    rsf_files = filter(lambda x:os.path.splitext(x)[1]==ivlad.ext, names)
    rsf_files.sort()

    for f in rsf_files:
        f = os.path.abspath(os.path.join(dirname,f))
        com_out = commands.getoutput(sfin+'n '+f)
        msg = None
        if com_out[:6] == 'sfin: ':
            msg = f
            if(verb):
                msg += '\t'
                err_msg = com_out.split(':')[2].strip()
                if err_msg[:15] == 'No in= in file ':
                    # Check if header is zero
                    header_sz = os.path.getsize(f)
                    if header_sz == 0:
                        msg += 'Empty header'
                    else:
                        msg += 'Missing in='
                elif err_msg[:22] == 'Cannot read data file ':
                    # Check the binary
                    bnr = err_msg[22:]
                    (bnr_path, bnr_file) = os.path.split(bnr)
                    if not os.access(bnr_path,os.X_OK):
                        msg += 'Missing +x permissions to dir: ' + bnr_path
                    elif not os.path.isfile(bnr):
                        msg += 'Missing binary: ' + bnr
                    elif not os.access(bnr,os.R_OK):
                        msg += 'Missing read permissions to file: ' + bnr
                else: # Error not classified in this script
                    msg += err_msg
        else:
            # Check for incomplete binaries
            com_out = commands.getoutput(sfin+'y '+f)
            line_list = com_out.split('\n')

            # Using a grep equivalent below because the "Actually..." error
            # message in sfin is not always on last line.

            # Clumsy interpretation of grep below, with a touch of awk
            # Using Unix grep in original command will not work
            # because errors come via stderr and grep works on stdout

            err_list = []
            for line in line_list:
                if line[:6] == 'sfin: ':
                    err_msg = line[6:]
                    if err_msg != 'This data file is entirely zeros.':
                        if err_msg[:10] != 'The first ':
                            err_list.append(err_msg.strip())
            if err_list != []:
                msg = f
                if verb:
                    # Showing only the first error. Probably the only one.
                    msg += '\t' + err_list[0]

        if msg != None:
            print msg

###############################################################################

def main(argv=sys.argv):

    par = rsf.Par(argv)
    verb = par.bool('verb', False)     # Display what is wrong with the dataset
    mydir = par.string('dir')          # directory with files
    recursive = par.bool('rec', False) # Whether to go down recursively

    # I should really use exceptions here. No time now.
    if ivlad.chk_dir(mydir) == ivlad.unix_error:
        return ivlad.unix_error

    if recursive:
        os.path.walk(mydir, visit, verb)
    else:
        files = filter(lambda x:os.path.isfile(x),os.listdir(mydir))
        visit(verb, mydir, files)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main

