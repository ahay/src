#! /usr/bin/env python
'Finds RSF files with missing or incomplete binaries or headers'

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

def main(argv=sys.argv):

    par = rsf.Par(argv)
    verb = par.bool('verb', False) # verbosity flag
    mydir = par.string('dir') # directory with files
    if mydir == None:
        rsfprog.selfdoc()
        return ivlad.unix_error
    else:
        if not os.path.isdir(mydir):
            print mydir + ' is not a valid directory'
            return ivlad.unix_error
        if not os.access(mydir,os.X_OK):
            print mydir + ' lacks +x permissions for ' + os.getlogin()
            return ivlad.unix_error
    if not os.access(mydir,os.R_OK):
        print mydir + ' lacks read permissions for ' + os.getlogin()
        return ivlad.unix_error

    if os.path.abspath(mydir) == os.getcwd():
        mydir = '' # for a more readable display

    rsf_list = glob.glob(os.path.join(mydir,'*.rsf'))

    # Filter out metafile directories (yes, some people use them)
    rsf_files = filter(lambda x:os.path.isfile(x),rsf_list)
    rsf_files.sort()

    sfin = 'sfin info='

    for f in rsf_files:
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
            if line_list[-1][:6] == 'sfin: ':
                msg = f
                if(verb):
                    msg += '\tBinary '
                    msg += line_list[-1][6:].strip().split(',')[1].rstrip('.').lstrip()
        if msg != None:
            print msg

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
