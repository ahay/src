#! /usr/bin/env python
'''Finds RSF files with missing or incomplete binaries or headers.
Delete them all with shell constructs like: rm -f `sfinvalid dir=.`'''

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

# This program is dependent on the output of sfin and sfattr

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
    verb = par.bool('verb', False)      # Display what is wrong with the dataset
    mydir = par.string('dir')           # Directory with files
    recursive = par.bool('rec', False)  # Whether to go down recursively
    chk4nan = par.bool('chk4nan',False) # Check for NaN values. Expensive!!

    if mydir == None:
        rsfprog.selfdoc()
        return ivlad.unix_error

    invalid_files_list = [] # will contain tuples: (file, msg)

    if recursive:
        for root, dirs, files in mydir:
            invalid_files_list += \
            ivlad.list_invalid_rsf_files(root, files, chk4nan)
    else:
        files = filter(lambda x:os.path.isfile(x),os.listdir(mydir))
        invalid_files_list += \
        ivlad.list_invalid_rsf_files(mydir, files, chk4nan)

    myline = ''

    for entry in invalid_files_list:
        myline = entry[0]
        if verb:
            myline += ': ' + entry[1]
        print myline

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main

