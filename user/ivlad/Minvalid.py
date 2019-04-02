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

from __future__ import print_function
import os, sys

try: # Give precedence to local version
    import ivlad, m8rex
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex

###############################################################################

def main(par):

    verb = par.bool('verb', False)      # Display what is wrong with the dataset
    mydir = par.string('dir', '.')      # Directory with files
    recursive = par.bool('rec', False)  # Whether to go down recursively
    chk4nan = par.bool('chk4nan',False) # Check for NaN values. Expensive!!

    invalid_files_list = [] # will contain tuples: (file, msg)

    if recursive:
        for root, dirs, files in os.walk(mydir):
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
        print(myline)# This program's stdout is a list of files. Do not use msg

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main)
