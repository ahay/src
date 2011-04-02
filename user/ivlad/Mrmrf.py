#! /usr/bin/env python
'''Recursively removes all RSF headers in a directory (associated binaries too)
Missing binaries do not cause failure.'''

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

import sys, os, glob

try: # Give precedence to local version
    import ivlad, m8rex, sf
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex
    import rsf.user.sf as sf

###############################################################################

def main(par):

    verb = par.bool('verb', False) # Display headers and binaries being deleted
    mydir = par.string('dir')           # Directory with files
    recursive = par.bool('rec', False)  # Whether to go down recursively

    # Clean up headers with existing binaries

    valid_files_list = []

    if recursive:
        for root, dirs, files in mydir:
            valid_files_list += \
            ivlad.list_valid_rsf_files(root, files, chk4nan)
    else:
        files = filter(lambda x:os.path.isfile(x),os.listdir(mydir))
        valid_files_list += \
        ivlad.list_valid_rsf_files(mydir, files, chk4nan=False)

    for f in valid_files_list:
        ivlad.msg(f + ': ' + ivlad.getout('sfin',['info=n',f]), verb)
        sf.rm(f, verb)

    # Clean up headers with no binaries

    if recursive:
        hdr_str = ivlad.getout('find',[mydir, '-type', 'f', '-name', '"*.rsf"'])
        hdr_list = hdr_str.split('\n')
    else:
        hdr_list = filter(lambda x:os.path.isfile(x),
                          glob.glob(os.path.join(mydir,'*.rsf')))

    for f in hdr_list:
        ivlad.msg(f)
        os.remove(f)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main, ['dir'])
