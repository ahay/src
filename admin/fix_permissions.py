#!/usr/bin/env python
'''
Sets permissions properly for files that is known for sure to not be executable

Call this program from the upper directory (RSFSRC) with no arguments, i.e.:
    ./admin/fix_permissions.py
'''
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

from rsf.user.ivlad import unix_success, unix_error
import os, sys

if not hasattr(os, 'walk'):
    print 'Please get Python 2.3 or greater'
    sys.exit(unix_error)

###############################################################################

def main():

    dirs_to_check = 'api book framework pens plot su system user .'
    nonexec_ext = 'c cc f f90 h i m tex txt'
    nonexec_files = 'SConstruct Makefile'

    nonexec_ext_list   = nonexec_ext.split()
    nonexec_files_list = nonexec_files.split()

    for d in dirs_to_check.split():
        for root,dirs,files in os.walk(d):
            # Avoid Subversion bookkeeping directories
            if '.svn' not in root:
                for file in files:
                    full_path = os.path.join(root, file)
                    if not os.path.islink(full_path):
                        nonexec = False
                        if file in nonexec_files_list:
                            nonexec = True
                        (shortname, extension) = os.path.splitext(file)
                        if extension.lstrip('.') in nonexec_ext_list:
                            nonexec = True
                        if nonexec:
                            os.system('chmod a-x ' + full_path)

    return unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
