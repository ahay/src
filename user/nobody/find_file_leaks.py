#!/usr/bin/env python
'''
Finds C source code files that call sf_input for a given file, but forget to
call sf_close. Call this program from the upper directory (RSFSRC) with no
arguments, i.e.:
    ./admin/find_file_leaks.py > results.asc
It will print a list of affected files, sorted in a few categories.

To eliminate the file leak, make sure all file pointers are initialized with
NULL, i.e.:
    sf_file in=NULL;

Then call just before exit():
    sf_close();

This script checks that each file that has at least one sf_input also has a
sf_close call. It can be fooled by multi-procedure files, or by procedures
specially created to handle input. Flag such files by including them in the
'exceptions' variable in main().

Make sure that the list of directories to check for *.c files (variable
dirs_to_check at the beginning of this script) is not out of date.

This script has grep as a dependency.
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

try:
    import subprocess
except:
    print 'Please get Python 2.4 or greater, or just the subprocess module'
    sys.exit(unix_error)

###############################################################################

def grep_from_py(filename, string):
    'Runs "grep string filename", and returns the result'
    proc = subprocess.Popen('grep ' + string + ' ' + filename,
                            shell=True,
                            stdout=subprocess.PIPE)
    stdout_val = proc.communicate()[0]
    return stdout_val

###############################################################################

def print_list(file_list, what_is_list):
    'Formated display of list of suspect files'

    if file_list != []:
        msg = '\n-----------------------\nFound %d ' % len(file_list)
        msg += what_is_list + ':\n'
        print msg
        file_list.sort()
        for file in file_list:
            print file

###############################################################################

def main():

    dirs_to_check = 'api pens plot su system user'

    exceptions = '''
        api/matlab/rsf_write.c
        system/main/attr.c
        system/main/in.c
        system/main/mpi.c
        system/main/omp.c
        user/trip/wavefun.c
    '''

    exceptions_list = exceptions.split()

    primary_suspects    = [] # Main progs with sf_input, but no sf_fileclose
    secondary_suspects  = [] # other files with sf_input, but no sf_fileclose

    for d in dirs_to_check.split():
        for root,dirs,files in os.walk(d):
            # Avoid Subversion bookkeeping directories
            if '.svn' not in root:
                for file in files:
                    (shortname, extension) = os.path.splitext(file)
                    if extension == '.c':
                        c_file = os.path.join(root, file)
                        if c_file not in exceptions_list:
                            grep_4_sf_input = grep_from_py(c_file, 'sf_input')
                            if grep_4_sf_input != '': # file_c has sf_input
                                grep_4_sf_close = \
                                grep_from_py(c_file, 'sf_close')
                                if grep_4_sf_close == '': # no sf_fileclose
                                    if shortname[0] == 'M' or 'system/main/' in c_file:
                                        primary_suspects.append(c_file)
                                    else:
                                        secondary_suspects.append(c_file)

    print_list(primary_suspects,
    'main program C files containing sf_input, but not sf_fileclose')

    print_list(secondary_suspects,
    'other C files containing sf_input, but not sf_fileclose')

    return unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
