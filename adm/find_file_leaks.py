#!/usr/bin/env python
'''
Finds C source code files that call sf_input for a given file, but forget to
call sf_fileclose. Call this program from the upper directory (RSFSRC) with no arguments, i.e.:
    ./adm/find_file_leaks.py > results.asc
It will print a list of affected files, sorted in a few categories.

To eliminate the file leak, make sure all file pointers are initialized with
NULL, i.e.:
    sf_file foo=NULL;
Then close the file only if the pointer has been allocated:
    if (foo != NULL) sf_fileclose(foo);
Keep in mind that this script checks for occurences of sf_input and
sf_fileclose in a file, not in a procedure. Multi-procedure files, or
procedures created especially to handle input, can easily "fool" this filter.

Make your life easier by placing already-checked files in the similarly-named
variable in main().

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

unix_success = 0
unix_error   = 1

import os, sys

if not hasattr(os, 'walk'):
    print 'Please get Python 2.3 or greater'
    sys.exit(unix_error)

try:
    import subprocess
except:
    print 'Please get Python 2.4 or greater, or just the subprocess module'
    sys.exit(unix_error)

def grep_from_py(filename, string):
    'Runs "grep string filename", and returns the result'
    proc = subprocess.Popen('grep ' + string + ' ' + filename,
                            shell=True,
                            stdout=subprocess.PIPE)
    stdout_val = proc.communicate()[0]
    return stdout_val

def print_list(file_list, what_is_list):
    'Formated display of list of suspect files'

    if file_list != []:
        msg = '\n-----------------------\nFound %d ' % len(file_list)
        msg += what_is_list + ':\n'
        print msg
        file_list.sort()
        for file in file_list:
            print file

def main():

    dirs_to_check = 'api pens plot su system user'

    already_checked_files = '''
    plot/main/graph.c
    plot/main/grey.c
    plot/main/grey3.c
    system/generic/Mfft1.c
    system/generic/Mfft3.c
    system/generic/Mnoise.c
    system/generic/Mremap1.c
    system/generic/Msmooth.c
    system/generic/Munif3.c
    system/main/add.c
    system/main/attr.c
    system/main/dd.c
    system/main/disfil.c
    system/main/mask.c
    system/main/pad.c
    system/main/real.c
    system/main/reverse.c
    system/main/rotate.c
    system/main/rtoc.c
    system/main/spray.c
    system/main/stack.c
    system/main/transp.c
    system/main/window.c
    system/seismic/Mdepth2time.c
    system/seismic/Mai2refl.c
    system/seismic/Mricker.c
    system/seismic/Mricker1.c
    system/seismic/Mtime2depth.c
    system/seismic/Mzomig.c
    user/jennings/Mclip2.c
    user/jennings/Mlistminmax.c
    user/jennings/Mminmax.c
    user/psava/Msrmig3.c
    user/psava/Msrmod3.c
    '''

    already_checked_list = already_checked_files.split()

    primary_suspects    = [] # Main progs with sf_input, but no sf_fileclose
    secondary_suspects  = [] # other files with sf_input, but no sf_fileclose
    tertiary_suspects   = [] # files w. mismatch between input and close calls
    quaternary_suspects = [] # Files with both sf_input and sf_fileclose

    for d in dirs_to_check.split():
        for root,dirs,files in os.walk(d):
            # Avoid Subversion bookkeeping directories
            if '.svn' not in root:
                for file in files:
                    (shortname, extension) = os.path.splitext(file)
                    if extension == '.c' and file not in already_checked_list:
                        c_file = os.path.join(root, file)
                        grep_4_sf_input = grep_from_py(c_file, 'sf_input')
                        if grep_4_sf_input != '': # file_c has sf_input
                            grep_4_sf_fileclose = \
                            grep_from_py(c_file, 'sf_fileclose')
                            if grep_4_sf_fileclose == '': # no sf_fileclose
                                if shortname[0] == 'M' or 'system/main/' in c_file:
                                    primary_suspects.append(c_file)
                                else:
                                    secondary_suspects.append(c_file)
                            else: # file has both sf_input and sf_fileclose
                                # Assumption: only one sf_input per line and
                                #             only one sf_fileclose per line
                                if grep_4_sf_input.count    ('\n') == \
                                   grep_4_sf_fileclose.count('\n'):
                                    quaternary_suspects.append(c_file)
                                else:
                                    tertiary_suspects.append(c_file)

    print_list(primary_suspects,
    'main program C files containing sf_input, but not sf_fileclose')

    print_list(secondary_suspects,
    'other C files containing sf_input, but not sf_fileclose')

    print_list(tertiary_suspects,
    'C files with a mismatch between the number of sf_input and sf_fileclose calls')

    print_list(quaternary_suspects,
    'C files with as many sf_input as sf_fileclose calls')

    print_list(already_checked_list,
    'already checked files')

    return unix_success

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
