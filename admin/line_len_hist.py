#!/usr/bin/env python
'''
Outputs histogram with line length in Madagascar source code base for one type of file,
OR names of files which contains line longer than a given threshold

Line length does not include the newline character
Trailing whitespaces at the end of a line are not considered

Call this program from the upper directory (RSFSRC), i.e.:
    ./adm/line_len_hist.py ext=[py/c/cc/m/f90/f/txt] [pflt=<positive integer>] >results.asc

ext=py will also scan SConstruct files in book/ , (and nothing else in book/).
Other extensions will only go through source code directories.

pflt stands for "Print Files with Lines longer Than"
If a pflt value is specified, a list of filenames will be given
If no pflt is specified, output will have two columns: line length, nr of occurences
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

from rsf.user.ivlad import unix_error, unix_success

try:
    import rsf.api as rsf
except: # Madagascar's Python API not installed
    import rsf.apibak as rsf

import os, sys, glob

if not hasattr(os, 'walk'):
    print 'Please get Python 2.3 or greater'
    sys.exit(unix_error)

try:
    import subprocess
except:
    print 'Please get Python 2.4 or greater, or just the subprocess module'
    sys.exit(unix_error)

###########################################################################

def main():

    par  = rsf.Par(sys.argv)
    ext  = par.string('ext')  # Filetype by extension
    pflt = par.int('pflt', 0) # Print Files with Lines longer Than

    if ext == None:
        ext = 'c'

    dirs_to_search = 'adm api framework pens plot su system user'.split()

    if ext == 'py':
        dirs_to_search.append('book')

    if pflt > 0:
        guilty_files = []
    else:
        # This list will have at position i the number of lines of that length
        line_len_list = []

    for d in dirs_to_search:
        for root,dirs,files in os.walk(d):
            # Avoid Subversion bookkeeping directories
            if '.svn' not in root:
                myfiles = glob.glob(os.path.join(root,'*.'+ext))
                if ext == 'py':
                    local_sconstruct = os.path.join(root,'SConstruct')
                    if os.path.isfile(local_sconstruct):
                        myfiles.append(local_sconstruct)
                if myfiles != []:
                    for file in myfiles:
                        f = open(file,'r')
                        for line in f:
                            lilen = len(line.rstrip())
                            if pflt > 0:
                                if lilen > pflt:
                                    guilty_files.append(file)
                                    break
                            else:
                                # Make sure line_len_list is long enough
                                l4 = len(line_len_list)
                                ldiff = lilen+1-l4
                                if ldiff > 0:
                                    line_len_list.extend(ldiff*[0])
                                line_len_list[lilen] += 1
                        f.close()

    if pflt > 0:
        guilty_files.sort()
        for file in guilty_files:
            print file
    else:
        for i in range(len(line_len_list)):
            print '%d\t%d' % (i, line_len_list[i])

    return unix_success

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
