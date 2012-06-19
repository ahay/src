#!/usr/bin/env python
'''Inventory test results of Madagascar example book directories.
Scan a directory tree (or list of directory trees) to a given depth,
and inventory the contents and test results.

The inventory occurs only at the specified depth (default 3), not the
intervening depths.  Directories named .svn are skipped.  Only
directories containing an SConstruct file are listed.

Examples (from within $RSFSRC):

sftestlist book                         # inventory of book
sftestlist levels=2 book/geostats       # inventory of book/geostats
'''

# Copyright (C) 2009 James W. Jennings Jr.
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

import os, sys
import rsf.prog as rsfprog

if not hasattr(os,'walk'):
    sys.stderr.write('testlist needs Python 2.3 or greater\n')
    sys.exit(unix_error)

import rsf.api as rsf

sfprefix = 'sf'                 # prefix for rsf commands
plprefix = 'vp'                 # prefix for vpl commands
sfsuffix = '.rsf'               # suffix for rsf files
vpsuffix = '.vpl'               # suffix for vplot files

###############################################################################

def size_string(size):
    'Make a human-readable size string'

    units = [' B','KB','MB','GB','TB']
    p = 1
    while 1024**p <= size : p = p+1
    return "%4.2f %s" % (float(size)/1024**(p-1),units[p-1])

def int_string(value):
    'Format an integer into a string'

    if value == None:   return '-'
    else:               return '%d' % value

def read_rsfproj(root,files):
    'Read contents of an .rsfproj file'

                                    # initialize
    values = {  'error':    0, 
                'exist':    False,
                'uses':     [],
                'data':     'unknown',
                'size':     None}

    if '.rsfproj' in files:         # get info from .rsfproj
        values['exist'] = True

        g = {}                      # read file
        l = {}
        execfile(os.path.join(root,'.rsfproj'),g,l)

                                    # get uses
        if 'uses' in l: values['uses']  = l['uses']
        else:           values['error'] = 1

                                    # get data type
        if 'data' in l:
            if len(l['data']) == 0:     values['data'] = 'none   '
            if len(l['data']) >  0:     values['data'] = 'public '
            if 'PRIVATE' in l['data']:  values['data'] = 'private'
            if 'LOCAL'   in l['data']:  values['data'] = 'local'
        else:
            values['error'] = 1

                                    # get data size
        if 'size' in l: values['size']  = l['size']
        else:           values['error'] = 1

    return values

def read_rsftest(root,files):
    'Read contents of an .rsftest file'

                                    # initialize
    items = ['time','exist_fig','exist_lock','miss','extra','diff','same']

    values = {'error':0, 'exist':False}
    
    for item in items:  values[item] = None
               
    if '.rsftest' in files:         # get info from .rsftest
        values['exist'] = True

        g = {}                      # read file
        l = {}
        execfile(os.path.join(root,'.rsftest'),g,l)
                
                                    # load info into values dictionary
        for item in items:
            if item in l:   values[item]    = l[item]
            else:           values['error'] = 1

    return values

###############################################################################

def main(argv=sys.argv):

    unix_success = 0
    unix_error   = 1

################    get user parameters

    if len(argv) < 2 :              # print selfdoc and exit if no parameters
        rsfprog.selfdoc()
        return unix_error

    par = rsf.Par(argv)             # get parameters

    levels        = par.int('levels',3)
    # directory search depth
    
    outfile_name  = par.string('outfile')
    # file name for detailed inventory table, default none

    untested      = par.bool('untested',False)
    # list untested examples?

################    build list of search directories

    books = []
    for item in argv[1:]:
        if '='  in item: continue       # skip parameter items
        if not os.path.exists(item):    # does item exist?
            sys.stderr.write("Directory '%s' does not exist.\n" % item)
            continue
        if not os.path.isdir(item):     # is item a directory?
            sys.stderr.write("File '%s' is not a directory.\n" % item)
            continue
        books.append(os.path.normpath(item))

    if len(books) == 0:                 # nothing to search
        sys.stderr.write("No directories to search.\n")
        return unix_error

    sys.stdout.write("Searching in: %s\n\n" % books)

################    open outut file

    if outfile_name == None:
        outfile = None
    else:
        outfile = open(outfile_name,'w')
        outfile.write('proj error  ')
        outfile.write('test error  time                                ')
        outfile.write('size  type     fig lock miss  ext diff same  dir\n')

################    search directory tree

    sbool  = {True:'yes', False:'no ', None:'-  '}
    serror = {0:'ok   ', 1:'error'}

    rsftest_errors  = 0
    rsftest_no      = 0
    figdir_no       = 0
    lockdir_no      = 0
    figs_missing    = 0
    figs_extra      = 0
    figs_diff       = 0

    rsftest_error_list  = []
    rsftest_no_list     = []
    figdir_no_list      = []
    lockdir_no_list     = []
    figs_missing_list   = []
    figs_extra_list     = []
    figs_diff_list      = []

    for bookdir in books:
        for root,dirs,files in os.walk(bookdir):

            if '.svn' in dirs: dirs.remove('.svn')      # don't visit .svn dirs

            reldir = root.replace(bookdir,"dummy",1)    # get relative dir name
            level  = len(reldir.split(os.sep))-1        # get dir level

            if level <  levels: continue        # scan to depth 'levels' ...
            if level == levels: del dirs[:]     # ... but not beyond 'levels'

                                                # skip dirs without SConstruct
            if 'SConstruct' not in files: continue

            rsfproj_vals = read_rsfproj(root,files)     # read rsfproj file
            rsftest_vals = read_rsftest(root,files)     # read rsftest file

                                                        # count errors
            if rsftest_vals['error'] == 1:
                rsftest_errors  = rsftest_errors+1
                rsftest_error_list.append(root)

            if not rsftest_vals['exist']:
                rsftest_no  = rsftest_no+1
                rsftest_no_list.append(root)

            if rsftest_vals['exist_fig'] == False:
                figdir_no  = figdir_no+1
                figdir_no_list.append(root)

            if rsftest_vals['exist_lock'] == False:
                lockdir_no  = lockdir_no+1
                lockdir_no_list.append(root)

            if rsftest_vals['miss'] > 0:
                figs_missing = figs_missing+1
                figs_missing_list.append(root)

            if rsftest_vals['extra'] > 0:
                figs_extra = figs_extra+1
                figs_extra_list.append(root)

            if rsftest_vals['diff'] > 0:
                figs_diff = figs_diff+1
                figs_diff_list.append(root)

                                                # write summary table
            if outfile != None:
                outfile.write('%s  '      % sbool[rsfproj_vals['exist']])
                outfile.write('%s  '      % serror[rsfproj_vals['error']])
                outfile.write('%s  '      % sbool[rsftest_vals['exist']])
                outfile.write('%s  '      % serror[rsftest_vals['error']])
                outfile.write('"%-24s"  ' % rsftest_vals['time'])
                outfile.write('%12s  '    % int_string(rsfproj_vals['size']))
                outfile.write('%s  '      % rsfproj_vals['data'])
                outfile.write('%s  '      % sbool[rsftest_vals['exist_fig']])
                outfile.write('%s  '      % sbool[rsftest_vals['exist_lock']])
                outfile.write('%3s  '     % int_string(rsftest_vals['miss']))
                outfile.write('%3s  '     % int_string(rsftest_vals['extra']))
                outfile.write('%3s  '     % int_string(rsftest_vals['diff']))
                outfile.write('%3s  '     % int_string(rsftest_vals['same']))
                outfile.write('%s\n'      % root)
            
################    write summary

    sys.stdout.write("Examples without an .rsftest file (%d):\n" % rsftest_no)
    if untested:
        rsftest_no_list.sort()
        for item in rsftest_no_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Examples with .rsftest file errors (%d):\n" % rsftest_errors)
    rsftest_error_list.sort()
    for item in rsftest_error_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Tested examples without a fig directory (%d):\n" % figdir_no)
    figdir_no_list.sort()
    for item in figdir_no_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Tested examples without a lock directory (%d):\n" % lockdir_no)
    lockdir_no_list.sort()
    for item in lockdir_no_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Tested examples with missing figs (%d):\n" % figs_missing)
    figs_missing_list.sort()
    for item in figs_missing_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Tested examples with extra figs (%d):\n" % figs_extra)
    figs_extra_list.sort()
    for item in figs_extra_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    sys.stdout.write("Tested examples with non-matching figs (%d):\n" % figs_diff)
    figs_diff_list.sort()
    for item in figs_diff_list: sys.stdout.write("%s\n" % item)
    sys.stdout.write("\n")

    return unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main())

# $Id: Mbooklist.py 4965 2009-11-05 03:27:38Z jennings_jim $
