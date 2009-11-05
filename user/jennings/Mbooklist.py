#!/usr/bin/env python
'''List properties of Madagascar example book directories.
Scan a directory tree (or list of directory trees) to a given depth,
inventory the contents, and optionally execute a command in each leaf.

The inventory and optional command occurs only at the specified depth
(default 3), not the intervening depths.  Directories named .svn are
skipped.  Only directories containing an SConstruct file are listed or
executed.

Optional directory filters controlling inventory or command execution
may be specified based on existence of .rsfproj file, sf programs used,
type of external data required, and total rsf data-file size of the
completed example.

The optional command is executed in the users default shell.

Examples (from within $RSFSRC):

sfbooklist book                         # inventory of book
sfbooklist levels=2 book/geostats       # inventory of book/geostats
sfbooklist command=scons book           # build examples with default filters
sfbooklist size=5 command=scons book    # build examples smaller than 5MB
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

import os, copy, sys, signal, rsfprog

if not hasattr(os,'walk'):
    sys.stderr.write('booklist needs Python 2.3 or greater\n')
    sys.exit(unix_error)

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf

sfprefix = 'sf'                 # prefix for rsf commands
plprefix = 'vp'                 # prefix for vpl commands
sfsuffix = '.rsf'               # suffix for rsf files
vpsuffix = '.vpl'               # suffix for vplot files

###############################################################################

def handler(signum, frame):
    'signal handler for abortion [Ctrl-C]'

    sys.stderr.write('\n[Ctrl-C] Aborting...\n')
    if child:
        os.kill (child,signal.SIGINT)
    sys.exit(-1)

signal.signal(signal.SIGINT,handler) # handle interrupt

child = None
def command_wait(comm):
    'Interruptable system command'

    global child
    child = os.fork()
    if child:
        t0 = os.times()
        (pid,exit) = os.waitpid(child,0)
        t1 = os.times()
        child = 0
        dt_user = t1[2]-t0[2]
        dt_sys  = t1[3]-t0[3]
        dt_real = t1[4]-t0[4]
        return (exit,dt_user,dt_sys,dt_real)
    else:
        status = os.system(comm)
        if (status==0): os._exit(0)
        else:           os._exit(1)

def size_string(size):
    'Make a human-readable size string'

    units = [' B','KB','MB','GB','TB']
    p = 1
    while 1024**p <= size : p = p+1
    return "%4.2f %s" % (float(size)/1024**(p-1),units[p-1])

def read_rsfproj(root,files):
    'Read contents of an .rsfproj file'

    error           = 0             # defaults
    rsfproj_exist   = 'no'
    uses_list       = []
    data_type       = 'unknown'
    data_size       = 0

    if '.rsfproj' in files:         # get info from .rsfproj
        rsfproj_exist = 'yes'

        g = {}                      # read file
        l = {}
        execfile(os.path.join(root,'.rsfproj'),g,l)

                                    # get uses
        if 'uses' in l: uses_list = l['uses']
        else:           error = 1

                                    # get data type
        if 'data' in l:
            if len(l['data']) == 0:     data_type = 'none'
            if len(l['data']) >  0:     data_type = 'public'
            if 'PRIVATE' in l['data']:  data_type = 'private'
        else:
            error = 1

                                    # get data size
        if 'size' in l: data_size = l['size']
        else:           error = 1

    return (error,rsfproj_exist,uses_list,data_type,data_size)

def calc_filter(options,props):
    'Calculate command filter'

    filter = True
    (rsfproj,uses,size,fetch_none,fetch_public,fetch_private) = options
    (rsfproj_exist,uses_list,data_size,data_type) = props

    # rsfproj existence filter
    if (rsfproj != 'both') and (rsfproj != rsfproj_exist): filter = False

    # uses filter
    if (uses != 'any') and (uses not in uses_list): filter = False

    # total rsf data file size filter
    if (data_size > size): filter = False

    # external data type filter
    if (fetch_none    == False) and (data_type == 'none'):    filter = False
    if (fetch_public  == False) and (data_type == 'public'):  filter = False
    if (fetch_private == False) and (data_type == 'private'): filter = False

    return filter

###############################################################################

def main(argv=sys.argv):

    unix_success = 0
    unix_error   = 1

################    get user parameters

    if len(argv) < 2 :              # print selfdoc and exit if no parameters
        rsfprog.selfdoc()
        return unix_error

    par = rsf.Par(argv)             # get parameters

    levels = par.int('levels',3)
    # directory search depth

    list = par.string('list')
    # how much to list [all,filter,none], default = all

    rsfproj = par.string('rsfproj')
    # rsfproj filter [yes,no,both], default = yes

    size = par.int('size',1024**2)
    # max data size filter (MB)
    size = size*1024**2

    uses = par.string('uses')
    # uses filter, default = any

    fetch_none = par.bool('nofetch',True)
    # fetch-no-data filter

    fetch_public = par.bool('public',False)
    # fetch-public-data filter

    fetch_private = par.bool('private',False)
    # fetch-private-data filter

    command = par.string('command')
    # command to execute in each directory, default = none

    if list is None:    list='all'
    if list not in ['all','filter','none']:
        sys.stderr.write('Unknown list option: %s\n' % list)
        return unix_error

    if rsfproj is None: rsfproj='yes'
    if rsfproj not in ['yes','no','both']:
        sys.stderr.write('Unknown rsfproj option: %s\n' % rsfproj)
        return unix_error

    if uses is None:    uses='any'

    if fetch_none is None:      fetch_none    = True
    if fetch_public is None:    fetch_public  = False
    if fetch_private is None:   fetch_private = False

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
    sys.stdout.flush()

################    search directory tree

    if (list == 'all') or (list == 'filter'):
        sys.stdout.write('command   rsfproj     size    data        directory\n')
        sys.stdout.flush()

    total_list      = 0
    total_size      = 0
    pass_list       = 0
    pass_size       = 0
    rsfproj_yes     = 0
    rsfproj_no      = 0
    rsfproj_error   = 0
    command_error   = 0
    data_unknown    = 0
    data_none       = 0
    data_public     = 0
    data_private    = 0

    for bookdir in books:
        for root,dirs,files in os.walk(bookdir):

            if '.svn' in dirs: dirs.remove('.svn')      # don't visit .svn dirs

            reldir = root.replace(bookdir,"dummy",1)    # get relative dir name
            level  = len(reldir.split(os.sep))-1        # get dir level

            if level <  levels: continue        # scan to depth 'levels' ...
            if level == levels: del dirs[:]     # ... but not beyond 'levels'

                                                # skip dirs without SConstruct
            if 'SConstruct' not in files: continue

                                                # read rsfproj file
            tuple = read_rsfproj(root,files)
            (error,rsfproj_exist,uses_list,data_type,data_size) = tuple
            if error==1:
                rsfproj_error = rsfproj_error+1
                string = "   *********  .rsfproj error   *********  %s\n"
                sys.stdout.write(string % root)
                sys.stdout.flush()

                                                # calculate directory filter
            options = (rsfproj,uses,size,fetch_none,fetch_public,fetch_private)
            props   = (rsfproj_exist,uses_list,data_size,data_type)
            filter  = calc_filter(options,props)
            if filter==True:
                pass_list = pass_list+1
                pass_size = pass_size+data_size

                                        # print data for each directory
            if (list == 'all') or ((list == 'filter') and (filter == True)):

                total_list = total_list+1
                total_size = total_size+data_size

                if filter:  filter_command = 'yes'
                else:       filter_command = 'no '

                if rsfproj_exist == 'yes' : rsfproj_yes   = rsfproj_yes+1
                if rsfproj_exist == 'no'  : rsfproj_no    = rsfproj_no+1

                if data_type == 'unknown' : data_unknown = data_unknown+1
                if data_type == 'none'    : data_none    = data_none+1
                if data_type == 'public'  : data_public  = data_public+1
                if data_type == 'private' : data_private = data_private+1

                tuple = (filter_command,rsfproj_exist,
                         size_string(data_size),data_type,root)
                sys.stdout.write('%s       %-3s     %10s  %-7s     %s\n' % tuple)
                sys.stdout.flush()

                                        # execute command in directory
            if (filter is True) and (command is not None):
                string = "   +++++++++  running command  +++++++++  %s\n"
                sys.stdout.write(string % root)
                sys.stdout.flush()

                tuple = command_wait(' '.join(['cd',root,';',command]))
                (exit,dt_user,dt_sys,dt_real) = tuple

                string = "   user %6.2f   sys %6.2f  real %6.2f  %s\n"
                sys.stdout.write(string % (dt_user,dt_sys,dt_real,root))
                if (exit==0):
                    string = "   ---------  command success  ---------  %s\n"
                else:
                    string = "   *********   command error   *********  %s\n"
                    command_error = command_error+1
                sys.stdout.write(string % root)
                sys.stdout.write("\n")
                sys.stdout.flush()

    sys.stdout.write("\n")
    sys.stdout.write("Directories listed : %3d\n" % total_list)
    sys.stdout.write("Total data size    : %s\n" % size_string(total_size))
    sys.stdout.write("\n")
    sys.stdout.write("Directories with .rsfproj    : %3d\n" % rsfproj_yes)
    sys.stdout.write("Directories without .rsfproj : %3d\n" % rsfproj_no)
    sys.stdout.write(".rsfproj errors              : %3d\n" % rsfproj_error)
    sys.stdout.write("\n")
    sys.stdout.write("Directories for command     : %3d\n" % pass_list)
    sys.stdout.write("Total data size for command : %s\n"  % size_string(pass_size))
    sys.stdout.write("Command errors              : %3d\n" % command_error)
    sys.stdout.write("\n")
    sys.stdout.write("Directories using unknown external data : %3d\n" % data_unknown)
    sys.stdout.write("Directories using no external data      : %3d\n" % data_none)
    sys.stdout.write("Directories using public external data  : %3d\n" % data_public)
    sys.stdout.write("Directories using private external data : %3d\n" % data_private)
    sys.stdout.write("\n")

    return unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main())

# $Id$
