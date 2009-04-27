#!/usr/bin/env python
'''Compare Vplot files in Fig and Lock directories
Parameter figdir is path to Fig directory, default is ./Fig.
Parameter lockdir is path to Lock directory:
    If figdir is in $RSFSRC/book/[book]/[chapter]/[section],
        then default lockdir is $RSFFIGS/[book]/[chapter]/[section].
    If figdir is not in $RSFSRC/book/[book]/[chapter]/[section],
        then default lockdir is $RSFALTFIGS/[book]/[chapter]/[section].

Parameter list controls files to list, default is all.
Parameter show controls files to flip with xtpen, default is none.

list|show = none    No files, print only summary.
list|show = diff    Files that are different, determined by sfvplotdiff.
list|show = miss    Files missing from figdir or lockdir, and different files.
list|show = all     All files.

File list codes:

space   indicates files that are the same.
  -     indicates file in lockdir that is missing from figdir.
  +     indicates extra file in figdir that is missing from lockdir.
number  is return code from sfvplotdiff indicating different files.'''

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

import os, copy, sys, rsfprog

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf
    
sfprefix = 'sf'                 # prefix for rsf commands
plprefix = 'vp'                 # prefix for vpl commands
sfsuffix = '.rsf'               # suffix for rsf files
vpsuffix = '.vpl'               # suffix for vplot files

###############################################################################

def book_path_split(figdir):
    '''
    split rsf fig path into book location, book tree, and fig directory
    '''
    
    tuple4 = os.path.split(figdir)
    tuple3 = os.path.split(tuple4[0])
    tuple2 = os.path.split(tuple3[0])
    tuple1 = os.path.split(tuple2[0])
    
    fig_dir       = tuple4[1]
    book_location = tuple1[0]
    book_tree     = os.path.join(tuple1[1],tuple2[1],tuple3[1])

    return (book_location,book_tree,fig_dir)
    
def vpl_list(list):
    '''
    remove non-vplot filenames from a list
    '''
        
    newlist = []
    for item in list:
        if os.path.splitext(item)[1] == vpsuffix: newlist.append(item)
        
    return newlist

###############################################################################

def main(argv=sys.argv):

    unix_success = 0
    unix_error   = 1

################    get user parameters
    
    if len(argv) < 2 :              # print selfdoc and exit if no parameters
        rsfprog.selfdoc()
        return unix_error

    par     = rsf.Par(argv)         # get parameters
    figdir  = par.string('figdir')  # fig directory, default = ./Fig
    lockdir = par.string('lockdir') # lock directory, default = lock counterpart of figdir
    list    = par.string('list')    # how much to list [none,diff,miss,all], default = all
    show    = par.string('show')    # how much to show [none,diff,miss,all], default = none
        
                                    # check list and show parameters
    options = ['none','diff','miss','all']
    if list is None:    list='all'
    if options.count(list) == 0:
        print "Unknown list option: %s" % list
        return unix_error
    if show is None:    show='none'
    if options.count(show) == 0:
        print "Unknown show option: %s" % show
        return unix_error

################    get environment variables
    
                                    # get $RSFFIGS variable
    rsffigs = os.environ.get('RSFFIGS')
    
                                    # use $RSFROOT/figs if $RSFFIGS not defined
    if rsffigs == None:
        if os.environ.get('RSFROOT') == None:
            print "Environment variable $RSFROOT not defined."
            return unix_error
        else:
            rsffigs = os.path.expandvars(os.path.join('$RSFROOT','figs'))
    
                                    # get $RSFALTFIGS variable
    rsfaltfigs = os.environ.get('RSFALTFIGS')
    
                                    # use rsffigs if $RSFALTFIGS not defined
    if rsfaltfigs == None: rsfaltfigs = rsffigs

################    get directory paths

                                    # get fig directory path
    if figdir is None:  figpath = os.path.abspath(os.path.join('.','Fig'))
    else:               figpath = os.path.abspath(figdir)
            
                                    # get lock directory path
    if lockdir is None:
                                    # default to $RSFFIGS/book_tree
        book_paths = book_path_split(figpath)
        lockpath = os.path.expandvars(os.path.join(rsffigs,book_paths[1]))
        
                                    # if $RSFSRC is defined and we are not in
                                    # $RSFSRC/book/book_tree then
                                    # default to $RSFALTFIGS/book_tree
        if os.environ.get('RSFSRC') != None:
            if book_paths[0] != os.path.expandvars(os.path.join('$RSFSRC','book')):
                lockpath = os.path.expandvars(os.path.join(rsfaltfigs,book_paths[1]))
                    
                                    # use user-supplied lockdir
    else:
        lockpath = os.path.abspath(lockdir)
    
                                    # check if fig directory exists
    print ""
    if not os.path.exists(figpath):
        print "Fig directory does not exist:"
        print figpath
        return unix_error
    else:
        print "Fig directory:"
        print figpath
    
                                    # check if lock directory exists
    print ""
    if not os.path.exists(lockpath):
        print "Lock directory does not exist:"
        print lockpath
        return unix_error
    else:    
        print "Lock directory:"
        print lockpath

################    get file lists

                                    # get lists of vpl files    
    figlist  = vpl_list(os.listdir(figpath))
    locklist = vpl_list(os.listdir(lockpath))
    
    figlist.sort()
    locklist.sort()
    
    files = {}                      # merge vpl lists
    for item in figlist:  files[item] = '  '
    for item in locklist: files[item] = '  '
    filelist = files.keys()
    filelist.sort()
    
################    find missing and different files

                                    # find missing files
    n_fig_missing = 0
    n_fig_extra   = 0
    for item in filelist:
        if figlist.count(item)  == 0: 
            files[item]   = ' -'
            n_fig_missing = n_fig_missing+1
    for item in filelist:
        if locklist.count(item) == 0: 
            files[item]   = ' +'
            n_fig_extra   = n_fig_extra+1

                                    # find different files
    n_different = 0
    binpath = os.path.expandvars(os.path.join('$RSFROOT','bin'))
    command = os.path.join(binpath,sfprefix+'vplotdiff')
    for item in filelist:
        if files[item] == '  ':
            figfile  = os.path.join(figpath,item)
            lockfile = os.path.join(lockpath,item)        
            diff     = os.system(' '.join([command,figfile,lockfile,'&> /dev/null']))
            if diff != 0: 
                files[item] = '%2d' % (diff/256)
                n_different = n_different+1
    
################    print file list and show selected figs

    n_same = 0
    print ""
    command = os.path.join(binpath,'xtpen')
    for item in filelist:
        miss_check =     (files[item] != '  ')

        diff_check = (    (files[item] != '  ')
                      and (files[item] != ' -')
                      and (files[item] != ' +'))
                     
        list_check = (   ((list == 'all'))
                      or ((list == 'miss') and miss_check)
                      or ((list == 'diff') and diff_check))
                     
        show_check = (   ((show == 'all'))
                      or ((show == 'miss') and miss_check)
                      or ((show == 'diff') and diff_check))
        
        if files[item] == ' -': figfile  = ''
        else:                   figfile  = os.path.join(figpath,item)
        
        if files[item] == ' +': lockfile = ''
        else:                   lockfile = os.path.join(lockpath,item)
        
        if list_check: print " %s %s" % (files[item],item)
        if show_check: os.system(' '.join([command,figfile,lockfile]))
        if files[item] == '  ': n_same = n_same+1
        
    print ""
    print "Identical files:         %3d" % n_same
    print "Different files:         %3d" % n_different
    print "Files missing from Fig:  %3d" % n_fig_missing
    print "Extra files in Fig:      %3d" % n_fig_extra
    print "Total vplot files:       %3d" % len(filelist)
    print ""
   
    return unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main())

# $Id: Mminmax.c 4378 2009-04-25 04:30:26Z jennings_jim
