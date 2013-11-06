#!/usr/bin/env python

# Copyright (C) 2008 University of Texas at Austin
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

import os, re, sys

################################################################################

def datapath():
    'Path for binary and temporary files'
    path = os.environ.get('DATAPATH')
    if not path:
        try:
            pathfile = open('.datapath','r')
        except:
            try:
                pathfile = open(os.path.join(os.environ.get('HOME'),
                                             '.datapath'),'r')
            except:
                pathfile = None
        if pathfile:
            for line in pathfile.readlines():
                check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],
                                 line)
                if check:
                    path = check.group(1)
            pathfile.close()
    if not path:
        path = './' # the ultimate fallback
    return path

################################################################################

def dirtree(cwd=None):
    'hierarcical directory structure'
    if not cwd:
        cwd = os.getcwd()
    tree = (os.path.basename(os.path.dirname(os.path.dirname(cwd))),
            os.path.basename(os.path.dirname(cwd)),
            os.path.basename(cwd))
    return tree

################################################################################

def mkdir(dir):
    'Recursive directory making'
    while os.path.basename(dir) == '.':
        dir = os.path.dirname(dir)
    if dir and not os.path.isdir(dir):
        mkdir(os.path.dirname(dir))
        os.mkdir(dir)
    return dir

################################################################################

def sconsign(env):
    'SConsign database file'

    try: 	 
        import dbhash 	 
        env.SConsignFile(env.path+ '.sconsign.dbhash',dbhash)	 
    except: 	 
        try: 	 
            import gdbm 	 
            env.SConsignFile(env.path+ '.sconsign.gdbm',gdbm)
        except:
            env.SConsignFile(env.path+ '.sconsign')

################################################################################

def getpath(cwd):
    top = datapath()
    path = os.path.dirname(top)
    if top[:2] != './':
        # create a hierarchical structure
        tree = dirtree(cwd)
        for level in tree:
            if level:
                path = os.path.join(path,level)
        mkdir(path)
    path = os.path.join(path,os.path.basename(top))
    return path

