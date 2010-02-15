#! /usr/bin/env python
'''
NAME
	m8rex
DESCRIPTION
	Collection of user-defined exceptions
SOURCE
	user/ivlad/m8rex.py
'''
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

import os

########################################

class Error(Exception):
    'Base class for exceptions in this module.'
    def __str__(self):
        return self.msg

########################################

# Exceptions ordered alphabetically

class MissingArgument(Error):
    'For checking command-line arguments to Python main programs'
    def __init__(self, arg):
        self.msg = 'Missing argument: ' + arg

class NotAValidDir(Error):
    def __init__(self, dirnm):
        self.msg = 'Not a valid directory: ' + dirnm

class NoXPermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +x permissions for ' + os.getlogin()

class NoReadPermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +r permissions for ' + os.getlogin()

class NoWritePermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +w permissions for ' + os.getlogin()

class StringParamNotInAcceptableValueList(Error):
    def __init__(self, param, avl):
        self.msg = 'Parameter ' + param + ' not in: ' + ', '.join(avl)


