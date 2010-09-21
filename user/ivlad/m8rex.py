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
# Message should specify what the problem is

class ConflictingArgs(Error):
    def __init__(self, arg1, val1, arg2, val2):
        self.msg = '%s cannot be %s when %s=%s' % \
        (arg1, str(val1), arg2, str(val2))

class FailedExtCall(Error):
    def __init__(self, call):
        self.msg = 'Failed external call:\n' + call

class FailedWrite(Error):
    def __init__(self, fname):
        self.msg = 'Failed writing to file %s. Possibly disk full.' % fname

class MissingArgument(Error):
    def __init__(self, arg):
        self.msg = 'Missing argument: ' + arg

class MissingProgram(Error):
    def __init__(self, prog):
        self.msg = 'Missing executable: ' + prog

class NdimsMismatch(Error):
    def __init__(self, filenm, ndims):
        self.msg = 'File ' + filenm + ' must be %d-D' % ndims

class NotAValidDir(Error):
    def __init__(self, dirnm):
        self.msg = 'Not a valid directory: ' + dirnm

class NotAValidFile(Error):
    def __init__(self, filenm):
        self.msg = 'Not a valid file: ' + filenm

class NoXPermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +x permissions for ' + os.getlogin()

class NoReadPermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +r permissions for ' + os.getlogin()

class NoReturnFromExtProgram(Error):
    def __init__(self, prognm):
        self.msg = 'No return from program ' + prognm

class NoWritePermissions(Error):
    def __init__(self, name):
        self.msg = name + ' lacks +w permissions for ' + os.getlogin()

class ParBeyondLimit(Error):
    def __init__(self, parnm, limval, comp):
        self.msg = parnm + ' must be ' + comp + str(limval)

class ParamOutOfRange(Error):
    def __init__(self, param, minval, maxval):
        self.msg = 'Parameter %s not between %s and %s' % \
                (param, str(minval), str(maxval))

class StringParamNotInAcceptableValueList(Error):
    def __init__(self, param, avl):
        self.msg = 'Parameter ' + param + ' not in: ' + ', '.join(avl)

class StringParamInvalidFormat(Error):
    def __init__(self, param, msg):
        self.msg = 'Parameter ' + param + ' ' + msg

class TypeHandlingNotImplemented(Error):
    def __init__(self, typenm):
        self.msg = 'Handling for type %s not implemented yet' % typenm

class WrongPath(Error):
    def __init__(self, path):
        self.msg = path + ' does not exist'
