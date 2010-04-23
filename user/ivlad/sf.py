#! /usr/bin/env python
'''
NAME
	sf
DESCRIPTION
	Wrappers for calling m8r programs elegantly and correctly from Py scripts
SOURCE
	user/ivlad/sf.py
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

import os, sys, math, rsfprog, string, random

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

# The exe argument allows the public functions below to perform in three ways:
# 1. called as sf.prog(arglist[,exe='x']) to execute directly a RSF program
# 2. called as sf.prog(arglist, exe='g') to return stdout of program
# 2. called as sf.prog(arglist, exe='p') to return a 'prog args' string
#    ready to be chained into a pipe and executed by ivlad.pp (no output 
#    captured)
# 3. called as sf.prog(arglist, exe='l') to return a 
#    'prog <inp >out args' string that can be written to a script, log, etc

################################################################################

def __run(prog, args, inp, out, verb, exe, postproc=None):
    
    prog = prog.strip()
    if args == None:
        args = ''
    else:
        args = args.strip()
    cmd  = prog

    if exe == 'g': # get output
        assert out == None
        out_str = ivlad.getout(prog, args, inp, verb)
        if postproc == None:
            return out_str
        else:
            return postproc(out_str)
    elif exe == 'p': # call part of a pipe
        if args != '':
            cmd += ' ' + args
        return cmd
    else:
        if inp != None:
            cmd += ' <' + inp.strip()
        if out != None:
            cmd += ' >' + out.strip()
        if args != '':
            cmd +=  ' ' + args
        if exe == 'l': # log
            return cmd
        else: # exe=='x', execute
            ivlad.exe(cmd, verb)

################################################################################
    
def __parse(arg_dict):
    'Turn the list of args into a string (default vals must be None!)'  
    
    args = ''
    for arg_nm in arg_dict.keys():
	if arg_nm not in ('inp', 'out', 'verb', 'exe'):
	    arg_v = arg_dict[arg_nm]
	    args += ivlad.switch(arg_v==None, '', ' %s=%s' % (arg_nm,arg_v))

    return args

################################################################################

def bandpass(inp=None, out=None, fhi=None, flo=None, nphi=None, nplo=None, 
    phase=None, verb=False, exe='x'):
    
    return __run('sfbandpass', __parse(locals()), inp, out, verb, exe)

################################################################################

def clip(inp=None, out=None, clip=None, verb=False, exe='x'):

    return __run('sfclip', __parse(locals()), inp, out, verb, exe)

################################################################################

def cp(inp, out, verb=False, exe='x'):
    return __run('sfcp', inp + ' ' + out, None, None, verb, exe)

################################################################################

def create(out=None, n=None, verb=False, exe='x'):

    args = ''
    nlist = ivlad.mklist(n)
    for i in range(len(nlist)):
        args += ' n%d=%d' % (i+1,nlist[i])
    return __run('sfcreate', args, None, out, verb, exe)

################################################################################

def csv2rsf(inp=None, out=None, delimiter=None, dtype=None, debug=None,
    trunc=None, o1=None, o2=None, d1=None, d2=None, unit1=None, unit2=None,
    label1=None, label2=None, verb=False, exe='x'):

    return __run('sfcsv2rsf', __parse(locals()), inp, out, verb, exe)

################################################################################

def filedims(inp=None, out=None, large=None, parform=None, verb=False, exe='x'):

    if exe == 'x' and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sffiledims', __parse(locals()), inp, out, verb, exe)

################################################################################

def fileflush(inp=None, out=None, verb=False, exe='x'):

    return __run('sffileflush', None, inp, out, verb, exe)

################################################################################

def get(inp=None, par=None, parform=False, out=None, verb=False, exe='x'):

    args = ['parform=' + ivlad.switch(parform, 'y', 'n')] + ivlad.mklist(par)
    if exe == 'x' and out==None: # invalid combination, fix the call
        exe = 'g'
    def postproc(out_str):
        out = out_str.split()
        if len(out) == 1:
            return out[0]
        else:
            return out
    return __run('sfget', ' '.join(args), inp, out, verb, exe, postproc)

################################################################################

def gettype(inp=None, out=None, verb=False, exe='x'):

    if exe == 'x' and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfgettype', __parse(locals()), inp, out, verb, exe)

################################################################################

def invalid(out=None, chk4nan=None, dir=None, rec=None, verb=False, exe='x'):

    if exe == 'x' and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfinvalid', __parse(locals()), out, verb, exe)

################################################################################

def leftsize(inp=None, out=None, i=None, verb=False, exe='x'):

    if exe == 'x' and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfleftsize', __parse(locals()), inp, out, verb, exe)

################################################################################

def real(inp=None, out=None, verb=False, exe='x'):

    return __run('sfreal', None, inp, out, verb, exe)

################################################################################

def remap1(inp=None, out=None, d1=None, n1=None, o1=None, order=None, 
    pattern=None, verb=False, exe='x'):

    return __run('sfremap1', __parse(locals()), inp, out, verb, exe)

################################################################################

def rm(files, verb=False, exe='x'):

    if type(files) == str:
        args = files
    else:
        args = ' '.join(file_list)
    return __run('sfrm', args, None, None, verb, exe)

################################################################################

def rtoc(inp=None, out=None, verb=False, exe='x'):

    return __run('sfrtoc', None, inp, out, verb, exe)

################################################################################

def transp(inp=None, out=None, plane=None, memsize=None, verb=False, exe='x'):

    return __run('sftransp', __parse(locals()), inp, out, verb, exe)

################################################################################

def wuab(inp=None, prog=None, tpar=None, ipar=None, verb=False, exe='x'):

    return __run('sfwuab', __parse(locals()), None, None, verb, exe)
