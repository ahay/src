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

def bandpass(inp=None, out=None, fhi=None, flo=None, nphi=None, nplo=None, 
    phase=None, verb=False, exe='x'):
    fhi_str  = ivlad.switch( fhi==None,   '', ' fhi='  +str(fhi ))
    flo_str  = ivlad.switch( flo==None,   '', ' flo='  +str(flo ))
    nphi_str = ivlad.switch( nphi==None,  '', ' nphi=' +str(nphi))
    nplo_str = ivlad.switch( nplo==None,  '', ' nplo=' +str(nplo))
    phase_str = ivlad.switch(phase==None, '', ' phase='+str(phase))
    verb_str = ivlad.switch(verb, ' verb=y', '')
    args = fhi_str + flo_str + nphi_str + nplo_str + phase_str + verb_str
    return __run('sfbandpass', args, inp, out, verb, exe)

################################################################################

def clip(inp=None, out=None, clip=99, verb=False, exe='x'):
    return __run('sfclip', 'clip=' + str(clip), inp, out, verb, exe)

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

def real(inp=None, out=None, verb=False, exe='x'):
    return __run('sfreal', None, inp, out, verb, exe)

################################################################################

def remap1(inp=None, out=None, d1=None, n1=None, o1=None, order=None, 
    pattern=None, verb=False, exe='x'):
    d1_str = ivlad.switch(d1==None, '', ' d1='+str(d1))
    n1_str = ivlad.switch(n1==None, '', ' n1='+str(n1))
    o1_str = ivlad.switch(o1==None, '', ' o1='+str(o1))
    pattern_str = ivlad.switch(pattern==None, '', ' pattern='+pattern)
    args = d1_str + n1_str + o1_str + pattern_str
    return __run('sfremap1', args, inp, out, verb, exe)

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

def transp(inp=None, out=None, plane=12, memsize=None, verb=False, exe='x'):
    memsz_str = ivlad.switch(memsize==None, '', 'memsize=' + str(memsize))
    plane_str = ' plane=' + str(plane)
    args = memsz_str + plane_str
    return __run('sftransp', args, inp, out, verb, exe)

################################################################################

def wuab(inp=None, prog=None, tpar=None, ipar=None, verb=False, exe='x'):
    vflag = ivlad.switch(verb, 'y', 'n')
    args = 'inp=%s prog=%s tpar="%s" ipar="%s" verb=%s' %\
        (inp, prog, tpar, ipar, vflag)
    return __run('sfwuab', args, None, None, verb, exe)
