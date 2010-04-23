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

################################################################################

def bandpass(inp, out, fhi=None, flo=None, nphi=None, nplo=None, phase=None, 
    verb=False):
    fhi_str  = ivlad.switch(fhi==None,  '', ' fhi=' +str(fhi ))
    flo_str  = ivlad.switch(flo==None,  '', ' flo=' +str(flo ))
    nphi_str = ivlad.switch(nphi==None, '', ' nphi='+str(nphi))
    nplo_str = ivlad.switch(nplo==None, '', ' nplo='+str(nplo))
    phase_str = ivlad.switch(phase==None, '', ' phase='+phase)
    verb_str = ivlad.switch(verb, ' verb=y', '')
    cmd = 'sfbandpass <%s >%s' % (inp,out)
    cmd += fhi_str + flo_str + nphi_str + nplo_str + phase_str + verb_str
    ivlad.exe(cmd, verb)

def clip(inp, out, clip, verb=False):
    ivlad.exe('sfclip <%s clip=%s >%s' % (inp, clip, out), verb)
	
def cp(inp, out, verb=False):
    ivlad.exe('sfcp %s %s' % (inp, out), verb)

def create(out, n, verb=False):
    cmd = 'sfcreate >' + out
    nlist = ivlad.mklist(n)
    for i in range(len(nlist)):
        cmd += ' n%d=%d' % (i+1,nlist[i])
    ivlad.exe(cmd, verb)

def get(inp, par, parform=False, verb=False):
    args = ['parform=' + ivlad.switch(parform, 'y', 'n')] + ivlad.mklist(par)
    out = ivlad.getout('sfget', args, inp, verb).split()
    if len(out) == 1:
        return out[0]
    else:
        return out

def real(inp, out, verb=False):
    ivlad.exe('sfreal <%s >%s' % (inp, out), verb)

def remap1(inp, out, d1=None, n1=None, o1=None, order=None, pattern=None, 
    verb=False):
    d1_str = ivlad.switch(d1==None, '', ' d1='+str(d1))
    n1_str = ivlad.switch(n1==None, '', ' n1='+str(n1))
    o1_str = ivlad.switch(o1==None, '', ' o1='+str(o1))
    pattern_str = ivlad.switch(pattern==None, '', ' pattern='+pattern)
    cmd = 'sfremap1 <%s >%s' % (inp, out) 
    cmd += d1_str + n1_str + o1_str + pattern_str
    ivlad.exe(cmd, verb)

def rm(files, verb=False):
    if type(files) == str:
        f_str = files
    else:
        f_str = ' '.join(file_list)
    ivlad.exe('sfrm ' + f_str, verb)

def rtoc(inp, out, verb=False):
    ivlad.exe('sfrtoc <%s >%s' % (inp, out), verb)

def transp(inp, out, plane=12, memsize=None, verb=False):
    memsz_str = ivlad.switch(memsize==None, '', 'memsize=' + str(memsize))
    plane_str = ' plane=' + str(plane)
    cmd = 'sftransp <%s >%s' % (inp, out)
    ivlad.exe(cmd + memsz_str + plane_str, verb)

def wuab(inp, prog, tpar, ipar, verb=False):
    vflag = ivlad.switch(verb, 'y', 'n')
    cmd = 'sfwuab inp=%s prog=%s tpar="%s" ipar="%s" verb=%s' %\
        (inp, prog, tpar, ipar, vflag)
    ivlad.exe(cmd, verb)

