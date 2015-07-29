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

import os, sys, math, rsf.prog, string, random

try: # Give precedence to local version
    import ivlad, m8rex
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex

# The exe argument allows the public functions below to perform in three ways:
# 1. called as sf.prog(arglist[,exe=None]) to execute directly a RSF program
# 2. called as sf.prog(arglist, exe='g') to return stdout of program
# 2. called as sf.prog(arglist, exe='p') to return a 'prog args' string
#    ready to be chained into a pipe and executed by ivlad.pp (no output 
#    captured)
# 3. called as sf.prog(arglist, exe='l') to return a 
#    'prog <inp >out args' string that can be written to a script, log, etc
# Use sf._set_xmode to set exe to a certain value, to avoid having to specify it
# in every call
################################################################################

def __run(prog, args, inp, out, verb, exe, postproc=None):
    
    prog = prog.strip()
    if args == None:
        args = ''
    else:
        args = args.strip()
    cmd = os.path.join(rsf.prog.RSFROOT,'bin',prog)

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

glob_exe = 'x'

def _set_xmode(exec_mode=None):
    'Save programmer for typing exe= in every call'
    global glob_exe
    
    if exec_mode != None:
        ivlad.chk_par_in_list(exec_mode,['x','g','p','l'])
        glob_exe = exec_mode

def __x(exe, g_exe):
    if exe == None:
        return g_exe
    else:
        ivlad.chk_par_in_list(exe, ['x','g','p','l'])
        return exe

################################################################################
# Wrappers for sf programs
################################################################################

def add(out, files=None, scale=None, add=None, sqrt=None, abs=None, log=None,
    exp=None, mode=None, verb=False, exe=None):

    arg = __parse({'scale':scale, 'add':add, 'sqrt':sqrt, 'abs':abs, 'log':log,
        'exp':exp, 'mode':mode})
    # file names are not in key=val format, so __parse will not work
    filestr = ' '.join(files[1:])

    return __run('sfadd', arg+' '+filestr, 
        files[0], out,
        verb, __x(exe,glob_exe))    

################################################################################

def pad2nextfastsize(n, out=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'

    #def postproc(out_str):
    #    return int(out_str)

    return __run('sfpad2nextfastsize', __parse(locals()), 
        None, out, 
        verb, __x(exe,glob_exe), lambda x: int(x))

################################################################################

def attr(inp=None, out=None, lval=None, want=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'

    def postproc(out_str):
        if want in ('rms', 'min', 'max'):
            return float(out_str.split()[2])
        else:
            return out_str

    return __run('sfattr', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe), postproc)

################################################################################

def bandpass(inp=None, out=None, fhi=None, flo=None, nphi=None, nplo=None, 
    phase=None, verb=False, exe=None):
    
    return __run('sfbandpass', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def bar(inp=None, out=None, verb=False, exe=None):
  
    return __run('sfbar', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def byte(inp=None, out=None, allpos=None, bias=None, clip=None, gainpanel=None,
        pclip=None, verb=False, exe=None):

    return __run('sfbyte', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def cat(out, files=None, axis=None, o=None, d=None, verb=False, exe=None):

    aod = __parse({'axis':axis,'o':o,'d':d})
    # file names are not in key=val format, so __parse will not work
    filestr = ' '.join(files[1:])

    return __run('sfcat', aod+' '+filestr, 
        files[0], out,
        verb, __x(exe,glob_exe))    

################################################################################

def clip(inp=None, out=None, clip=None, verb=False, exe=None):

    return __run('sfclip', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def cp(inp, out, verb=False, exe=None):
        
    return __run('sfcp', inp + ' ' + out, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def create(out=None, n=None, verb=False, exe=None):

    args = ''
    nlist = ivlad.mklist(n)
    for i in range(len(nlist)):
        args += ' n%d=%d' % (i+1,nlist[i])
    return __run('sfcreate', args, 
        None, out, 
        verb, __x(exe,glob_exe))

################################################################################

def csv2rsf(inp=None, out=None, delimiter=None, dtype=None, debug=None,
    trunc=None, o1=None, o2=None, d1=None, d2=None, unit1=None, unit2=None,
    label1=None, label2=None, verb=False, exe=None):

    return __run('sfcsv2rsf', __parse(locals()), 
        inp,out,
        verb, __x(exe,glob_exe))

################################################################################

def filedims(inp=None, out=None, large=None, parform=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sffiledims', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def fileflush(inp=None, out=None, verb=False, exe=None):

    return __run('sffileflush', None, 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def get(inp=None, par=None, parform=False, out=None, verb=False, exe=None):

    args = ['parform=' + ivlad.switch(parform, 'y', 'n')] + ivlad.mklist(par)
    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'
    def postproc(out_str):
        out = out_str.split()
        if len(out) == 1:
            return out[0]
        else:
            return out
    return __run('sfget', ' '.join(args), 
        inp, out, 
        verb, __x(exe,glob_exe), postproc)

################################################################################

def gettype(inp=None, out=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfgettype', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def invalid(out=None, chk4nan=None, dir=None, rec=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfinvalid', __parse(locals()), 
        None, out, 
        verb, __x(exe,glob_exe))

################################################################################

def leftsize(inp=None, out=None, i=None, verb=False, exe=None):

    if exe == None and out==None: # invalid combination, fix the call
        exe = 'g'
    return __run('sfleftsize', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def mv(inp, out, verb=False, exe=None):
    return __run('sfmv', inp + ' ' + out, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def pclip(inp, out, pclip=None, verb=False, exe=None):

    arg_dict = locals()
    args = __parse(arg_dict) + ' inp='+inp + ' out=' + out
    return __run('sfpclip', args, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def real(inp=None, out=None, verb=False, exe=None):

    return __run('sfreal', None, 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def remap1(inp=None, out=None, d1=None, n1=None, o1=None, order=None, 
    pattern=None, verb=False, exe=None):

    return __run('sfremap1', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def rm(files, verb=False, exe=None):

    if type(files) == str:
        args = files
    else:
        args = ' '.join(files)
    return __run('sfrm', args, None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def rmrf(mydir, rec=False, verb=False, exe=None):

    arg = 'dir=' + mydir + ' rec=' + ivlad.switch(rec, 'y', 'n')

    return __run('sfrmrf', args, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def rtoc(inp=None, out=None, verb=False, exe=None):

    return __run('sfrtoc', None, 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def seekwin(inp, out=None, nread=None, nseek=None, whence=None, verb=False,
exe=None):
    # Takes stdin, but input cannot be from pipe

    return __run('sfseekwin', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def slant(inp=None, out=None, adj=False, rho=None, anti=None, np=None, dp=None,
p0=None, x0=None, dx=None, nx=None, p1=None, verb=False, exe=None):

    if adj:
        if np == None:
            raise m8rex.ConflictingArgs('np','None','adj','y')
        if p0 == None:
            raise m8rex.ConflictingArgs('p0','None','adj','y')
        if dp == None:
            raise m8rex.ConflictingArgs('dp','None','adj','y')
        adj = 'y'
    else:
        if nx == None:
            raise m8rex.ConflictingArgs('nx','None','adj','y')
        if x0 == None:
            raise m8rex.ConflictingArgs('x0','None','adj','y')
        if dx == None:
            raise m8rex.ConflictingArgs('dx','None','adj','y')
        adj = 'n'

    return __run('sfslant', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def split(inp, outdir=None, nthick=None, verb=False, exe=None):

    arg_dict = locals()
    args = __parse(arg_dict) + ' inp='+inp

    return __run('sfsplit', args, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def squeeze(inp=None, out=None, verb=False, exe=None):

    return __run('sfwindow', 'squeeze=y', 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def touch(mydir=None, rec=None, chk4nan=None, verb=False, exe=None):

    return __run('sftouch', __parse(locals()), 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def transp(inp=None, out=None, plane=None, memsize=None, verb=False, exe=None):

    return __run('sftransp', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def wiki2static(outdir=None, verb=False, exe=None):

    return __run('sfwiki2static', __parse(locals()), 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def window(inp=None, out=None, squeeze=None, j1=None, j2=None, j3=None, j4=None,
j5=None, j6=None, j7=None, j8=None, j9=None, d1=None, d2=None, d3=None, d4=None,
d5=None, d6=None, d7=None, d8=None, d9=None, f1=None, f2=None, f3=None, f4=None,
f5=None, f6=None, f7=None, f8=None, f9=None, min1=None, min2=None, min3=None,
min4=None, min5=None, min6=None, min7=None, min8=None, min9=None, n1=None, 
n2=None, n3=None, n4=None, n5=None, n6=None, n7=None, n8=None, n9=None, 
max1=None, max2=None, max3=None, max4=None, max5=None, max6=None, max7=None, 
max8=None, max9=None, verb=False, exe=None):

    return __run('sfwindow', __parse(locals()), 
        inp, out, 
        verb, __x(exe,glob_exe))

################################################################################

def wuab(inp, prog=None, tpar=None, ipar=None, verb=False, exe=None):

    arg_dict = locals()
    args = __parse(arg_dict) + ' inp='+inp

    return __run('sfwuab', args, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def ximage(inp, par=None, verb=False, exe=None):

    arg_dict = locals()
    args = 'inp='+inp + ivlad.switch(par==None, '', ' par="' + par + '"')

    return __run('sfximage', args, 
        None, None, 
        verb, __x(exe,glob_exe))

################################################################################

def zcp(inp, out, verb=False, exe=None):
    return __run('sfzcp', inp + ' ' + out, 
        None, None, 
        verb, __x(exe,glob_exe))
