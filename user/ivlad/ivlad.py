#! /usr/bin/env python
"""
NAME
	ivlad
DESCRIPTION
	Utilities for python metaprograms in RSFSRC/user/ivlad
SOURCE
	user/ivlad/ivlad.py
"""
# Copyright (C) 2007-2010 Ioan Vlad
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

import os, sys, math, string, random
import rsf.prog as rsfprog

try: # Give precedence to local version
    import m8rex
except: # Use distributed version
    import rsf.user.m8rex as m8rex

import rsf.api as rsf

if sys.version_info[0] >= 3:
    from subprocess import getoutput
else:
    from commands import getoutput

try:
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False
    import commands

# Operating system return codes
unix_success = 0
unix_error   = 1

# Other constants
pipe = ' | '
ext  = '.rsf'
fmt = {'int':'i', 'float':'f'} # For usage with struct.pack

# Horizontal ruler (for logs, screen messages):
hr = '#' + 79 * '-'

################################################################################

def pp(call_list, inp=None, out=None):
    'Print a pipe command. Output can be fed to ivlad.exe'
    # Example:
    # ivlad.pp([sf.rtoc(exe='p'),
    #           sf.real(exe='p'),
    #           sf.get(par='n1',exe='p')],
    #          'junk.rsf','file.asc')
    # Outputs:
    # '<junk.rsf sfrtoc | sfreal | sfget parform=n n1'

    cmd = ''
    if inp != None:
        cmd += '<' + inp + ' '
    cmd += pipe.join(call_list)
    if out != None:
        cmd += ' >' + out
    return cmd

################################################################################

def mklist(x):
    if type(x) == str:
        return x.split()
    elif type(x) == list:
        return x
    elif type(x) == int:
        return [x]

################################################################################

def append(string, logfile):
    'Equivalent of echo string >>logfile'
    newline = '\n'
    if string[-1] != newline:
        string += newline
    handle = open(logfile.strip(), 'a')
    handle.write(string)
    handle.close()

################################################################################

def show_man_and_out(condition):
    'Display self-doc (man page) and exit'

    if condition:
        rsfprog.selfdoc() # show the man page
        sys.exit(unix_error)

################################################################################

def run(main_func, cpar=[], nminarg=None, nmaxarg=None):
    'For eliminating boilerplate code in main programs'

    # main_func must take a single argument -- par
    # cpar is compulsory parameter list
    # nminarg: minimum number of command-line arguments to the program
    # nmaxarg: max nr of CLI args to prog

    if nminarg != None:
        show_man_and_out(len(sys.argv) < nminarg+1)
    if nmaxarg != None:
        show_man_and_out(len(sys.argv) > nmaxarg+1)
    par = rsf.Par(sys.argv) # Parse arguments into a parameter table
    show_man_and_out(par.bool('help', False))
    for p in cpar:
        q = par.string(p)
        show_man_and_out(q == None)
    try:
        status = main_func(par) # run the program
    except m8rex.Error as e:
        msg(e.msg, True)
        status = unix_error

    sys.exit(status)

################################################################################

def exe(cmd, verb=False, dryrun=False, logfile=None):
    'Echoes a string command to screen via stderr, then executes it'

    if type(cmd) == str:
        cmd_list = [cmd]
    elif type(cmd) == list:
        cmd_list = cmd
    else:
        raise

    if logfile != None:
        assert type(logfile) == str

    for c in cmd_list:
        assert type(c) == str
        msg(c, verb)

    if have_subprocess:
        system_call = lambda x: subprocess.call(x, shell=True)
    else:
        system_call = lambda x: os.system(x)

    if not dryrun:
        for c in cmd_list:
            if logfile != None:
                append(c,logfile)
            system_call(c)

################################################################################

def which(prog):
    'Same functionality as Unix which'

    def __is_x(fname):
        return os.path.exists(fname) and os.access(fname, os.X_OK)

    fpath, fname = os.path.split(prog)

    if fpath == None or fpath.strip()=='':
        for path in os.environ["PATH"].split(os.pathsep):
            x_file = os.path.join(path, prog)
            if __is_x(x_file):
                return x_file
    else:
        if __is_x(prog):
            return prog

################################################################################

def chk_prog(prog):
    if which(prog) == None:
        raise m8rex.MissingProgram(prog)

################################################################################

def get1st(file, type='float'):
    'Get the first element in file'

    out_str = getppout(['sfwindow','sfdisfil'],
                       [['n1=1','n2=1'],['number=n']],
                       file)

    out_str = out_str.strip()

    if type == 'float':
        return float(out_str)
    elif type == 'int':
        return int(out_str)

################################################################################

def getppout(proglist, arglist=[[]], stdin=None):
    'Get pipe output'

    if stdin == None:
        finp = None
    else:
        if os.path.isfile(stdin):
            finp = open(stdin,'r')
        else:
            raise m8rex.NotAValidFile(stdin)

    c = [proglist[0]] + arglist[0]
    clist = [subprocess.Popen(c, stdin=finp, stdout=subprocess.PIPE)]

    for i in range(1,len(proglist)):
        c = [proglist[i]]+arglist[i]
        clist.append(subprocess.Popen(c,
                                      stdin=clist[i-1].stdout,
                                      stdout=subprocess.PIPE))

    return clist[-1].communicate()[0]

################################################################################

def getout(prog, arg=None, stdin=None, verb=False, raiseIfNoneOut=False):
    '''Replacement for commands.getoutput. Arguments:
    - prog. Executable to be run. STRING. The only non-optional argument.
    - arg. List of strings with program arguments, or just a string
    - stdin. Filename STRING.
    - verb: whether to display command before executing it. BOOL
    Returned value: stdout of command, with stripped newlines'''

    assert type(prog) == str
    assert type(verb) == bool
    if stdin != None:
        stdin = stdin.strip()

    def cat_cmd(cmdlist, stdin):
        cmd = ' '.join(cmdlist)
        if stdin != None:
            cmd += ' <' + stdin
        return cmd

    prog = os.path.join(rsfprog.RSFROOT,'bin',prog)

    #if which(prog) == None:
    #   raise m8rex.MissingProgram(prog)

    # Build the [prog, args] list
    cmdlist = [prog]
    if arg != None:
        cmdlist += mklist(arg)

    # Build command string for printing or Python < 2.4
    if verb or not have_subprocess:
        cmd = cat_cmd(cmdlist, stdin)
        msg(cmd, verb)

    if have_subprocess:
        if stdin == None:
            finp = None
        else:
            if os.path.isfile(stdin):
                finp = open(stdin,'r')
            else:
                raise m8rex.NotAValidFile(stdin)
        try:
            s = subprocess.Popen(cmdlist,stdin=finp,stdout=subprocess.PIPE)
            output = s.communicate()[0]
        except:
            raise m8rex.FailedExtCall(cat_cmd(cmdlist, stdin))
    else: # no subprocess module present
        try:
            output = getoutput(cmd)
        except:
            raise m8rex.FailedExtCall(cmd)

    if output == None:
        if raiseIfNoneOut:
            raise m8rex.NoReturnFromExtProgram(prog)
    else:
        output = output.rstrip('\n')

    return output

################################################################################

def getcataxis(filenm):
    'Returns sffiledims(filenm)+1'

    ndims_str = getout('sffiledims', 'parform=n', filenm)
    nds = int(ndims_str.split(':')[0])
    return str(nds+1)

################################################################################

def readaxis( inp, axisnr, verb=False ):
    '''Reads n,o,d for one axis of a Madagascar hypercube
    - inp: filename. STRING.
    - axisnr: INT
    - verb: BOOL'''

    ax = str(axisnr)

    n = getout('sfget',['parform=n','n'+ax], inp, verb)
    o = getout('sfget',['parform=n','o'+ax], inp, verb)
    d = getout('sfget',['parform=n','d'+ax], inp, verb)

    return( int(n), float(o), float(d))

################################################################################

def ndims(filename):
    'Returns nr of dims in a file, ignoring trailing length 1 dimensions'

    max_dims = 9
    nlist = [1] * max_dims

    for dim in range(1,max_dims+1):
        sfget_out = getout('sfget',['parform=n','n'+str(dim)],filename)
        if len(sfget_out) > 0:
            curr_n = int(sfget_out)
            if curr_n > 1:
                nlist[dim-1] = curr_n

    for dim in range(max_dims,0,-1):
        if nlist[dim-1] > 1:
            break

    return dim

################################################################################

def add_zeros(i, n):
    '''Generates string of zeros + str(i)'''

    ndigits_n = int(math.floor(math.log10(n-1)))
    if i == 0:
        nzeros = ndigits_n
    else:
        nzeros = ndigits_n - int(math.floor(math.log10(i)))

    return nzeros*'0'+str(i)

################################################################################

def chk_dir_rx(dir_nm):
    'Checks a directory exists and has +rx permissions'

    if not os.path.isdir(dir_nm):
        raise m8rex.NotAValidDir(dir_nm)

    if not os.access(dir_nm,os.X_OK):
        raise m8rex.NoXPermissions(dir_nm)

    if not os.access(dir_nm,os.R_OK):
        raise NoReadPermissions(dir_nm)

################################################################################

def chk_dir_wx(dir_nm):
    'Checks a directory exists and has +rx permissions'

    if not os.path.isdir(dir_nm):
        raise m8rex.NotAValidDir(dir_nm)

    if not os.access(dir_nm,os.X_OK):
        raise m8rex.NoXPermissions(dir_nm)

    if not os.access(dir_nm,os.W_OK):
        raise NoWritePermissions(dir_nm)

################################################################################

def chk_file_r(filename):
    'Checks if a file exists, is a regular file and is readable'

    if not os.path.exists(filename):
        raise m8rex.WrongPath(filename)

    if not os.path.isfile(filename):
        raise m8rex.NotAValidFile(filename)

    if not os.access(filename, os.R_OK):
        raise m8rex.NoReadPermissions(filename)

################################################################################

def isvalid(f,chk4nan=False):
    'Determines whether f is a valid RSF file'
    # Also returns msg with invalidity reason, or says if first x bytes are zero
    # Depends on sfin and sfattr output and formatting. To enhance robustness,
    # most parsing of sfin output messages is wrapped in try... except clauses

    msg = None
    is_rsf_ok = True

    sfin = 'sfin info='
    com_out = getout('sfin',['info=n', f])

    if com_out[:6] == 'sfin: ':
        is_rsf_ok = False

        try:
            err_msg = com_out.split(':')[2].strip()
        except:
            err_msg = ''

        if err_msg[:15] == 'No in= in file ':
            # Check if header is zero
            if os.path.getsize(f) == 0:
                msg = 'Empty header'
            else:
                msg = 'Missing in='

        elif err_msg[:22] == 'Cannot read data file ':
            # Check the binary
            bnr = err_msg[22:]
            (bnr_path, bnr_file) = os.path.split(bnr)
            if not os.access(bnr_path,os.X_OK):
                msg = 'Missing +x permissions to dir: ' + bnr_path
            elif not os.path.isfile(bnr):
                msg = 'Missing binary: ' + bnr
            elif not os.access(bnr,os.R_OK):
                msg = 'Missing read permissions to file: ' + bnr
        else: # Error not classified in this script
            msg = err_msg

    else: # No error message from sfin
        # Check for incomplete binaries
        com_out = getout('sfin',['info=y',f])
        line_list = com_out.split('\n')

        # Using a grep equivalent below because the "Actually..." error
        # message in sfin is not always on last line.

        # Clumsy implementation of grep below, with a touch of awk
        # Using Unix grep in original command will not work
        # because errors come via stderr and grep works on stdout

        err_list = []
        for line in line_list:
            if line[:6] == 'sfin: ':
                err_msg = line[6:]
                if err_msg != 'This data file is entirely zeros.':
                    if err_msg[:10] != 'The first ':
                        err_list.append(err_msg.strip())
        if err_list != []:
            is_rsf_ok = False
            # Showing only the first error. Probably the only one.
            msg = err_list[0]

    if is_rsf_ok: # no errors found by sfin
        if chk4nan:
            cmd_out = getout('sfattr', None, f)
            line_list = cmd_out.split('\n')
            for l in line_list:
                if '=' in l:
                    val = l.split('=')[1].strip().split('at')[0].strip()
                    if val in ('inf', '-inf', 'nan'):
                        is_rsf_ok = False
                        msg = 'Data contains ' + val

    return (is_rsf_ok, msg)

################################################################################

def list_invalid_rsf_files(dirname, flist, chk4nan=False):
    'Not recursive. Returns list of (file,msg) tuples.'

    chk_dir_rx(dirname)

    invalid_files = []

    rsf_files = [x for x in flist if os.path.splitext(x)[1]==ext]
    rsf_files.sort()

    for f in rsf_files:
        f = os.path.abspath(os.path.join(dirname,f))
        (is_rsf_ok, msg) = isvalid(f, chk4nan)
        if not is_rsf_ok:
            invalid_files.append((f, msg))

    return invalid_files

################################################################################

def list_valid_rsf_files(dirname, flist, chk4nan=False):
    'Not recursive'

    chk_dir_rx(dirname)

    valid_files = []

    rsf_files = [x for x in flist if os.path.splitext(x)[1]==ext]
    rsf_files.sort()

    for f in rsf_files:
        f = os.path.abspath(os.path.join(dirname,f))
        (is_rsf_ok, msg) = isvalid(f, chk4nan)
        if is_rsf_ok:
            valid_files.append(f)

    return valid_files

################################################################################

def msg(message, verb=True):
    'Print message using stderr stream'

    if type(verb) == int:
        if verb == 0:
            verb = False
        elif verb > 0:
            verb = True

    if verb:
        sys.stderr.write(str(message) + os.linesep)

################################################################################

def stdout_is_pipe():
    'True if writing to a pipe, false if writing to file (prog>file.rsf)'

    try:
        pos = sys.stdout.tell() # pos is offset in file
        return False # I am writing to a file
    except:
        return True

################################################################################

def trunc_or_append(n, ilist, append_val=None, verb=False):
    'Make a list be n elements regardless if it is longer or shorter'

    assert type(ilist) == list

    ndiff = n - len(ilist)

    if ndiff < 0:
        msg('Truncated %d elements' % ndiff, verb)
        return ilist[:n]
    elif ndiff > 0:
        msg('Appended %d elements' % ndiff, verb)
        return ilist + ndiff * [append_val]
    else: # n = len(ilist)
        return ilist

################################################################################

def gen_random_str(strlen):
    chars = string.letters + string.digits
    char_list = []
    for i in range(strlen):
        char_list += random.choice(chars)
    return ''.join(char_list)

################################################################################

def get_stdout_nm():
    'Find the name of the file that the stdout stream is writing to'
    # The kernel of this function initially provided by Jeff Godwin

    found_stdout = False

    for f in os.listdir('.'):
        # Comparing the unique file ID stored by the OS for the file stream
        # stdout with the known entries in the file table:
        if os.fstat(1)[1] == os.stat(f)[1]:
            found_stdout = True
            break

    if found_stdout:
        return f
    else:
        return None

################################################################################

def data_file_nm():
    'Get stdout file name, or if not possible, generate a random name'

    stdout = get_stdout_nm()

    if stdout == None:
        nm_len = 16 # Arbitrary random file name length
        return gen_random_str(nm_len)+'.rsf'
    else:
        return stdout

################################################################################

def chk_file_dims(filenm, ndims, verb=False):

    if int(getout('sfleftsize', 'i='+str(ndims), filenm, verb)) > 1:
        raise m8rex.NdimsMismatch(filenm, ndims)

################################################################################

def chk_par_in_list(par,avl):
    'Checks whether parameter is in acceptable values list'

    if par not in avl:
        if type(par) == str:
            raise m8rex.StringParamNotInAcceptableValueList(par, avl)
        # Fill in other types here

################################################################################

def chk_param_limit(parval, parnm, limval=0, comp='>'):
    '"Assert" uility for the CLI'
    if (comp == '>'  and parval <= limval) or \
       (comp == '>=' and parval <  limval) or \
       (comp == '<=' and parval >  limval) or \
       (comp == '<'  and parval >= limval):
        raise m8rex.ParBeyondLimit(parnm,limval,comp)

################################################################################

def valswitch(var,val1,val2):
    if var == val1:
        return val2
    else:
        return var

################################################################################

def switch(condition, val1, val2):
    if condition:
        return val1
    else:
        return val2
