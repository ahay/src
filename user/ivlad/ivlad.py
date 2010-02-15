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

import os, sys, commands, math, rsfprog, string, random

try: # Give precedence to local version
    import m8rex
except: # Use distributed version
    import rsfuser.m8rex as m8rex

# Operating system return codes
unix_success = 0
unix_error   = 1

# Other constants
pipe = ' | '
ext  = '.rsf'
fmt = {'int':'i', 'float':'f'} # For usage with struct.pack

# Horizontal ruler (for logs, screen messages):
hr = 80 * '-'

# Try to avoid old platform-dependent modules that will be deprecated.
# See http://docs.python.org/dev/lib/module-subprocess.html
try:
    import subprocess
    have_subprocess=True
except: # Python < 2.4
    have_subprocess=False
    from os import system
    from commands import getoutput

###############################################################################

def send_to_os(prog, arg=None, stdin=None, stdout=None, want=None, verb=False):
    '''Sends command to the operating system. Arguments:
    - prog. Executable to be run. STRING. The only non-optional argument.
    - arg. List of strings with program arguments. LIST (can be STRING for only 1 arg)
    - stdin. Filename STRING.
    - stdout. Filename STRING. Must not be specified if want=stdout
    - want: What to return. STRING: 'stdout' or 'stderr'. [Stderr not implemented yet]
    - verb: whether to print command before executing it. BOOL
    If stdout= is given, then the function writes to that file and returns None
    If want=stdout, then a STRING with stripped newlines is returned.'''

    if want:
        if want not in ('stdout','stderr'):
            sys.stderr.write('The "want" argument to send_to_os must be "stdout" or "stderr"')
            sys.exit(error)
        if stdout and want == 'stdout':
            sys.stderr.write('No stdout= should be given to send_to_os when want="stdout"')
            sys.exit(error)

    # Build the [prog, args] list
    if arg:
        if str(arg) == arg: # it's a string
            arg = [arg]     # make it a list
        arg.insert(0,prog)
        cmdlist = arg
    else:
        cmdlist = [prog]

    # Build command string for printing or Python < 2.4
    if verb or not have_subprocess:
        cmd4print = cmdlist[:]
        if stdin:
            cmd4print.append('<')
            cmd4print.append(stdin)
        if stdout:
            cmd4print.append('>')
            cmd4print.append(stdout)
        command = ' '.join(cmd4print)
        if verb:
            print command

    if have_subprocess:
        if stdin:
            finp = open(stdin,'r')
        else:
            finp = None
        if stdout:
            fout = open(stdout,'w')
        else:
            fout = None
        if want == 'stdout':
            s = subprocess.Popen(cmdlist,stdin=finp,stdout=subprocess.PIPE)
            output = s.communicate()[0]
            return output.rstrip('\n') # Remove the newline character
        else:
            subprocess.Popen(cmdlist,stdin=finp,stdout=fout).wait()
            return None
    else: # no subprocess module present
        if want == 'stdout':
            s = getoutput( command )
            return s
        elif stdout:
            system(command)
            return None

###############################################################################

def readaxis( inp, axisnr, verb=False ):
    '''Reads n,o,d for one axis of a Madagascar hypercube
    - inp: filename. STRING.
    - axisnr: INT
    - verb: BOOL'''

    ax = str(axisnr)

    n = send_to_os('sfget',
                   arg = ['parform=n','n'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    o = send_to_os('sfget',
                   arg = ['parform=n','o'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    d = send_to_os('sfget',
                   arg = ['parform=n','d'+ax],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
    return( int(n), float(o), float(d))

###############################################################################

def ndims(filename):
    'Returns nr of dims in a file, ignoring trailing length 1 dimensions'

    max_dims = 9
    nlist = [1] * max_dims

    for dim in range(1,max_dims+1):
        command = '<'+filename+' sfget parform=n n' + str(dim)
        sfget_out = commands.getoutput(command)
        if sfget_out[:15] != 'sfget: No key n':
            curr_n = int(sfget_out)
            if curr_n > 1:
                nlist[dim-1] = curr_n

    for dim in range(max_dims,0,-1):
        if nlist[dim-1] > 1:
            break

    return dim

####################################################################

def add_zeros(i, n):
    '''Generates string of zeros + str(i)'''

    ndigits_n = int(math.floor(math.log10(n-1)))
    if i == 0:
        nzeros = ndigits_n
    else:
        nzeros = ndigits_n - int(math.floor(math.log10(i)))

    return nzeros*'0'+str(i)

####################################################################

def execute(command, verb=False):
    '''Echoes a command to screen, then executes it'''

    if verb:
        print command
    os.system(command)

################################################################################

def chk_dir(mydir):
    'Checks a directory exists and has +rx permissions'

    if mydir == None:
        rsfprog.selfdoc()
        return unix_error

    if not os.path.isdir(mydir):
        print mydir + ' is not a valid directory'
        return unix_error

    if not os.access(mydir,os.X_OK):
        print mydir + ' lacks +x permissions for ' + os.getlogin()
        return unix_error

    if not os.access(mydir,os.R_OK):
        print mydir + ' lacks read permissions for ' + os.getlogin()
        return unix_error

    return unix_success

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

def isvalid(f,chk4nan=False):
    'Determines whether f is a valid RSF file'
    # Also returns msg with invalidity reason, or says if first x bytes are zero
    # Depends on sfin and sfattr output and formatting. To enhance robustness,
    # most parsing of sfin output messages is wrapped in try... except clauses

    msg = None
    is_rsf_ok = True

    sfin = 'sfin info='
    com_out = commands.getoutput(sfin+'n '+f)
    
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
        com_out = commands.getoutput(sfin+'y '+f)
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
            cmd_out = commands.getoutput('sfattr <'+f)
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
    
    if chk_dir(dirname) == unix_error:
        return None # I should really use exceptions here. No time now.

    invalid_files = []

    rsf_files = filter(lambda x:os.path.splitext(x)[1]==ext, flist)
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
    
    if chk_dir(dirname) == unix_error:
        return None # I should really use exceptions here. No time now.

    valid_files = []

    rsf_files = filter(lambda x:os.path.splitext(x)[1]==ext, flist)
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

    if verb:
        sys.stderr.write(str(message)+'\n')

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
    return ''.join(random.choice(chars) for i in xrange(strlen))

