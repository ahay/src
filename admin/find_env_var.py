#!/usr/bin/env python
'''
Finds new environment variables used in Python and C programs
Looks for "environ.get" in py progs (including SConstruct files), 
and "getenv" in c progs
Run it from RSFSRC, i.e. admin/find_env_var.py
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

import rsf.user.ivlad as ivlad
import os, sys, string

if not hasattr(' ', 'partition'):
    ivlad.msg('Please get Python 2.5 or greater')
    sys.exit(ivlad.unix_error)

try:
    import subprocess
except:
    print 'Please get Python 2.4 or greater, or just the subprocess module'
    sys.exit(unix_error)

###############################################################################

def grep_from_py(filename, string):
    'Runs "grep string filename", and returns the result'
    proc = subprocess.Popen('grep ' + string + ' ' + filename,
                            shell=True,
                            stdout=subprocess.PIPE)
    stdout_val = proc.communicate()[0]
    return stdout_val

###############################################################################

def extract_env_var_nm(line, sep):

    out = []
    s1 = line.strip()
    s2 = s1.partition(sep)[2]
    s3 = s2.strip().lstrip('(')
    s4 = s3.lstrip('\'')
    if s3 == s4:
        # No single quote was present
        quotsign = '"'
        s5 = s4.lstrip(quotsign)
    else:
        quotsign = "'"
        s5 = s4
    env_var, x, rest = s5.partition(quotsign)
    # Assumption that a valid env variable name must have at least a letter
    valid_env_var_name = False
    for letter in string.ascii_letters:
        if letter in env_var:
            valid_env_var_name = True
            break
    if valid_env_var_name and env_var != "key+'PATH')" and env_var != 'env)':
        out.append(env_var)
    if sep in rest:
        out += extract_env_var_nm(rest, sep)
    return out

###############################################################################

known_vars = '''
RSF_DATASERVER RSFBOOK RSFFIGS RSFALTFIGS RSFMEMSIZE RSFSRC TMPDATAPATH
LATEX2HTML DEFAULT_PAPER_SIZE FATMULT GIFBORDER GIFDELAY IMAGE_TYPE PATTERNMULT
PLOTSTYLE PPI PPMSCALE PSBORDER PSPRINTER PSTEXPENOPTS VPLOTFONTDIR 
VPLOTSPOOLDIR WSTYPE CWPROOT DISPLAY HOME LD_LIBRARY_PATH MATLABPATH XAUTHORITY
PYTHONPATH RSFROOT DATAPATH PATH
'''.split()

py_files_with_no_ext = '''
SConstruct latex2wiki pscons sfdoc sfkill sftop sftour
'''.split()

###############################################################################

def main():

    dirs_to_check = 'api book framework pens plot site_scons su system user'
    func_names = {'c':'getenv', 'py':'environ.get'}
    env_var_list = []

    for d in dirs_to_check.split():
        for root,dirs,files in os.walk(d):
            # Avoid Subversion bookkeeping directories
            if '.svn' not in root:
                for file in files:
                    shortnm, ext = os.path.splitext(file)
                    fullnm = os.path.join(root, file)
                    if ext == '.c':
                        filetype = 'c'
                    elif ext == '.py' or shortnm in py_files_with_no_ext:
                        filetype = 'py'
                    else:
                        continue
                    grep4env = grep_from_py(fullnm, func_names[filetype])
                    grep4env = grep4env.strip()
                    if grep4env == '':
                        continue
                    else:
                        sys.stderr.write('.')
                        for line in grep4env.split('\n'):
                            vl = extract_env_var_nm(line,func_names[filetype])
                            if vl != []:
                                for var in vl:
                                    if var not in env_var_list:
                                        env_var_list.append(var)

    sys.stderr.write('\n')
    if env_var_list != []:
        env_var_list.sort()
        for v in env_var_list:
            if v not in known_vars:
                print v

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
