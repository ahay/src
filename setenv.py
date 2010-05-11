#!/usr/bin/env python

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

import os, sys

if sys.version_info[:2] < (2, 7):
    import distutils.sysconfig as sysconfig
else:
    import sysconfig

dpath = '/var/tmp'
root = os.environ.get('RSFROOT','.')
unix_success = 0

################################################################################

def get_local_site_pkgs(root=None, verb=False):
    'Get the directory that should be added to PYTHONPATH'

    central_site_pkgs = sysconfig.get_python_lib()

    # This is actually platform dependent (Mac vs. Linux), because Python is not
    # as cross-platform as it should be. When we need to install in a 
    # canonical location on Mac, platform detection will be included
    prefix = sysconfig.PREFIX
    
    if root == None:
        root = os.environ.get('RSFROOT',os.environ['HOME'])

    if central_site_pkgs[:len(prefix)] == prefix:
        local_site_pkgs = central_site_pkgs.replace(prefix,root,1)
    else:
        local_site_pkgs = os.path.join(root,'lib')

    if verb:
        print local_site_pkgs
        return unix_success
    else:
        return local_site_pkgs

################################################################################

def init_globs(rsfroot=None, user_dpath=None):
    'Initialize datapath to make for an easy CLI interface'
    global dpath, root

    if user_dpath != None:
        dpath = user_dpath

    if rsfroot != None:
        root=rsfroot

################################################################################

def mk_sh_script(script_nm, rsfroot):
    'Write the (ba)sh environment setup script'

    global dpath, root

    bsh = open(script_nm, 'w')
    bsh.write('''#!/bin/sh

export RSFROOT=%s
if [ -n "$PYTHONPATH" ]; then
export PYTHONPATH=${PYTHONPATH}:$RSFROOT%s
else
export PYTHONPATH=$RSFROOT/lib
fi
export PATH=$RSFROOT/bin:$PATH
export DATAPATH=%s
export MANPATH=$RSFROOT/share/man:$(manpath)
export LD_LIBRARY_PATH=$RSFROOT/lib:$LD_LIBRARY_PATH
''' % (rsfroot, get_local_site_pkgs(root=''), dpath))
    bsh.close()

    return unix_success

###############################################################################

def mk_csh_script(script_nm, rsfroot):
    'Write the (t)csh environments setup script'
    global dpath, root

    csh = open(script_nm, 'w')
    csh.write('''#!/bin/csh

setenv RSFROOT %s
if ($?PYTHONPATH) then
setenv PYTHONPATH ${PYTHONPATH}:$RSFROOT%s
else
setenv PYTHONPATH $RSFROOT/lib
endif
set path = ($RSFROOT/bin $path)
setenv DATAPATH %s
setenv MANPATH $RSFROOT/share/man:`manpath`
setenv LD_LIBRARY_PATH $RSFROOT/lib:$LD_LIBRARY_PATH
''' % (rsfroot, get_local_site_pkgs(root=''), dpath))
    csh.close()

    return unix_success

################################################################################

def get_pkgdir(root=None):
    'Return directory of the RSF Python package'
    
    return os.path.join(get_local_site_pkgs(root),'rsf')

################################################################################

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs == 1:
        sys.exit(get_local_site_pkgs(verb=True))
    elif nargs == 4:
        if sys.argv[1] == 'sh':
            sys.exit(mk_sh_script(sys.argv[2],sys.argv[3]))
        elif sys.argv[1] == 'csh':
            sys.exit(mk_csh_script(sys.argv[2],sys.argv[3]))
        else:
            sys.stderr.write('Second argument to setenv should be sh or csh\n')
            sys.exit(1)
    else:
        sys.stderr.write('Usage: setenv.py (c)sh script_name <RSFROOT>\n')
        sys.exit(1)
