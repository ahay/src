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
        return 0 # UNIX success
    else:
        return local_site_pkgs

################################################################################

def mk_sh_script(target, source=None, env=None): 
    'Write the (ba)sh environment setup script'
    # Needs this specific interface because it will be called by a SCons Command
    # in RSFSRC/SConstruct

    bsh = open(str(target[0]), 'w')
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
export RSFSRC=%s
''' % (env['RSFROOT'], get_local_site_pkgs(root=''), env['DATAPATH'], 
    env['RSFSRC']))

    bsh.close()

    return None

###############################################################################

def mk_csh_script(target, source=None, env=None):
    'Write the (t)csh environments setup script'
    
    csh = open(str(target[0]), 'w')
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
setenv RSFSRC %s
''' % (env['RSFROOT'], get_local_site_pkgs(root=''), env['DATAPATH'], 
    env['RSFSRC']))

    csh.close()

    return None

################################################################################

def get_pkgdir(root=None):
    'Return directory of the RSF Python package'
    
    return os.path.join(get_local_site_pkgs(root),'rsf')

################################################################################

if __name__ == '__main__':
    sys.exit(get_local_site_pkgs(verb=True))
