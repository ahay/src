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

def shell_script(target, source=None, env=None): 
    'Write the environment setup script'
    # Needs this specific interface because it will be called by a SCons Command
    # in RSFSRC/SConstruct

    shell = env['shell']
    user_shell = os.path.basename(os.environ.get('SHELL','/bin/sh'))

    rsfroot = env['RSFROOT']

    datapath = os.environ.get('DATAPATH','$RSFROOT/data')
    if datapath[-1] != '/':
        datapath += '/'

    pythonpath = os.environ.get('PYTHONPATH')
    mypath = get_local_site_pkgs(env['RSFROOT'])
    envpath = mypath.replace(rsfroot,'$RSFROOT',1)

    if pythonpath and not mypath in pythonpath.split(':'):
        pythonpath = ':'.join([envpath,pythonpath])
    else:
        pythonpath = envpath

    ldlibpath = os.environ.get('LD_LIBRARY_PATH')
    mypath = os.path.join(rsfroot,'lib')
    if ldlibpath and not mypath in ldlibpath.split(':'):
        ldlibpath = ':'.join(['$RSFROOT/lib',ldlibpath])
    else:
        ldlibpath = '$RSFROOT/lib'

    shrc = open(str(target[0]), 'w')
    shrc.write('#!/bin/%s\n' % shell)

    shenv = {
        'RSFROOT':[rsfroot,'Madagascar installation directory'],
        'PYTHONPATH':[pythonpath,'Making sure Python finds the rsf package'],
        'DATAPATH':[datapath,
                    'Default location for binary data files part of RSF datasets'],
        'LD_LIBRARY_PATH':[ldlibpath,
                           'Making sure Madagascar shared object files are found at runtime'],
        'MANPATH':['`manpath`:$RSFROOT/share/man',
                   'Making sure the "man" program finds Madagascar manual pages'],
        'PATH':['$RSFROOT/bin:$PATH',
                  'Making sure shell finds Madagascar executables']
        }

    keys = ('RSFROOT','PYTHONPATH','DATAPATH',
            'MANPATH','LD_LIBRARY_PATH','PATH')

    myrc = ''
    if shell == 'csh':
        for par in keys:
            if par == 'PATH':
                myrc += '\n# %s\nset path = ($RSFROOT/bin $path)' % \
                    shenv[par][1]
            else:
                myrc += '\n# %s\nsetenv %s %s\n' % \
                    (shenv[par][1],par,shenv[par][0])
    else:
        for par in keys:
            myrc += '\n# %s\nexport %s=%s\n' % \
                (shenv[par][1],par,shenv[par][0])

    if rsfroot != sysconfig.PREFIX and \
            user_shell in ({'sh':('bash','sh'),'csh':('tcsh','csh')}[shell]):
        sys.stderr.write('------------------------\n')
        sys.stderr.write('You may want to add the following to %s:\n' %
        {'sh':'.bashrc or .bash_profile','csh':'.cshrc'}[shell])
        sys.stderr.write(myrc)
        sys.stderr.write('------------------------\n')
    
    shrc.write(myrc)
    shrc.close()

    return 0

################################################################################

def get_pkgdir(root=None):
    'Return directory of the RSF Python package'
    
    return os.path.join(get_local_site_pkgs(root),'rsf')

################################################################################

if __name__ == '__main__':
    sys.exit(get_local_site_pkgs(verb=True))
