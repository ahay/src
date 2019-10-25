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

from __future__ import print_function, division, absolute_import
import os, sys, distutils.sysconfig

################################################################################

def get_local_site_pkgs(root=None, verb=False):
    'Get the directory that should be added to PYTHONPATH'

    central_site_pkgs = distutils.sysconfig.get_python_lib()

    # This is actually platform dependent (Mac vs. Linux), because Python is not
    # as cross-platform as it should be. When we need to install in a
    # canonical location on Mac, platform detection will be included
    prefix = distutils.sysconfig.PREFIX

    if root == None:
        root = os.environ.get('RSFROOT', prefix)

    if central_site_pkgs[:len(prefix)] == prefix:
        local_site_pkgs = central_site_pkgs.replace(prefix,root,1)
    else:
        local_site_pkgs = os.path.join(root,'lib')

    if verb:
        print(local_site_pkgs)
        return 0 # UNIX success
    else:
        return local_site_pkgs

################################################################################

def shell_script(target, source=None, env=None):
    'Write the environment setup script'
    # Needs this specific interface because it will be called by a SCons Command
    # in RSFSRC/SConstruct

    import configure
    import imp
    config = imp.load_source('config', 'config.py')

    shell = env['shell']
    rsfroot = config.RSFROOT

    pypath = get_local_site_pkgs(rsfroot)
    datapath = os.environ.get('DATAPATH','/var/tmp')
    if datapath[-1] != os.sep:  datapath += os.sep

    ldlibpath = os.path.join(rsfroot,'lib')

    manpath = '`manpath`:$RSFROOT/share/man'

    shrc = open(str(target[0]), 'w')
    shrc.write('#!/bin/%s\n' % shell)

    shenv = {
        'RSFROOT':(rsfroot,'Madagascar installation directory'),
        'PYTHONPATH':(pypath,'Python packages'),
        'JULIA_LOAD_PATH':(ldlibpath,'Julia API'),
        'DATAPATH':(datapath,'binary data files part of RSF datasets'),
        'LD_LIBRARY_PATH':(ldlibpath,'shared object files'),
        'MANPATH':(manpath,'manual pages'),
        'RSFSRC':(os.getcwd(),'Madagascar source directory'),
        'PATH':('$RSFROOT/bin:$PATH','executables')
        }

    keys = ('RSFROOT','RSFSRC','PYTHONPATH','JULIA_LOAD_PATH', 'DATAPATH',
            'MANPATH','LD_LIBRARY_PATH','PATH')

    myrc = ''
    for par in keys:
        comment = shenv[par][1]
        value = shenv[par][0]
        if par != 'RSFROOT':
            value = value.replace(rsfroot+os.sep,'$RSFROOT'+os.sep)
        redefine =  (par == 'PYTHONPATH' or par == 'JULIA_LOAD_PATH' or par == 'LD_LIBRARY_PATH')

        myrc += '\n# Path for %s\n' % comment
        if shell == 'csh':
            if par == 'PATH':
                myrc += 'set path = ($RSFROOT/bin $path)\n'
            else:
                if redefine:
                    myrc += 'if ($?%s) then\n' % par
                    if par.startswith('JULIA'):
                        myrc += 'setenv %s %s:${%s}:\n' % (par,value,par)
                    else:
                        myrc += 'setenv %s %s:${%s}\n' % (par,value,par)
                    myrc += 'else\n'
                myrc += 'setenv %s %s:\n' % (par,value)
                if redefine:
                    myrc += 'endif\n'
        else:
            if redefine:
                myrc += 'if [ -n "$%s" ]; then\n' % par
                if par.startswith('JULIA'):
                    myrc += 'export %s=%s:${%s}:\n' % (par,value,par)
                else:
                    myrc += 'export %s=%s:${%s}\n' % (par,value,par)
                myrc += 'else\n'
            if par == 'MANPATH':
                myrc += 'unset MANPATH\n' # necessary for recent bash shells
            myrc += 'export %s=%s\n' % (par,value)
            if redefine:
                myrc += 'fi\n'

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
