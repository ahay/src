#!/usr/bin/env python
'''
Generates a Debian/Ubuntu package.
'''
# Copyright (C) 2009 The University of Texas at Austin
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

import os, shutil

def makedir(path,name):
    dir = os.path.join(path,name)
    if not os.path.isdir(dir):
        os.mkdir(dir)
        return dir

PACKVERSION = '1' # package version

DEBNAME = 'madagascar'
DEBVERSION = '0.9.8'
DEBMAINT = 'Sergey Fomel [sergey.fomel@gmail.com]'
DEBARCH = os.uname()[4]

DEBDEPENDS = '' # what are we dependent on?
DEBDESC = 'Geophysical data processing and reproducible numerical experiments'

# This is the debian package we're going to create

name = '-'.join([DEBNAME, DEBVERSION, PACKVERSION, DEBARCH])
dir = makedir(os.environ.get('DATAPATH','.'),name)

# This copies the necessary files into place into place.

root = os.path.join(dir,'usr/local/RSFROOT')
rsfroot = os.environ.get('RSFROOT')

for subdir in ('bin','doc','include','lib','man'):
    print 'Packing %s...' % subdir
    shutil.copytree(os.path.join(rsfroot,subdir),
                    os.path.join(root,subdir))

# Now to create the control file:

CONTROL_TEMPLATE = '''
Package: %s
Priority: extra
Section: misc
Installed-Size: %s
Maintainer: %s
Architecture: %s
Version: %s-%s
Depends: %s
Description: %s

'''

installed_size = os.popen('du -bs %s' % root).read().split()[0]
control_info = CONTROL_TEMPLATE % (
    DEBNAME, installed_size, DEBMAINT, DEBARCH, DEBVERSION,
    PACKVERSION, DEBDEPENDS, DEBDESC)

debdir = makedir(dir,'DEBIAN')
DEBCONTROLFILE = os.path.join(debdir, 'control')
f = open(DEBCONTROLFILE, 'w')
f.write(control_info)
f.close()

command = 'fakeroot dpkg-deb -b ' + dir
print command
os.system(command)

