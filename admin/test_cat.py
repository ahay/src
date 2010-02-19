#! /usr/bin/env python
'''
Tests how many files can sfcat concatenate
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

import subprocess

nfiles = 511 # May fail for n > 510 on some systems
cmd_root = 'sfspike n1=10 > '
fname_root = 'test_junk_'
all_files = ''

for i in range(nfiles):
    fname = fname_root + str(i) + '.rsf'
    cmd = cmd_root + fname
    print cmd
    subprocess.call(cmd, shell=True)
    all_files += ' ' + fname

print '\n# Created all files\n'

cmd = 'sfcat axis=2' + all_files + ' > ' + fname_root + 'out.rsf'
print cmd

try:
    subprocess.call(cmd, shell=True)
except:
    print 'sfcat failed'
    subprocess.call('sfrm'+ all_files)

