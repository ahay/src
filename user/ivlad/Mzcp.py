#! /usr/bin/env python
''''Copies header of float file, fills output binary with zeros
Usage: sfzcp file1.rsf file2.rsf'''

# Copyright (C) 2009 Ioan Vlad
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

import commands, os, sys

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

###############################################################################

def main(argv=sys.argv):

    if len(sys.argv) != 3:
        sys.stderr.write('Usage: sfzcp file1.rsf file2.rsf\n')
        return ivlad.unix_error

    f1 = sys.argv[1] + ' '
    f2 = sys.argv[2]

    f1_type = commands.getoutput('<' + f1 + 'sfgettype')
    if f1_type != 'SF_FLOAT':
        sys.stderr.write('Other types not yet implemented\n')
        return ivlad.unix_error

    dims_str = commands.getoutput('<' + f1 + 'sffiledims parform=n')
    ndims = int(dims_str.split(':')[0])

    cmd = 'sfspike mag=0 '
    fields = 'n o d unit label'.split()

    for i in range(ndims):
        for z in fields:
            cmd_local = '<' + f1 + 'sfget parform=y ' + z+str(i+1)
            cmd += commands.getoutput(cmd_local) + ' '

    cmd += ' >' + f2

    os.system(cmd)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
