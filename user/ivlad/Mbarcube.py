#!/usr/bin/env python
'''<cube.rsf sfbarcube color=[i] verb=[n] dryrun=[n] >file.vpl
Like sfgrey3, but computes the colorbar using sfbar'''

# Copyright (C) 2011 Ioan Vlad
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

import random, subprocess, sys

try: # Give precedence to local version
    import ivlad, sf
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.sf as sf

################################################################################
# <inp sfbar >bar.rsf

#<i11.transp.w.rsf sfbyte gainpanel=a pclip=98 | sfgrey3 frame1=105 frame2=235 frame3=2 scalebar=y color=j bar=bar_noarg.rsf >junk.vpl

def main(par):

    verb   = par.bool('verb'  , False)
    dryrun = par.bool('dryrun', False)
    
    # sfbyte arguments:
    bytepars = {}
    for key in ['gainpanel', 'pclip', 'clip']:
        bytepars[key] = par.string(key, None)

    # sfgrey3 arguments:
    color  = par.string('color', None)
    
    lbl = int(100000000*random.random())
    suffix = '.' + str(lbl) + ivlad.ext
    prefix = sys.argv[0].rstrip('.py') + '.'
    
    tmp_filenm = prefix + 'inp' + suffix
    bar_filenm = prefix + 'bar' + suffix
    
    # Redirect the input
    if not dryrun:
        tmp_file = open(tmp_filenm, 'w')
        tmp_file.write(sys.stdin.read())
        tmp_file.close()
       
    sf._set_xmode('l')
       
    cmd_list = ['# Redirect input to ' + tmp_filenm,
                sf.bar(tmp_filenm, bar_filenm)]
    ivlad.exe(cmd_list, verb, dryrun)
    
    if dryrun:
        cmd = '<' + tmp_filenm + ' sfbyte | sfgrey3 scalebar=y bar='
        cmd+= bar_filenm + ' >stdout'
        ivlad.msg(cmd, verb)
    else:  
        p1 = subprocess.Popen(['sfbyte'],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE)
        tmp_file = open(tmp_filenm, 'r')
        p1.stdin.write(tmp_file.read())
        p2 = subprocess.Popen(['sfgrey3', 'scalebar=y', 'bar='+bar_filenm],
                              stdin=p1.stdout,
                              stdout=subprocess.PIPE)
        sys.stdout.write(p2.stdout.read())
        p1.stdin.close()
        tmp_file.close()
        p2.stdout.close()
    
    cmd = sf.rm([tmp_filenm, bar_filenm])
    
    ivlad.exe(cmd, verb, dryrun)
    
    return ivlad.unix_success

################################################################################

if __name__ == '__main__':
    ivlad.run(main)
