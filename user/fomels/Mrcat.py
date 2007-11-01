#!/usr/bin/env python
'Recursive sfcat (usefule for a large list of files)'

##   Copyright (C) 2007 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os,sys,string

sfcat = os.path.join(os.environ.get('RSFROOT'),'bin/sfcat')

def rcat(files,options,out):
    'Call sfcat recursively on a list of files'
    global sfcat
    if len(files) <= 2:
        stdout = os.dup(1)
        os.dup2(out,1)
        os.spawnv(os.P_WAIT,sfcat,['sfcat',]+options+files)
        os.dup2(stdout,1)
        os.close(out)
    else:
        middle = len(files)/2
        first = files[:middle]
        secon = files[middle:]
        ffile = string.join(first,'_')+'@'
        sfile = string.join(secon,'_')+'@'
        rcat(first,options,os.open(ffile,os.O_WRONLY | os.O_CREAT))
        if sfile != ffile:
            rcat(secon,options,os.open(sfile,os.O_WRONLY | os.O_CREAT))
        rcat([ffile,sfile],options,out)
        os.system('sfrm ' + ffile)
        if sfile != ffile:
            os.system('sfrm ' + sfile)

if __name__ == "__main__":
    options = []
    files = []
    for arg in sys.argv[1:]:
        if '=' in arg:
            options.append(arg)
        else:
            files.append(arg)
    rcat(files,options,1)

    sys.exit(0)


