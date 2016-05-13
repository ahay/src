#!/usr/bin/env python
'Convert LAS-2 well logs to RSF'

##   Copyright (C) 2010 University of Texas at Austin
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

import os, sys
import numpy, m8r
from las import LASReader

def las2rsf(lasf,rsff):
    las = LASReader(lasf)
    rsf = m8r.Output(rsff)
    data = las.data2d
    shape = data.shape
    rsf.put('n1',shape[1])
    rsf.put('n2',shape[0])
    rsf.put('o2',las.start)
    rsf.put('d2',las.step)
    rsf.put('null',las.null)
    k = 0
    for name in las.curves.names:
        k += 1
        key = 'key%d' % k
        rsf.put(key,name)
        item = las.curves.items[name]
        rsf.put(key+'_units',item.units)
        desc = ' '.join(item.descr.translate(None,'"').split()[1:])
        rsf.put(name,desc)
    for name in las.well.names:
        item = las.well.items[name]
        desc = item.data.translate(None,'"')
        rsf.put(name,desc)
    rsf.write(data)
    rsf.close()

if __name__ == "__main__":
    usage = '''
Usage: %s file.las [file.rsf]

Check the output using sfheaderattr < file.rsf segy=n
    '''

    if len(sys.argv) < 2:
       print usage
       sys.exit(1)
    
    lasfile = sys.argv[1]

    if len(sys.argv) < 3:
        rsffile=os.path.splitext(lasfile)[0]+'.rsf'
    else:
        rsffile = sys.argv[2]

    las2rsf(lasfile,rsffile)
    sys.exit(0)

        
