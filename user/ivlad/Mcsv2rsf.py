#!/usr/bin/env python
'''Convert a delimited-text ASCII file to RSF binary floating point or int.
Zeros will be added if number of elements is not the same in each row.
n1 and n2 are computed automatically. For consistency with sfdisfil and
sfmatmult, output is C-style order (row-first), i.e. rows in input file
become dimension-1 columns in output. Output encoding is native. If n2=1 in
output, the second dimension will not be written to the header.'''

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

from __future__ import print_function
import csv, struct, sys, os

try: # Give precedence to local version
    import ivlad, ooio
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.ooio  as ooio

################################################################################

def main(par):

    # Input parameters

    delim = par.string('delimiter',',') # Separator between values in input file
    numtype = par.string('dtype', 'float') # Input type

    ivlad.chk_par_in_list(numtype, ooio.dtype_avl)

    # Process parameters

    verb = par.bool('verb', False) # Whether to echo n1, n2, infill/truncation
    debug = par.bool('debug', False) # Extra verbosity for debugging
    truncate = par.bool('trunc', False)
    # Truncate or add zeros if nr elems in rows differs

    header = par.bool('header',False) # If the first line is a header

    # Output parameters

    o = [
    par.float('o1', 0.), # Origin of axis 1 in output (rows in input)
    par.float('o2', 0.)] # Origin of axis 2 in output (columns in input)

    d = [
    par.float('d1', 1.), # Axis 1 sampling
    par.float('d2', 1.)] # Axis 2 sampling

    unit = [
    par.string('unit1', 'unknown'),
    par.string('unit2', 'unknown')]

    lbl = [
    par.string('label1', 'unknown'),
    par.string('label2', 'unknown')]

    ##### End reading parameters #####

    stdin = csv.reader(sys.stdin, delimiter=delim)

    # Copy stdin so we can go back through it
    lines = []
    i = 0
    nr_cols_cst = True # Whether the nr of values in rows is constant

    if header: # the first line contains header keys
        line = next(stdin)
        k = 0
        for name in line:
            k += 1
            print('key%d=%s' % (k,name))

    # Find max nr of elements in a row
    for line in stdin:
        if line == []: # throw away blank lines
            continue
        curline = [float(x) for x in [x or '0' for x in line]]
        if numtype == 'int':
            curline = [int(x) for x in curline]
        lines.append(curline)
        llen = len(curline)
        i+=1 # We have successfully read line i
        if i==1:
            n2 = llen
        elif llen != n2:
            nr_cols_cst = False
            if (llen < n2 and truncate) or (llen > n2 and not truncate):
                n2 = llen

    n1 = len(lines)

    if not nr_cols_cst: # Truncate or append as necessary
        for i in range(n1):
            line = lines[i]
            lines[i] = ivlad.trunc_or_append(n2, line, 0, verb)

    # Avoiding to add a second dimension of length 1
    if n1 == 1:
        ndim_out = 1
        n = [n2]
        o = [o[0]]
        d = [d[0]]
        unit = [unit[0]]
        lbl = [lbl[0]]
    else:
        ndim_out = 2
        n = [n2, n1]

    out = ooio.RSFfile(ooio.stdout,par,ndim=ndim_out,intent='out',dtype=numtype)
    out.set_hdr_info(n, o, d, unit, lbl)

    if debug:
        out.print_self('out')
        out.hdr.print_self('out.hdr')
        out.dat.print_self('out.dat')
        ivlad.msg(ivlad.hr)

    for line in lines:
        for val in line:
            out.write(val)

    return ivlad.unix_success

##############################################

if __name__ == '__main__':
    if len(sys.argv)==1: # fix for python3
        sys.argv.append('-')
    ivlad.run(main)
