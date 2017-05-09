#! /usr/bin/env python
'''Convert legacy SEPlib complex datasets to RSF

I.e. from 

esize=8 data_format=xdr_float

to 

esize=8 data_format=xdr_complex 

This combination is tolerated by SEPlib versions released after 2011-01-20,
and required by all versions of Madagascar. Previous to the date above, it
was impossible to have a complex single-precision dataset that was valid both
in SEPlib and Madagascar

This program opens the SEPlib file in read-write mode!

Handles in=stdin case (header and data in one file)
'''

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

import os, stat, time

import rsf.user.ivlad as ivlad
import rsf.user.m8rex as m8rex
import rsf.user.ooio  as ooio
import rsf.user.sf    as sf

###############################################################################

def main(par):

    f_nm = par.string('file') # Name of file to process
    preserve_t = par.bool('preserve_t', True) # Whether to preserve timestamp
    verb = par.bool('verb', False) # Say if file was converted or unchanged

    ivlad.chk_file_r(f_nm)
    
    if file_is_not_old_seplib_cmplx(f_nm):
        ivlad.msg(f_nm + ': Conversion not necessary', verb)
        return ivlad.unix_success

    if not os.access(f_nm, os.W_OK):
        # Not using an exception because this program may be called on long
        # lists of files
        ivlad.msg(f_nm + ': NEED WRITE PERMISSION', True)
        return ivlad.unix_error

    if preserve_t: 
        tstamp = get_timestamp(f_nm)

    bptr = sf.get(f_nm, 'in') # Pointer to binary

    if bptr == 'stdin': # Header and binary in a single file
        __replace_last_dataformat_statement(f_nm)
    else:
        f = open(f_nm, 'a')
        f.write('\n'+ooio.first_line_rsf_hdr_entry()) 
        f.write('\tdata_format="xdr_complex"\n')
        f.close()

    if preserve_t:
        os.utime(f_nm, tstamp)

    ivlad.msg(f_nm + ': Converted', verb)

    return ivlad.unix_success

################################################################################

def __replace_last_dataformat_statement(f_nm):
    'Turns the last xdr_float data_format into xdr_complex'

    # We will overwrite the last data_format="xdr_float" with
    # data_format=xdr_complex (without quotation marks, to keep the number
    # of characters the same, as to avoid overwriting any other info!

    xdrfloat_line_beg, xdrfloat_line = __get_xdrfloat_line_info(f_nm)
    
    f = open(f_nm, 'r+')
    f.seek(xdrfloat_line_beg)
    new_line = xdrfloat_line.replace('"xdr_float"','xdr_complex')
    f.write(new_line)
    f.close() 

################################################################################

def __get_xdrfloat_line_info(f_nm):

    f = open(f_nm, 'r')

    # Read the file line by line (not the entire file in memory),
    # to avoid reading the binary

    old_fpos = 0
    line = f.readline()
    while line:
        fpos = f.tell()
        if 'data_format="xdr_float"' in line:
            # There may be more than one such entry in the header
            # We should replace only the last one. So we first go through
            # header once, find which is the last one, then return to it
            xdrfloat_line = line
            xdrfloat_line_beg = old_fpos
        # Is the next "line" the beginning of the binary?
        teststr = f.read(4)
        if ooio.RSF_hdr_stream_end in teststr:
            break # Yes, the binary is "after the corner"
        else:
            f.seek(fpos) # Return to beginning of new line
        old_fpos = fpos
        line = f.readline()
    f.close()
    return xdrfloat_line_beg, xdrfloat_line 

################################################################################

def get_timestamp(f_nm):
    'Returns a tuple with access time, modification time'
    
    st = os.stat(f_nm)

    return (st[stat.ST_ATIME], st[stat.ST_MTIME]) 

################################################################################

def file_is_not_old_seplib_cmplx(f_nm):
    'Checks if a file has esize=8 and data_format=xdr_float'

    esize = sf.get(f_nm, 'esize')
    data_format = sf.get(f_nm, 'data_format')

    if esize == '8' and data_format == 'xdr_float':
        return False
    else:
        return True

################################################################################

if __name__ == '__main__':
    ivlad.run(main, ['file'])
