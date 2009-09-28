#!/usr/bin/env python
"""
Splits file into slices along the last dimension

Usage:

sfsplit inp=file.rsf outdir=[file_split.rsf] nthick=[1]

Parameter nthick gives the maximum thickness of a slice. The last slice may
be thinner."""

import os, sys, rsfprog, commands, math
from rsfuser.ivlad import unix_error, unix_success, ext, ndims, add_zeros, execute

try:
    import rsf
except:
    import rsfbak as rsf

################################################################################

def main(argv=sys.argv):

    par = rsf.Par(argv)

    verb = par.bool('verb', False)

    inp = par.string('inp')
    if inp == None:
        rsfprog.selfdoc()
        return unix_error

    outdir = par.string('outdir')
    if outdir == None:
        out_basenm = os.path.basename(inp).rstrip(ext) + '_split'
        outdir = out_basenm + ext

    nthick = par.int('nthick', 1)
    if nthick < 1:
        print 'Slice thickness must be 1 or higher'
        return unix_error

    lastdim = str(ndims(inp))

    n = int(commands.getoutput('sfget <' + inp + ' parform=n n' + lastdim))
    if nthick >= n:
        print 'Slice thickness (%d) should be < axis length (%d)' % (nthick, n)
        return unix_error

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Slices will be called 01.rsf, etc. Avoid datapath aliasing:
    dpath = 'datapath=$DATAPATH/' + out_basenm + '_'

    nslices_whole = int(math.floor(n/nthick))

    f = 0

    for i in range(nslices_whole):

        i_str = add_zeros(i,nslices_whole)
        i_slc = os.path.join(outdir, i_str + '_stdout' + ext)

        command = '<%s sfwindow f%s=%d ' % (inp, lastdim,f)
        command += 'n%s=%d %s > %s' % (lastdim, nthick, dpath, i_slc)
        execute(command, verb)

        f += nthick

    ist = n - nslices_whole * nthick # Incomplete Slice Thickness

    if ist > 0:
        i_str = qsub.add_zeros(i+1,nslices_whole)
        i_slc = os.path.join(outdir, i_str + '_stdout' + ext)

        command = '<%s sfwindow f%s=%d ' % (inp, lastdim, f)
        command += 'n%s=%d > %s' % (lastdim, ist, i_slc)
        execute(command, verb)

    return unix_success

##############################################

if __name__ == '__main__':

    try:
        status = main()
    except:
        status = unix_error

    sys.exit(status)
