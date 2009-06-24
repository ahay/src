#! /usr/bin/env python
'''Resamples a 2-D dataset to the desired picture resolution, with antialias
Only one of the h and w parameters needs to be specified.
If prar=n, no action will be taken on axis for which h/w was not specified
If prar=y and only one par (h or w) is specified, the picture will scale
along both axes until it is of the specified dimension.'''

# Copyright (C) 2007, 2009 Ioan Vlad
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

import os, sys, rsfprog

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

###############################################################################

def main(argv=sys.argv):

    par = rsf.Par(argv)

    inp = par.string('inp') # input file
    out = par.string('out') # output file
    if None in (inp, out):
        rsfprog.selfdoc()
        return ivlad.unix_error

    verb = par.bool('verb', False) # if y, print system commands, outputs

    leftsize = int(ivlad.send_to_os('sfleftsize', arg='i=2', stdin=inp,
                              want='stdout', verb=verb))

    if leftsize > 1:
        sys.stderr.write('Input must be 2-D file\n')
        return ivlad.unix_error
    del leftsize

    n1 = ivlad.send_to_os('sfget',
                    arg = ['parform=n','n1'],
                    stdin = inp,
                    want = 'stdout',
                    verb = verb)
    n2 = ivlad.send_to_os('sfget',
                    arg = ['parform=n','n2'],
                    stdin = inp,
                    want = 'stdout',
                    verb = verb)
    n1 = int(n1)
    n2 = int(n2)

    # Read and check the h and w parameters
    # The debacle below could be avoided if rsf module had a function specifying
    # whether a parameter for whom a default value is required was given at the
    # command line or not.
    impossible=-100000
    none=impossible #even uglier hack, for self-doc to show user-friendly message
    h = par.int('h', none) # output height
    w = par.int('w', none) # output width

    if w==impossible:
        if h==impossible:
            sys.stderr.write('Either w or h must be specified!!\n')
            return ivlad.unix_error
        else: # h was read from the command line
            w = None
    else: # w was read from the command line
        if h==impossible:
            h = None
        else: # both w and h were read. Check for sanity:
            if w <= 0:
                sys.stderr.write('w is a width, must be >0\n')
                return ivlad.unix_error
            if h <= 0:
                sys.stderr.write('h is a height, must be >0\n')
                return ivlad.unix_error
            if (h,w)==(n1,n2):
                sys.stderr.write('Change h or w if you want out!=inp\n')
                ivlad.send_to_os('sfcp', arg=[inp, out], verb=verb)
                return ivlad.unix_success
            if h==n1:
                h = None # No interp needed on axis 1
            if w==n2:
                w = None # Same for axis 2
    del impossible
    # Now h and w are either None, or have a value worth doing interpolation on.

    # Transform h and w to pixels, if they are not
    # No default value for par.string -- Quirk of rsf, replicated in rsfbak
    unit = par.string('unit') # unit of h and w. Can be: px(default), mm, cm, in
    if unit: # something was read from command line
        if unit not in ('mm','cm','in','px'):
            sys.stderr.write('unit must be: mm/cm/in/px\n')
            return ivlad.unix_error
        if unit != 'px':
            ppi = par.int('ppi') # output resolution (px/in). Necessary when unit!=px
            if ppi <= 0:
                sys.stderr.write('ppi must be > 0\n')
                return ivlad.unix_error
            # Transform w and h to px
            if unit == 'in':
                scale = 1
            elif unit == 'mm':
                scale = 254
            elif unit == 'cm':
                scale = 25.4
            if w:
                w *= ppi / float(scale)
            if h:
                h *= ppi / float(scale)
            # Don't worry, we'll convert to int after prar block
            del scale
    else: # unit not in prog arguments 
        unit = 'px'
    # Now h and w are either None, or in pixels

    # If prar=y, then h and/or w define a bounding box.
    # Find the dimensions of the image inside this box
    prar = par.bool('prar', True) # if y, PReserve Aspect Ratio of input
    if prar: # preserve aspect ratio
        if h and not w:
            w = n2 * float(h) / n1
        elif w and not h:
            h = n1 * float(w) / n2
        else: # Full bounding box specified
            hscale = float(h) / n1
            wscale = float(w) / n2
            if hscale < wscale:
                w = n2 * hscale
            else:
                h = n1 * wscale

    h = int(h)
    w = int(w)
    if h == n1:
        h = None # No interp needed on axis 1
    if w == n2:
        w = None # Same for axis 2
    if h < 2 or w < 2:
        sys.stderr.write('sanity check failed\n')
        return ivlad.unix_error

    # Put tmp files together with the binaries,
    # so that if prep4plot crashes, user is not
    # left with junk files all over his dir
    tmp = os.path.join(os.environ['DATAPATH'],
                       os.path.split(inp)[1],
                       '.prep4plot_junk_')

    # Interpolation and, if needed, bandpass 
    if h:
        d1 = ivlad.send_to_os('sfget',
                   arg = ['parform=n','d1'],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
        d1 = float(d1) * (n1-1)/float(h-1)
        if h < n1:
            ready_for_remap_1 = tmp + '1'
            ivlad.send_to_os('sfbandpass',
                       arg='fhi='+str(0.5/d1),
                       stdin=inp,
                       stdout=ready_for_remap_1,
                       verb=verb)
            rem2del_junk1 = True
        else:
            ready_for_remap_1 = inp
            rem2del_junk1 = False
        if w:
            out_remap1 = tmp + '2'
            rem2del_junk2 = True
        else:
            out_remap1 = out
            rem2del_junk2 = False
        ivlad.send_to_os('sfremap1',
                   arg=['n1='+str(h), 'd1='+str(d1)],
                   stdin=ready_for_remap_1,
                   stdout=out_remap1,
                   verb=verb)
        if rem2del_junk1:
            ivlad.send_to_os('sfrm',
                       arg=ready_for_remap_1,
                       verb=verb)
    else: # no action on axis 1
        out_remap1 = inp
        rem2del_junk2 = False

    if w:
        d2 = ivlad.send_to_os('sfget',
                   arg = ['parform=n','d2'],
                   stdin = inp,
                   want = 'stdout',
                   verb = verb)
        d2 = float(d2) * (n2-1)/float(w-1)
        out_transp1 = tmp + '3'
        ivlad.send_to_os('sftransp',
                   stdin=out_remap1,
                   stdout=out_transp1,
                   verb=verb)
        if rem2del_junk2:
            ivlad.send_to_os('sfrm',
                       arg=out_remap1,
                       verb=verb)
        if w < n2:
            ready_for_remap_2 = tmp + '4'
            ivlad.send_to_os('sfbandpass',
                       arg='fhi='+str(0.5/d2),
                       stdin=out_transp1,
                       stdout=ready_for_remap_2,
                       verb=verb)
            rem2del_junk4 = True
        else:
            ready_for_remap_2 = out_transp1
            rem2del_junk4 = False
        ready_for_transp2 = tmp + '5'
        ivlad.send_to_os('sfremap1',
                   arg=['n1='+str(w), 'd1='+str(d2)],
                   stdin=ready_for_remap_2,
                   stdout=ready_for_transp2,
                   verb=verb)
        ivlad.send_to_os('sfrm',
                   arg=out_transp1,
                   verb=verb)
        if rem2del_junk4:
            ivlad.send_to_os('sfrm',
                       arg=ready_for_remap_2,
                       verb=verb)
        ivlad.send_to_os('sftransp',
                   stdin=ready_for_transp2,
                   stdout=out,
                   verb=verb)
        ivlad.send_to_os('sfrm',
                   arg=ready_for_transp2,
                   verb=verb)
    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main
