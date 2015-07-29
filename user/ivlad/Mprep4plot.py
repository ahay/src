#! /usr/bin/env python
'''Resamples a 2-D dataset to the desired picture resolution, with antialias
Only one of the h and w parameters needs to be specified.
If prar=n, no action will be taken on axis for which h/w was not specified
If prar=y and only one par (h or w) is specified, the picture will scale
along both axes until it is of the specified dimension.'''

# Copyright (C) 2007-2010 Ioan Vlad
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

import os, sys, tempfile
import rsf.path

try: # Give precedence to local version
    import ivlad, m8rex, sf
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex
    import rsf.user.sf as sf

###############################################################################

def main(par):

    inp = par.string('inp') # input file
    out = par.string('out') # output file

    verb = par.bool('verb', False) # if y, print system commands, outputs

    ivlad.chk_file_dims(inp, 2)

    n1 = ivlad.getout('sfget',['parform=n','n1'], inp, verb)
    n2 = ivlad.getout('sfget',['parform=n','n2'], inp, verb)
    n1 = int(n1)
    n2 = int(n2)

    h = par.int('h', 768)  # output height
    w = par.int('w', 1024) # output width
    has_h = par.string('h')
    has_w = par.string('w')

    if has_h and has_w: # both w and h were read. Check for sanity:
        ivlad.chk_param_limit(w, 'w')
        ivlad.chk_param_limit(h, 'h')
        if (h,w) == (n1,n2):
            ivlad.msg('Change h or w if you want out != inp')
            sf.cp(inp, out, verb)
            return ivlad.unix_success

    # Transform h and w to pixels, if they are not
    # No default value for par.string -- Quirk of rsf, replicated in rsfbak
    unit = par.string('unit', 'px') # unit of h and w. Can be: px, mm, cm, in
    ivlad.chk_par_in_list(unit,['mm','cm','in','px'])

    if unit != 'px':
        ppi = par.int('ppi') # outp. resolution (px/in). Necessary when unit!=px
        ivlad.chk_param_limit(ppi, 'ppi')
        # Transform w and h to px
        if unit == 'in':
            scale = 1
        elif unit == 'mm':
            scale = 254
        elif unit == 'cm':
            scale = 25.4
        w *= ppi / float(scale)
        h *= ppi / float(scale)
        # Don't worry, we'll convert to int after prar block
        del scale

    # Now h and w are in pixels

    # If prar=y, then h and/or w define a bounding box.
    # Find the dimensions of the image inside this box
    prar = par.bool('prar', True) # if y, PReserve Aspect Ratio of input
    if prar: # preserve aspect ratio
        if has_h and not has_w:
            w = n2 * float(h) / n1
        elif has_w and not has_h:
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

    assert h > 1
    assert w > 1

    h = ivlad.valswitch(h, n1, None)
    w = ivlad.valswitch(w, n2, None)
    tmp = tempfile.mktemp(dir=rsf.path.datapath())

    # Interpolation and, if needed, bandpass 
    if h != None:
        d1 = ivlad.getout('sfget', ['parform=n','d1'], inp, verb)
        d1 = float(d1) * (n1-1)/float(h-1)
        if h < n1:
            ready_for_remap_1 = tmp + '1'
            sf.bandpass(inp, ready_for_remap_1, fhi=0.5/d1, verb=verb)
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
        sf.remap1(ready_for_remap_1, out_remap1, n1=h, d1=d1, verb=verb)
        if rem2del_junk1:
            sf.rm(ready_for_remap_1, verb)
    else: # no action on axis 1
        out_remap1 = inp
        rem2del_junk2 = False

    if w != None:
        d2 = ivlad.getout('sfget', ['parform=n','d2'], inp, verb)
        d2 = float(d2) * (n2-1)/float(w-1)
        out_transp1 = tmp + '3'
        sf.transp(out_remap1, out_transp1, verb=verb)
        if rem2del_junk2:
            sf.rm(out_remap1, verb)
        if w < n2:
            ready_for_remap_2 = tmp + '4'
            sf.bandpass(out_transp1, ready_for_remap_2, fhi=0.5/d2, verb=verb)
            rem2del_junk4 = True
        else:
            ready_for_remap_2 = out_transp1
            rem2del_junk4 = False
        ready_for_transp2 = tmp + '5'
        sf.remap1(ready_for_remap_2,ready_for_transp2,n1=w,d1=d2,verb=verb)
        sf.rm(out_transp1, verb)
        if rem2del_junk4:
            sf.rm(ready_for_remap_2, verb)
        sf.transp(ready_for_transp2, out, verb=verb)
        sf.rm(ready_for_transp2, verb)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main, ['inp','out'])
