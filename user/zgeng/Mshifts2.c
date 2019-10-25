/* Apply linear time shifts on multiple traces. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, n3, i2, i3, s0, ds, pad;
    float *trace, *zero;
    sf_file in, shift;

    sf_init(argc,argv);
    in = sf_input("in");
    shift = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2) || n2 < 2) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getint("s0",&s0)) sf_error("Need s0=");
    /* first shift (in number of samples along 1st axis) */
    if (!sf_getint("ds",&ds)) sf_error("Need ds=");
    /* shift sampling */

    trace = sf_floatalloc(n1);
    zero = sf_floatalloc(n1);
    memset(zero,0,n1);

    for (i3=0; i3 < n3; i3++) {
        for (i2=0; i2 < n2; i2++) {
            sf_floatread(trace,n1,in);
            pad = s0 + ds*i2;
            if (pad < 0 && pad > -n1)
            {
                sf_floatwrite(&trace[-pad],n1+pad,shift);
                sf_floatwrite(zero,-pad,shift);
            }
            else if (pad < n1)
            {
                sf_floatwrite(zero,pad,shift);
                sf_floatwrite(trace,n1-pad,shift);
            }
            else 
                sf_floatwrite(zero,n1,shift);
        }
    }

    exit(0);
}
