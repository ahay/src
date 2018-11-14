/* Estimate smoothing radius (min = 1) */
/*
Copyright (C) 2018 University of Texas at Austin

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
#include <math.h>

int main (int argc, char* argv[])
{
    int n1, n1f, n2, n2f, i, n12, n12f;
    float *rect, *fr, maxrad, c, *rad;
    sf_file in, out, freq;
	
    sf_init (argc,argv);

    in = sf_input("in");
    freq = sf_input("freq");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(freq,"n1",&n1f)) sf_error("No n1= in frequency difference.");

    n2 = sf_leftsize(in,1);    
    n2f = sf_leftsize(freq,1);    
    
    n12 = n1*n2;
    n12f = n1f*n2f;

    if (n1 != n1f) sf_error("Need matching n1");
    if (n2 != n2f) sf_error("Need matching n2");

    if (!sf_getfloat("c",&c)) c=1.;
    if (!sf_getfloat("maxrad",&maxrad)) maxrad=1000.;

    rect = sf_floatalloc(n12);
    sf_floatread(rect,n12,in);

    fr = sf_floatalloc(n12f);
    sf_floatread(fr,n12,freq);

    rad = sf_floatalloc(n12);
	
    for (i=0; i < n12; i++) {

        /* update radius */
	    rad[i] = rect[i]+c*fr[i];

        /* set constraint conditions: [1, maxrad] */
        if (rad[i] > maxrad)
                rad[i] = maxrad;
        if (rad[i] < 1.0)
                rad[i] = 1.0;
    }

    sf_floatwrite(rad,n12,out);
    exit(0);
}


