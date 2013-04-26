/* Compute hilbert transform using different methods. 
type=t -> time domain convolution
type=f -> frequency domain multiplication and use FFT
type=m -> closed-form design of maximally flat FIR hilbert transformer
*/
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include "hilbert.h"

int main (int argc, char* argv[])
{
    int n1,n2,i2, n;
    float *data, *hilb, c, a;
    char *type;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);

    if(NULL==(type=sf_getstring("type"))) type="m";
    /* Choosing hilbert transform method, type=t means time domain, 
       type=f means freqency domain, type=m means using more robust algorithm
       the default is type=m */
    if (!sf_getint("order",&n)) n=100;
    /* Hilbert transformer order if type=m */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) if type=m */

    if (!sf_getfloat("phase",&a)) a=90.;
    /* phase shift (in degrees) */
    a *= SF_PI/180.;
    /* convert to radian */

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	switch(type[0]){
	    case 't':
		 hilbert_t(data,hilb,n1);
		 break;
	    case 'f':
		 hilbert_f(data,hilb,n1);
		 break;
	    case 'm':
		 hilbert_m(data,hilb,n1,n,c);
		 break;
	    default:
	 	 sf_error("%s: Unknown Hilbert type", __FILE__);
	}

	sf_floatwrite(hilb,n1,out);
    }

    exit(0);
}


