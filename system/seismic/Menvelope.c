/* Compute data envelope or phase rotation. 

November 2011 program of the month:
http://ahay.org/rsflog/index.php?/archives/274-Program-of-the-month-sfenvelope.html
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

int main (int argc, char* argv[])
{
    bool hlb;
    int n1,n2, i1,i2, n;
    float *data, *hilb, c, a=90;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);

    if (!sf_getint("order",&n)) n=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    if (!sf_getbool("hilb",&hlb)) hlb=false;
    /* if y, compute Hilbert transform */
    
    if (hlb && !sf_getfloat("phase",&a)) a=90.;
    /* phase shift (in degrees) to use with hilb=y */
    a *= SF_PI/180.;
    /* convert to radian */

    sf_hilbert_init(n1, n, c);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	sf_hilbert(data,hilb);

	if (!hlb) {
	    for (i1=0; i1 < n1; i1++) {
		hilb[i1] = hypotf(data[i1],hilb[i1]);
	    }
	} else {
	    for (i1=0; i1 < n1; i1++) {
		hilb[i1] = data[i1]*cosf(a) + hilb[i1]*sinf(a);
	    }
	}

	sf_floatwrite(hilb,n1,out);
    }

    exit(0);
}

/* 	$Id: Menvelope.c 12370 2014-05-04 11:49:34Z sfomel $	 */
