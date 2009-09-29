/* Mono-frequency wavelet estimation.
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

#include <float.h>
#include <math.h>
#include <rsf.h>
#include "monof.h"

int main(int argc, char* argv[])
{
    int n2, i2, nk, ik, niter, i0;
    float k0, dk, k, a0, f, *data=NULL, a;

    bool verb;
    sf_file in=NULL, out=NULL, ma=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    ma = sf_output("ma");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");

    if (!sf_getfloat("a0",&a0)) a0=1.;
    /* starting sharpness */
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(ma,"n1",2);
    sf_putint(ma,"nk",nk);
    sf_putfloat(ma,"dk",dk);
    sf_putfloat(ma,"k0",k0);
    sf_fileflush(ma,in);

    data = sf_floatalloc(nk);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,nk,in);

	a = 0.;
	i0 = 0;
	for (ik=0; ik < nk; ik++) {
	    f = data[ik];
	    if (f > a) {
		a = f;
		i0 = ik;
	    }
	}

	a = monof(data,i0,niter,a0,nk,2.*SF_PI*dk,verb);

	k = (float) i0;

	sf_floatwrite(&a,1,ma);
	sf_floatwrite(&k,1,ma);

	sf_floatwrite (data,nk,out);
    }

    exit (0);
}
