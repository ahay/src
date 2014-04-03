/* Cosine taper around the borders (N-D). 

April 2014 program of the month:
http://ahay.org/rsflog/index.php?/archives/381-Program-of-the-month-sfcostaper.html 
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

#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int dim, dim1, nw[SF_MAX_DIM], s[SF_MAX_DIM], n[SF_MAX_DIM];
    int i, j, iw, n1, n2, i0, i2;
    float *data, *w[SF_MAX_DIM], wi;
    char key[4];
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,4,"nw%d",i+1);
	if (!sf_getint(key,nw+i)) nw[i]=0;
	/*( nw#=0 tapering on #-th axis )*/
	if (nw[i] > 0) {
	    dim1 = i;
	    w[i] = sf_floatalloc(nw[i]);
	    for (iw=0; iw < nw[i]; iw++) {
		wi = sinf(0.5*SF_PI*(iw+1.)/(nw[i]+1.));
		w[i][iw] = wi*wi;
	    }
	} else {
	    w[i] = NULL;
	}
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	for (i=0; i <= dim1; i++) {
	    if (nw[i] <= 0) continue;

	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		
		for (iw=0; iw < nw[i]; iw++) {
		    wi = w[i][iw];
		    data[i0+iw*s[i]]          *= wi;
		    data[i0+(n[i]-1-iw)*s[i]] *= wi;
		}
	    }
	}

	sf_floatwrite(data,n1,out);
    }


    exit(0);
}
