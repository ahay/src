/* Nonstationary autoregressive shaping */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "mmmult.h"

int main(int argc, char* argv[])
{
    int i, n1, n2, n12, nf1, nf2, nf3, nf4, nf1234, ir, nr;
    float *data, *tmp1, *tmp2, *filt, eps, r;
    sf_file inp, out, fil;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getint("rect",&nr)) nr=1;
    /* shaping radius */

    /* input data, output model */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    /* input filter */
    if (!sf_histint(fil,"n1",&nf1)) sf_error("No n1= in filter");
    if (!sf_histint(fil,"n2",&nf2)) sf_error("No n2= in filter");
    if (!sf_histint(fil,"n3",&nf3)) sf_error("No n3= in filter");
    if (!sf_histint(fil,"n4",&nf4)) sf_error("No n4= in filter");
    
    if (nf3!=n1 || nf4!=n2) sf_error("need n1==nf3 && n2==nf4");
    nf1234 = nf1*nf2*nf3*nf4;
    
    filt = sf_floatalloc(nf1234);
    data = sf_floatalloc(n12);
    
    tmp1 = sf_floatalloc(n12);
    tmp2 = sf_floatalloc(n12);

    mmmult_init (filt, nf1, nf2, nf3, nf4);	

    sf_floatread(filt,nf1234,fil);
    sf_floatread(data,n12,inp);

    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = -0.5*eps/(r*r);

	for (i=0; i < n12; i++) {
	    tmp1[i] = r*data[i];
	}
	    
	mmmult_lop(false,false,n12,n12,tmp1,tmp2);
	mmmult_lop(true,true,n12,n12,data,tmp2);
    }
	
    sf_floatwrite(data,n12,out);
    
    exit(0);
}

