/* Non-stationary smoothing by shaping regularization. */
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

int main(int argc, char* argv[])
{
    int i, dim, nm, niter, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    char key[6];
    float *rough, *smooth, *limit, lsum, li;
    sf_file inp, out, lim;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    lim = sf_input("limit");
    /* limiter */

    if (SF_FLOAT != sf_gettype(inp) ||
	SF_FLOAT != sf_gettype(lim)) sf_error("Need float input");
    dim = sf_filedims (inp,n);

    nm = 1;
    for (i=0; i < dim; i++) {
		nm *= n[i];

		if (n[i] > 1) {
			snprintf(key,6,"rect%d",i+1);
			if (!sf_getint(key,rect+i)) rect[i]=1;
			/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
		} else {
			rect[i]=1;
		}
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    
    sf_divn_init(dim,nm,n,rect,niter,true);
    
    rough = sf_floatalloc(nm);
    smooth = sf_floatalloc(nm);
    limit = sf_floatalloc(nm);

    sf_floatread(rough,nm,inp);
    sf_floatread(limit,nm,lim);

    lsum = 0.;
    for (i = 0; i < nm; i++) {
		li = limit[i];
		lsum += li*li;
    }
    lsum = sqrtf (lsum/nm);

    for (i=0; i < nm; i++) {
		limit[i] /= lsum;
		rough[i] *= limit[i];
    }
    
    sf_divn(rough,limit,smooth);

    sf_floatwrite(smooth,nm,out);

    exit(0);
}
