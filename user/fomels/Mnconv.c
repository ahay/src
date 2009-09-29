/* Non-stationary convolution */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
    bool adj;
    int it, nt, itau, ntau, i2, n2, lag;
    float *filt, *p, *q;
    sf_file inp, flt, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    flt = sf_input("filt");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(flt,"n2",&ntau)) sf_error("No n2= in flt");
    n2 = sf_leftsize(inp,1);

    if (!sf_getint("lag",&lag)) lag=(ntau+1)/2;
    /* filter lag */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    p = sf_floatalloc(nt);
    q = sf_floatalloc(nt);
    filt = sf_floatalloc(nt);
    
    for (i2=0; i2 < n2; i2++) {
	if (adj) {
	    sf_floatread(q,nt,inp);
	    for (it=0; it < nt; it++) {
		p[it] = 0.;
	    }
	} else {
	    sf_floatread(p,nt,inp);
	    for (it=0; it < nt; it++) {
		q[it] = 0.;
	    }
	}

	for (itau=0; itau < ntau; itau++) {
	    sf_floatread(filt,nt,flt);
	    
	    for (it = ntau-lag; it <= nt-lag; it++) {
		if (adj) {
		    p[it-itau+lag-1] += q[it]*filt[it];
		} else {
		    q[it] += p[it-itau+lag-1]*filt[it];
		} 
	    }
	}

	if (adj) {
	    sf_floatwrite(p,nt,out);
	} else {
	    sf_floatwrite(q,nt,out);
	}
    }

    exit(0);
}
