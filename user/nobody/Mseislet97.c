/* CDF 9/7 biorthogonal seislet transform */
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

#include "seislet97.h"

int main(int argc, char *argv[])
{
    int i1, n1, i2, n2, i3, n3, n12, niter;
    bool inv, adj, unit;
    char *type;
    float *pp, *qq, **ww, *hilb, **dd, eps;
    sf_file in, out, dip;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    pp = sf_floatalloc(n12);
    qq = sf_floatalloc(n12);
    dd = sf_floatalloc2(n1,n2);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */

    if (adj && 0 != niter) {
	ww = sf_floatalloc2(n1,n2);
	hilb = sf_floatalloc(n1);
	sf_weight_init(ww[0]);
	sf_hilbert_init(n1, 6, 1.);
    } else {
	ww = NULL;
	hilb = NULL;
    }

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [biorthogonal] wavelet type, the default is biorthogonal  */

    seislet97_init(n1,n2,inv,unit,eps,type[0]);
    seislet97_set(dd);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(pp,n12,in);
	sf_floatread(dd[0],n12,dip);

	if (adj) {
	    if (0==niter) {
		seislet97_lop(adj,false,n12,n12,qq,pp);
	    } else {
		sf_solver (seislet97_lop,sf_cgstep,
			   n12,n12,qq,pp,niter,"verb",true,"end");
		sf_cgstep_close();

		/* find envelope */
		for (i2=0; i2 < n2; i2++) {
		    sf_hilbert(qq+i2*n1,hilb);
		    for (i1=0; i1 < n1; i1++) {
			ww[i2][i1] = hypotf(qq[i1+i2*n1],hilb[i1]);
		    }
		}

		sf_solver_prec (seislet97_lop,sf_cgstep,sf_weight_lop,n12, 	      
				n12,n12,qq,pp,niter,0.,"verb",true,"end");
		sf_cgstep_close();
	    }
	} else {
	    seislet97_lop(adj,false,n12,n12,pp,qq);
	}
	sf_floatwrite(qq,n12,out);
    }

    exit(0);
}
