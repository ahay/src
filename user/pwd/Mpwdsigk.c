/* Signal component separation using plane-wave destruction. 

The program works with 2-D data.
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
    
#include "copyk.h"
#include "predk.h"

int main(int argc, char* argv[])
{
    int i, n1, n2, n12, n3, nk, n12k, niter, nliter, iter, i3, order;
    float eps, *d, *s, ***pp, *w=NULL, *p=NULL;
    bool verb;
    sf_file in, out, dips, weight=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    dips = sf_input("dips");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histint(dips,"n1",&i) || i != n1) sf_error("Wrong n1= in dips");
    if (!sf_histint(dips,"n2",&i) || i != n2) sf_error("Wrong n2= in dips");
    if (sf_histint(dips,"n4",&i) && i != n3) sf_error("Wrong n4= in dips");

    if (!sf_histint(dips,"n3",&nk)) nk=1;
    sf_putint (out,"n3",nk);
    sf_putint (out,"n4",n3);
    
    n12 = n1*n2;
    n12k = n12*nk;
    
    if (!sf_getint ("niter",&niter)) niter=50;
    /* maximum number of iterations */

    if (!sf_getint ("nliter",&nliter)) nliter=1;
    /* number of reweighting iterations */

    if (!sf_getfloat ("eps",&eps)) eps=0.;
    /* regularization parameter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    s = sf_floatalloc(n12k);
    d = sf_floatalloc(n12);
    pp = sf_floatalloc3(n1,n2,nk);

    if (nliter > 1) {
	w = sf_floatalloc(n12k);
	p = sf_floatalloc(n12k);

	if (NULL != sf_getstring("weight")) {
	    weight = sf_output("weight"); 
	    sf_putint(weight,"n3",nk);
	    sf_putint(weight,"n4",n3);
	}
    }

    predk_init(nk,n1,n2,0.0001,order,pp);
    copyk_init(nk,n12);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread (d,n12,in);
	sf_floatread (pp[0][0],n12k,dips);
	
	if (1 == nliter) {
	    sf_solver_prec (copyk_lop,sf_cgstep,predk_lop,
			    n12k,n12k,n12,s,d,niter,eps,"verb",verb,"end");
	    sf_cgstep_close();
	} else {
	    for (i=0; i < n12k; i++) {
		w[i] = 1.;
	    }
	    for (iter=0; iter < nliter; iter++) {
		sf_solver_prec (copyk_lop,sf_cgstep,predk_lop,
				n12k,n12k,n12,s,d,niter,eps,
				"verb",verb,"mwt",w,"xp",p,"end");
		sf_cgstep_close();
		
		if (iter < nliter-1) {
		    for (i=0; i < n12k; i++) {
			w[i] = fabsf(p[i]); /* "Cauchy" weight */
		    }	    
		} else {
		    for (i=0; i < n12k; i++) {
			w[i] *= p[i];
		    }
		}
	    }
	    
	    if (NULL != weight) sf_floatwrite(w,n12k,weight);
	}
	
	sf_floatwrite(s,n12k,out);
    }

    exit(0);
}
