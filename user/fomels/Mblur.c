/* 2-D blurring and deblurring */
/*
  Copyright (C) 2005 University of Texas at Austin
   
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
    bool adj, inv, spk;
    int n, n1, n2, niter, rect, nc, ic, nr;
    float *x, *y, *p, eps, perc, fact, rect2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n = n1*n2;

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */
    if (!sf_getbool("spk",&spk)) spk=true;
    /* spiky inversion */

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* scaling */

    if (!sf_getint("rect",&rect)) sf_error("Need rect=");
    /* blurring radius */

    if (!sf_getfloat("rect2",&rect2)) rect2=1.0;
    /* smoothing radius */

    if (!sf_getint("repeat",&nr)) nr=1;
    /* repeat smoothing */

    if (!sf_getint("ncycle",&nc)) nc=1;
    /* number of nonlinear cycles */
    
    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (!sf_getfloat("fact",&fact)) fact=0.5;
    /* factor for sharpening */
    
    sf_sharpen_init(n,perc,fact);

    x = sf_floatalloc(n);
    y = sf_floatalloc(n);
    if (adj && inv) {
	p = sf_floatalloc(n);
    } else {
	p = NULL;
    }

    sf_triangle2_init(rect,rect,n1,n2,nr);
    sf_floatread(x,n,in);

    if (adj) {
	if (inv) {
	    if (spk) {
		sf_conjgrad_init(n,n,n,n,eps,1.e-6,true,false);
		
		for (ic=0; ic < nc; ic++) {
		    sf_conjgrad(NULL,sf_triangle2_lop,sf_weight_lop,p,y,x,niter);
		    sf_sharpen(y);
		} 
	    } else {
		sf_impl2_init (rect2, rect2, n1, n2, 0.1, 50., 
			       false, true, NULL, 1, NULL);
		sf_solver_prec (sf_triangle2_lop,sf_cgstep,sf_impl2_lop,
				n,n,n,y,x,niter,eps,"verb",true,"end");
		sf_cgstep_close();
	    }
	} else {
	    sf_triangle2_lop(true,false,n,n,y,x);
	} 
    } else {
	sf_triangle2_lop(false,false,n,n,x,y);
    }

    sf_floatwrite(y,n,out);

    exit(0);
}

