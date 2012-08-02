/* Least square imaging condition with pwc regularization. 

*/
/*
  Copyright (C) 2009 China University of Petroleum-Beijing
                 and University of Texas at Austin
  
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
#include <rsfpwd.h>
    
#include "iccopyk.c"


int main(int argc, char* argv[])
{   
    int i, n1, n2, n12, n3, nk, reg, n12k, niter, nliter, iter, i3, order;
    float eps, maxweight, *d, *s, *dwave, ***pp, *w=NULL, *p=NULL, *ww=NULL;
    bool verb, sparse, cut_p;
    sf_file in, out, dips, down, weight=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    dips = sf_input("dips");
    down = sf_input("down");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histint(dips,"n1",&i) || i != n1) sf_error("Wrong n1= in dips");
    if (!sf_histint(dips,"n2",&i) || i != n2) sf_error("Wrong n2= in dips");
    if (sf_histint(dips,"n4",&i) && i != n3) sf_error("Wrong n4= in dips");

    if (!sf_histint(dips,"n3",&nk)) nk=1;
    
    if (!sf_histint(down,"n1",&i) || i != n1) sf_error("Wrong n1= in down");
    if (!sf_histint(down,"n2",&i) || i != n2) sf_error("Wrong n2= in down");
    if (sf_histint(down,"n3",&i) && i != n3) sf_error("Wrong n4= in down");
    
    if (!sf_getbool("sparse",&sparse)) sparse = true;
    /* if sparse = ture   sparse deconvolution cauchy-norm
          if reg = 0: regularization A = |I|
          if reg = 1:  regularization A = |PWD|
       if sparse = false  2-norn deconvolution regularization A = ||I||
    */
    if (!sf_getint("reg",&reg)) reg = 0;
    /* cut off value of precondition */
    

    if (!sf_getbool("cut_p",&cut_p)) cut_p = false;
    /* cut off value of precondition */
    
 
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

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    s = sf_floatalloc(n12k);
    d = sf_floatalloc(n12k);
    pp = sf_floatalloc3(n1,n2,nk);
    dwave = sf_floatalloc(n12k);
    
    if (sparse) {  /* sparse least square imaging condition */

        if (nliter > 1) {
	    w = sf_floatalloc(n12k);
            ww = sf_floatalloc(n12k);
          
	    p = sf_floatalloc(n12k);
          
	    if (NULL != sf_getstring("weight")) {
	        weight = sf_output("weight"); 
	        sf_putint(weight,"n3",nk);
    	        sf_putint(weight,"n4",n3);
	    }
        }
        if (reg == 1){
            predk_init(nk,n1,n2,0.0001,order,pp);
        }
    
        sf_floatread(dwave,n12k,down);
        divk_init(dwave,nk,n2,n1);

        for (i3=0; i3 < n3; i3++) {
	    if (verb) sf_warning("slice %d of %d",i3+1,n3);

	    sf_floatread (d,n12,in);
            if (reg == 1) {
	        sf_floatread (pp[0][0],n12k,dips);
	    }
	    if (1 == nliter) {
                if (reg == 1) {
	            sf_solver_prec (divk_lop,sf_cgstep,predk_lop,
		        	    n12k,n12k,n12,s,d,niter,eps,"verb",verb,"end");
	            sf_cgstep_close();
                }

                if (reg == 0) {
	            sf_solver_prec (divk_lop,sf_cgstep,sf_copy_lop,
		        	    n12k,n12k,n12,s,d,niter,eps,"verb",verb,"end");
	            sf_cgstep_close();
                }
	    } else {
	        for (i=0; i < n12k; i++) {
		    w[i] = 1.;
	        }
	        for (iter=0; iter < nliter; iter++) { sf_warning("nliter:-llllllllllll--%d of %d",iter,nliter);
		    if (reg == 0) {
                        sf_solver_prec (divk_lop,sf_cgstep,sf_copy_lop,
			     	        n12k,n12k,n12,s,d,niter,eps,
		    		        "verb",verb,"mwt",w,"xp",p,"end");
                    }
                    if (reg == 1) {
                       sf_solver_prec (divk_lop,sf_cgstep,predk_lop,
			     	       n12k,n12k,n12,s,d,niter,eps,
		    		       "verb",verb,"mwt",w,"xp",p,"end");
                    }
                    
		    sf_cgstep_close();
		
		    if (iter < nliter-1) {
                        
		        for (i=0; i < n12k; i++) {
			    w[i] = fabsf(p[i]); /* "Cauchy" weight */
                            ww[i] = w[i];
		        }
                        
                        if ((reg == 1) && cut_p) {
                            maxweight = sf_quantile(n12k-20,n12k,ww);
                            /* minweight = sf_quantile(n12k-40,n12k,ww); */
                            for (i=0; i < n12k; i++) {
			        if (w[i] > maxweight) w[i]=maxweight;
				/* if (w[i] > minweight) w[i]=minweight; */
                           }
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

    } else { /* 2-norm least square imaging condition */
        sf_warning("*********sparse==0**********");
        sf_floatread (d,n12,in);
        predk_init(nk,n1,n2,0.0001,order,pp);
        sf_floatread(dwave,n12k,down);
        divk_init(dwave,nk,n2,n1);
        sf_warning("sparse=false");
        sf_solver_reg(divk_lop,sf_cgstep,sf_copy_lop,
                      n12k,n12k,n12,s,d,niter,eps,
                      "verb",verb,"end");
        sf_cgstep_close();
        sf_floatwrite(s,n12k,out);
 
    }

    exit(0);
}
