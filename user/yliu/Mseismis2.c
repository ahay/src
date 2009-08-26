/* Missing data interpolation in 2-D using seislet transform. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include "seisletoper.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n12, i1, i2, i3, n3, iter, ibreg, nbreg, cnum, cut, num; 
    float *mm, *dd, *dd2=NULL, *dd3=NULL, **pp, *m=NULL, eps, perc1, perc2, ordert, iperc, orderc, inum;
    char *type, *oper;
    bool verb, *known;
    sf_file in, out, dip, mask=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("perc1",&perc1)) perc1=98.;
    /* percentage for shrinkage and Bregman iteration */


    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    if (NULL == (oper=sf_getstring("oper"))) oper="thresholding";
    /* [destruction,preconditioning,thresholding,shaping,bregman,cutting] method, the default is thresholding  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */


    pp = sf_floatalloc2(n1,n2);
    mm = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    known = sf_boolalloc(n12);
    m = sf_floatalloc(n12);

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }
 
    switch (oper[0]) {
	case 'd':
	    for (i1=0; i1 < n12; i1++) {
		mm[i1] = 0.;
	    }
	    seislet_init(n1,n2,true,true,eps,type[0]);
	    break;
	case 'p':
	    seislet_init(n1,n2,true,true,eps,type[0]);
	    sf_mask_init(known);
	    break;

	case 't':
	    if (!sf_getfloat("perc2",&perc2)) perc2=90.;
	    /* percentage for output in model space shrinkage*/

	    if (!sf_getfloat("ordert",&ordert)) ordert=1.;
	    /* Curve order for thresholding operator, default is linear */
	    
	    sf_sharpen_init(n12,perc1);
	    seislet_init(n1,n2,true,true,eps,type[0]);
	    dd2 = sf_floatalloc(n12);
	    break;

	case 'c':
	    if (!sf_getint("cnum",&cnum)) cnum=n2;
	    /* Max cutting in cutting operator, default is n2 */

	    if (!sf_getfloat("orderc",&orderc)) orderc=1.;
	    /* Curve order for cutting operator, default is linear */

	    seislet_init(n1,n2,true,true,eps,type[0]);
	    dd2 = sf_floatalloc(n12);
	    if (cnum > n2 || (cnum-1) <0) sf_error("need cnum in [1,n2].");
	    break;

	case 's':
	    sf_sharpen_init(n12,perc1);
	    seislet_init(n1,n2,true,true,eps,type[0]);
	    dd2 = sf_floatalloc(n12);
	    break;

	case 'b':
	    if (!sf_getint("nbreg",&nbreg)) nbreg=100;
	    /* number of iterations for Bregman iteration */
	    
	    sf_sharpen_init(n12,perc1);
	    seislet_init(n1,n2,true,true,eps,type[0]);
	    dd2 = sf_floatalloc(n12);
	    dd3 = sf_floatalloc(n12);
	    break;

	default:
	    sf_error("Unknown operator \"%s\"",oper);

    }

    seislet_set(pp);
    
    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);
	sf_floatread(dd,n12,in);
	sf_floatread(pp[0],n12,dip);
	
	if (NULL != mask) {
	    sf_floatread(m,n12,mask);
	    
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (m[i] != 0.);
	    }
	} else {
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (dd[i] != 0.);
	    }
	}
	
        switch (oper[0]) {
	    case 'd':
		/* Seislet destruct */
		sf_solver(seislet_destruct, sf_cgstep, n12, n12, dd, mm, niter,
			  "known", known, "x0", dd, "verb", verb, "end"); 
		for(i=0; i < n12; i++) {
		    if ( !known[i] ) {
			dd[i] = - dd[i];
		    }
		} 
		break;
	    case 'p': 
		/* Seislet construct */
		sf_solver_prec(sf_mask_lop, sf_cgstep, seislet_construct, n12, n12, n12, dd, dd, niter,
			       0., "verb", verb, "end"); 
		for(i=0; i < n12; i++) {
		    if ( !known[i] ) dd[i] = - dd[i];
		} 
		break;
		
	    case 't':
		/* Thresholding for shaping */
		for (iter=0; iter < niter; iter++) {
		    if (verb)
			sf_warning("Model space shrinkage iteration %d of %d",iter+1,niter);
		    seislet_lop(false,false,n12,n12,mm,dd2);
		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) dd2[i1]=dd[i1];
		    }
		    
		    seislet_lop(true,false,n12,n12,mm,dd2);
		    if(ordert==0.) {
			iperc = perc1;
		    } else {
			iperc = perc1-((perc1-1)*pow(iter,ordert)*1.)/pow(niter,ordert);
			if(iperc<0.) iperc=0.;
		    }
		    /* Thresholding */
		    sf_sharpen_init(n12,iperc);
		    sf_sharpen(mm);
		    sf_weight_apply(n12,mm);
		    sf_sharpen_close();
		}
		
		seislet_lop(false,false,n12,n12,mm,dd2);
		for (i1=0; i1 < n12; i1++) {
		    if (known[i1]) dd2[i1]=dd[i1];
		}
		sf_sharpen_init(n12,perc2);
		seislet_lop(true,false,n12,n12,mm,dd2);
		sf_sharpen(mm);
		sf_weight_apply(n12,mm);
		sf_sharpen_close();
		seislet_lop(false,false,n12,n12,mm,dd2);
		
		for (i1=0; i1 < n12; i1++) {
		    if (0.==perc2) {
			if (!known[i1]) dd[i1] = dd2[i1];
		    } else {
			dd[i1] = dd2[i1];
		    }
		}
		
		break;
		
	    case 'c':
		/* Cutting for linear inversion */
		for (iter=0; iter < niter; iter++) {
		    if (verb)
			sf_warning("Model space cutting iteration %d of %d",iter+1,niter);
		    seislet_lop(false,false,n12,n12,mm,dd2);
		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) dd2[i1]=dd[i1];
		    }
		    
		    seislet_lop(true,false,n12,n12,mm,dd2);
		    
		    if(orderc==0.) {
			inum = (float)cnum;
		    } else {
			inum = cnum*1. + ((n2*1.-cnum*1.)*pow((iter+1),orderc)*1.)/pow(niter,orderc) + 1.;
			if(inum<0.) inum = 0.;
		    }
		    cut = (int)inum;
		    
		    for (num=1; num < cut; num *= 2) ;
		    num /= 2;		
		    
		    for (i1=0; i1 < n1; i1++) {
			for (i2=num; i2 < n2; i2++) {
			    mm[i2*n1+i1] = 0.;
			}
		    }
		}
		
		seislet_lop(false,false,n12,n12,mm,dd2);
		
		for (i1=0; i1 < n12; i1++) {
		    if (!known[i1]) dd[i1] = dd2[i1];
		}
		
		break;
		
	    case 's':
		for (i1=0; i1 < n12; i1++) {
		    dd2[i1]= dd[i1];
		}
		for (iter=0; iter < niter; iter++) {
		    if (verb)
			sf_warning("Data space shrinkage iteration %d of %d",iter+1,niter);
		    seislet_lop(true,false,n12,n12,mm,dd2);
		    sf_sharpen(mm);
		    sf_weight_apply(n12,mm);
		    seislet_lop(false,false,n12,n12,mm,dd2);
		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) dd2[i1]= dd[i1];
		    }
		}
		for (i1=0; i1 < n12; i1++) {
		    dd[i1] = dd2[i1];
		}
		break;
		
	    case 'b':
		
		for (i1=0; i1 < n12; i1++) {
		    dd2[i1] = dd[i1];
		}
		for (ibreg=0; ibreg < nbreg; ibreg++) {
		    if (verb)
			sf_warning("Bregman iteration %d of %d",ibreg+1,nbreg);
		    for (i1=0; i1 < n12; i1++) {
			mm[i1]= 0.;
		    }		
		    for (iter=0; iter < niter; iter++) {
			if (verb)
			    sf_warning("Shrinkage iteration %d of %d",iter+1,niter);
			seislet_lop(false,false,n12,n12,mm,dd3);
			for (i1=0; i1 < n12; i1++) {
			    if (known[i1]) dd3[i1]=dd2[i1];
			}
			seislet_lop(true,false,n12,n12,mm,dd3);
			if (iter != 0) {
			    sf_sharpen(mm);
			    sf_weight_apply(n12,mm);
			}
		    }
		    
		    seislet_lop(false,false,n12,n12,mm,dd3);
		    
		    for (i1=0; i1 < n12; i1++) {
			if (!known[i1]) dd3[i1]= 0.;
		    }
		    for (i1=0; i1 < n12; i1++) {
			dd2[i1]= dd[i1]+dd2[i1]-dd3[i1];
		    }
		}
		seislet_lop(false,false,n12,n12,mm,dd2);
		for (i1=0; i1 < n12; i1++) {
		    if (!known[i1]) dd[i1] = dd2[i1];
		}
		
		break;
	} 
	sf_cgstep_close();
	sf_floatwrite (dd,n12,out);
    }
    
    exit(0);
}

/* 	$Id$	 */
