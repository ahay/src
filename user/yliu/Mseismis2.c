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
#include <rsfpwd.h>

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n12, i1, i2, i3, n3;
    int iter, cnum, cutting, num, order; 
    float *mm, *dd, *dd2=NULL, *dd3=NULL, **pp, *m=NULL;
    float eps, perc, ordert, iperc, orderc, inum;
    char *type, *oper;
    bool verb, *known, cut;
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

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */

    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    if (NULL == (oper=sf_getstring("oper"))) oper="shaping";
    /* [destruction,construction,shaping,pocs,bregman] method, the default is shaping  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool("cut",&cut)) cut = false;
    /* cutting flag, the default is shaping */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

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

    if (cut) {
	if (!sf_getint("cnum",&cnum)) cnum=n2;
	/* Max cutting in cutting operator, default is n2 */
	
	if (!sf_getfloat("orderc",&orderc)) orderc=1.;
	/* Curve order for cutting operator, default is linear */
	
	seislet_init(n1,n2,true,true,eps,order,type[0]);
	if (cnum > n2 || (cnum-1) <0) sf_error("need cnum in [1,n2].");
    } else {
	if (!sf_getfloat("perc",&perc)) perc=99.;
	/* percentage for soft-thresholding */ 
    }
 
    switch (oper[0]) {
	case 'd':
	    for (i1=0; i1 < n12; i1++) {
		mm[i1] = 0.;
	    }
	    seislet_init(n1,n2,true,true,eps,order,type[0]);
	    break;

	case 'c':
	    seislet_init(n1,n2,true,true,eps,order,type[0]);
	    sf_mask_init(known);
	    break;

	case 's':
	    if (!sf_getfloat("ordert",&ordert)) ordert=1.;
	    /* Curve order for thresholding operator, default is linear */
	    
	    sf_sharpen_init(n12,perc,0.5);
	    seislet_init(n1,n2,true,true,eps,order,type[0]);
	    dd2 = sf_floatalloc(n12);
	    break;

	case 'p':
	    if (!sf_getfloat("ordert",&ordert)) ordert=1.;
	    /* Curve order for thresholding operator, default is linear */

	    sf_sharpen_init(n12,perc,0.5);
	    seislet_init(n1,n2,true,true,eps,order,type[0]);
	    dd2 = sf_floatalloc(n12);
	    break;

	case 'b':
	    sf_sharpen_init(n12,perc,0.5);
	    seislet_init(n1,n2,true,true,eps,order,type[0]);
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
		    if (!known[i]) dd[i] = - dd[i];
		} 
		break;

	    case 'c': 
		/* Seislet construct */
		sf_solver_prec(sf_mask_lop, sf_cgstep, seislet_construct, 
			       n12, n12, n12, dd, dd, niter,
			       0., "verb", verb, "end"); 
		for(i=0; i < n12; i++) {
		    if (!known[i]) dd[i] = - dd[i];
		} 
		break;
		
	    case 's':
		for (i1=0; i1 < n12; i1++) {
		    dd2[i1] = dd[i1];
		}
		for (iter=0; iter < niter; iter++) { /* Outer iteration */
		    if (verb)
			sf_warning("Shaping iteration %d of %d",iter+1,niter);

		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) dd2[i1] = 0.;
		    }
		    for (i1=0; i1 < n12; i1++) {
			dd2[i1] += dd[i1];
		    }
		    
		    seislet_lop(true,false,n12,n12,mm,dd2);

		    if (cut) {
			/* Cutting */
			if(orderc==0.) {
			    inum = (float)cnum;
			} else {
			    inum = cnum*1. + 
				((n2*1.-cnum*1.)*pow((iter+1),orderc)*1.)/
				pow(niter,orderc) + 1.;
			    if(inum<0.) inum = 0.;
			}
			cutting = (int)inum;
			
			for (num=1; num < cutting; num *= 2) ;
			num /= 2;		
			
			for (i1=0; i1 < n1; i1++) {
			    for (i2=num; i2 < n2; i2++) {
				mm[i2*n1+i1] = 0.;
			    }
			}
		    } else {
			/* Thresholding */
			if(ordert==0.) {
			    iperc = perc;
			} else {
			    iperc = perc-((perc-1)*pow(iter,ordert)*1.)/
				pow(niter,ordert);
			    if(iperc<0.) iperc=0.;
			}
			sf_sharpen_init(n12,iperc,0.5);
			sf_sharpen(mm);
			sf_weight_apply(n12, mm);
			sf_sharpen_close();
		    }

		    seislet_lop(false,false,n12,n12,mm,dd2);		    
		} /* End outer interation */
		
		for (i1=0; i1 < n12; i1++) {
		    dd[i1] = dd2[i1];
		}
		break;

	    case 'p':
		for (i1=0; i1 < n12; i1++) {
		    dd2[i1] = dd[i1];
		}
		for (iter=0; iter < niter; iter++) { /* Outer iteration */
		    
		    if (verb)
			sf_warning("POCS iteration %d of %d",iter+1,niter);

		    seislet_lop(true,false,n12,n12,mm,dd2);

		    if (cut) {
			/* Cutting */
			if(orderc==0.) {
			    inum = (float)cnum;
			} else {
			    inum = cnum*1. + 
				((n2*1.-cnum*1.)*pow((iter+1),orderc)*1.)/
				pow(niter,orderc) + 1.;
			    if(inum<0.) inum = 0.;
			}
			cutting = (int)inum;
			
			for (num=1; num < cutting; num *= 2) ;
			num /= 2;		
			
			for (i1=0; i1 < n1; i1++) {
			    for (i2=num; i2 < n2; i2++) {
				mm[i2*n1+i1] = 0.;
			    }
			}
		    } else {
			/* Thresholding */
			if(ordert==0.) {
			    iperc = perc;
			} else {
			    iperc = perc-((perc-1)*pow(iter,ordert)*1.)/
				pow(niter,ordert);
			    if(iperc<0.) iperc=0.;
			}
			sf_sharpen_init(n12,iperc,0.5);
			sf_sharpen(mm);
			sf_weight_apply(n12, mm);
			sf_sharpen_close();
		    }
		    
		    seislet_lop(false,false,n12,n12,mm,dd2);

		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) dd2[i1] = 0.;
		    }
		    for (i1=0; i1 < n12; i1++) {
			dd2[i1] += dd[i1];
		    }
		} /* End outer interation */
		
		for (i1=0; i1 < n12; i1++) {
		    dd[i1] = dd2[i1];
		}
		break;

	    case 'b':
		for (i1=0; i1 < n12; i1++) {
		    dd2[i1] = dd[i1];
		    dd3[i1] = dd[i1];
		}
		for (iter=0; iter < niter; iter++) { /* Bregman iteration */
		    if (verb)
			sf_warning("Bregman iteration %d of %d",iter+1,niter);
		    
		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) {
			    dd2[i1] = dd3[i1];
			}
		    }
		    
		    seislet_lop(true,false,n12,n12,mm,dd2);
		    
		    if (cut) {
			/* Cutting */
			if(orderc==0.) {
			    inum = (float)cnum;
			} else {
			    inum = cnum*1. + 
				((n2*1.-cnum*1.)*pow((iter+1),orderc)*1.)/
				pow(niter,orderc) + 1.;
			    if(inum<0.) inum = 0.;
			}
			cutting = (int)inum;
			
			for (num=1; num < cutting; num *= 2) ;
			num /= 2;		
			
			for (i1=0; i1 < n1; i1++) {
			    for (i2=num; i2 < n2; i2++) {
				mm[i2*n1+i1] = 0.;
			    }
			}
		    } else {
			/* Thresholding */
			sf_sharpen(mm);
			sf_weight_apply(n12, mm);
		    }
		    
		    seislet_lop(false,false,n12,n12,mm,dd2);
		    
		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) {
			    dd3[i1]= dd[i1]+dd3[i1]-dd2[i1]; 
			}
		    }
		} /* End Linear Bregman interation */
		
		for (i1=0; i1 < n12; i1++) {
		    dd[i1] = dd2[i1];
		}
		break;
	} 
	sf_cgstep_close();
	sf_floatwrite (dd,n12,out);
    }

    exit(0);
}

/* 	$Id$	 */
