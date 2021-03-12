/* Missing data interpolation in 2-D using wavelet transform and compressive sensing. */
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
    int i, niter, n1, n2, n12, i1, i2, i3, n3, iter, ibreg, nbreg, j; 
    float *dd, *dd2=NULL, *d1, *d2, *m1, *m2, *dd3=NULL, *m, perc, ordert, iperc, *mm;
    char *oper, *type;;
    bool verb, *known;
    sf_file in, out, mask=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("perc",&perc)) perc=99.;
    /* percentage for soft-thresholding */

    if (NULL == (oper=sf_getstring("oper"))) oper="shaping";
    /* [shaping,pocs,bregman] method, the default is shaping  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    dd = sf_floatalloc(n12);
    mm = sf_floatalloc(n12);
    known = sf_boolalloc(n12);
    m = sf_floatalloc(n12);

    d1 = sf_floatalloc(n1);
    d2 = sf_floatalloc(n2);
    m1 = sf_floatalloc(n1);
    m2 = sf_floatalloc(n2);

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }
   
    switch (oper[0]) {
	case 's':
	    if (!sf_getfloat("ordert",&ordert)) ordert=1.;
	    /* Curve order for thresholding parameter, default is linear */

	    sf_sharpen_init(n12,perc,0.5);
	    dd2 = sf_floatalloc(n12);
	    break;
	case 'p':
	    if (!sf_getfloat("ordert",&ordert)) ordert=1.;
	    /* Curve order for thresholding parameter, default is linear */

	    sf_sharpen_init(n12,perc,0.5);
	    dd2 = sf_floatalloc(n12);
	    break;
	case 'b':
	    if (!sf_getint("nbreg",&nbreg)) nbreg=1;
	    /* number of iterations for Bregman iteration */
	    
	    sf_sharpen_init(n12,perc,0.5);
	    dd2 = sf_floatalloc(n12);
	    dd3 = sf_floatalloc(n12);
	    break;
	default:
	    sf_error("Unknown operator \"%s\"",oper);
    }

    for (i3=0; i3 < n3; i3++) {
	if (verb) {
	    sf_warning("slice %d of %d",i3+1,n3);
	} else {
	    sf_warning("slice %d of %d;",i3+1,n3);
	}
	sf_floatread(dd,n12,in);

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
		    
		    /* Forward 2-D DWT */
		    for (i2=0; i2 < n2; i2++) {
			/* Vertical direction */
			sf_wavelet_init(n1,true,true,type[0]);
			for (j=0; j < n1; j++) {
			    d1[j] = dd2[i2*n1+j];
			}
			sf_wavelet_lop(false,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    mm[i2*n1+j] = m1[j];
			}
			sf_wavelet_close();
		    } 
		    for (i1=0; i1 < n1; i1++) { 
			/* Horizontal direction */
			sf_wavelet_init(n2,true,true,type[0]);
			for (j=0; j < n2; j++) {
			    d2[j] = mm[j*n1+i1];
			}
			sf_wavelet_lop(false,false,n2,n2,d2,m2);
			for (j=0; j < n2; j++) {
			    mm[j*n1+i1] = m2[j];
			}
			sf_wavelet_close();
		    } /* Forward 2-D DWT end */	

		    /* Thresholding */
		    if(ordert==0.) {
			iperc = perc;
		    } else {
			iperc = perc-((perc-1)*pow(iter,ordert)*1.)/pow(niter,ordert);
			if(iperc<0.) iperc=0.;
		    }
		    /* Thresholding */
		    sf_sharpen_init(n12,iperc,0.5);
		    sf_sharpen(mm);
		    sf_weight_apply(n12,mm);
		    sf_sharpen_close();
		
		    /* Inverse 2-D DWT */
		    for (i1=0; i1 < n1; i1++) { 
			/* Horizontal direction */
			sf_wavelet_init(n2,true,true,type[0]);
			for (j=0; j < n2; j++) {
			    m2[j] = mm[j*n1+i1];
			}
			sf_wavelet_lop(true,false,n2,n2,d2,m2);
			for (j=0; j < n2; j++) {
			    dd2[j*n1+i1] = d2[j];
			}
			sf_wavelet_close();
		    }
		    for (i2=0; i2 < n2; i2++) {
			/* Vertical direction */
			sf_wavelet_init(n1,true,true,type[0]);
			for (j=0; j < n1; j++) {
			    m1[j] = dd2[i2*n1+j];
			}
			sf_wavelet_lop(true,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    dd2[i2*n1+j] = d1[j];
			}
			sf_wavelet_close();
		    } /* Inverse 2-D DWT end */	
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

		    /* Forward 2-D DWT */
		    for (i2=0; i2 < n2; i2++) {
			/* Vertical direction */
			sf_wavelet_init(n1,true,true,type[0]);
			for (j=0; j < n1; j++) {
			    d1[j] = dd2[i2*n1+j];
			}
			sf_wavelet_lop(false,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    mm[i2*n1+j] = m1[j];
			}
			sf_wavelet_close();
		    } 
		    for (i1=0; i1 < n1; i1++) { 
			/* Horizontal direction */
			sf_wavelet_init(n2,true,true,type[0]);
			for (j=0; j < n2; j++) {
			    d2[j] = mm[j*n1+i1];
			}
			sf_wavelet_lop(false,false,n2,n2,d2,m2);
			for (j=0; j < n2; j++) {
			    mm[j*n1+i1] = m2[j];
			}
			sf_wavelet_close();
		    } /* Forward 2-D DWT end */	

		    /* Thresholding */
		    if(ordert==0.) {
			iperc = perc;
		    } else {
			iperc = perc-((perc-1)*pow(iter,ordert)*1.)/pow(niter,ordert);
			if(iperc<0.) iperc=0.;
		    }
		    /* Thresholding */
		    sf_sharpen_init(n12,iperc,0.5);
		    sf_sharpen(mm);
		    sf_weight_apply(n12,mm);
		    sf_sharpen_close();

		    /* Inverse 2-D DWT */
		    for (i1=0; i1 < n1; i1++) { 
			/* Horizontal direction */
			sf_wavelet_init(n2,true,true,type[0]);
			for (j=0; j < n2; j++) {
			    m2[j] = mm[j*n1+i1];
			}
			sf_wavelet_lop(true,false,n2,n2,d2,m2);
			for (j=0; j < n2; j++) {
			    dd2[j*n1+i1] = d2[j];
			}
			sf_wavelet_close();
		    }
		    for (i2=0; i2 < n2; i2++) {
			/* Vertical direction */
			sf_wavelet_init(n1,true,true,type[0]);
			for (j=0; j < n1; j++) {
			    m1[j] = dd2[i2*n1+j];
			}
			sf_wavelet_lop(true,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    dd2[i2*n1+j] = d1[j];
			}
			sf_wavelet_close();
		    } /* Inverse 2-D DWT end */	

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
		for (iter=0; iter < niter; iter++) { /* Outer iteration */
		    
		    if (verb)
			sf_warning("Outer iteration %d of %d",iter+1,niter);
		    
		    for (ibreg=0; ibreg < nbreg; ibreg++) { /* Inner iteration */
			if (verb)
			    sf_warning("Inner iteration %d of %d",ibreg+1,nbreg);

			for (i1=0; i1 < n12; i1++) {
			    if (known[i1]) dd2[i1] = 0.;
			}
			for (i1=0; i1 < n12; i1++) {
			    dd2[i1] += dd3[i1];
			}

			/* Forward 2-D DWT */
			for (i2=0; i2 < n2; i2++) {
			    /* Vertical direction */
			    sf_wavelet_init(n1,true,true,type[0]);
			    for (j=0; j < n1; j++) {
				d1[j] = dd2[i2*n1+j];
			    }
			    sf_wavelet_lop(false,false,n1,n1,d1,m1);
			    for (j=0; j < n1; j++) {
				mm[i2*n1+j] = m1[j];
			    }
			    sf_wavelet_close();
			} 
			for (i1=0; i1 < n1; i1++) { 
			    /* Horizontal direction */
			    sf_wavelet_init(n2,true,true,type[0]);
			    for (j=0; j < n2; j++) {
				d2[j] = mm[j*n1+i1];
			    }
			    sf_wavelet_lop(false,false,n2,n2,d2,m2);
			    for (j=0; j < n2; j++) {
				mm[j*n1+i1] = m2[j];
			    }
			    sf_wavelet_close();
			} /* Forward 2-D DWT end */	
			

			/* Thresholding */
			sf_sharpen(mm);
			sf_weight_apply(n12,mm);

			/* Inverse 2-D DWT */
			for (i1=0; i1 < n1; i1++) { 
			    /* Horizontal direction */
			    sf_wavelet_init(n2,true,true,type[0]);
			    for (j=0; j < n2; j++) {
				m2[j] = mm[j*n1+i1];
			    }
			    sf_wavelet_lop(true,false,n2,n2,d2,m2);
			    for (j=0; j < n2; j++) {
				dd2[j*n1+i1] = d2[j];
			    }
			    sf_wavelet_close();
			}
			for (i2=0; i2 < n2; i2++) {
			    /* Vertical direction */
			    sf_wavelet_init(n1,true,true,type[0]);
			    for (j=0; j < n1; j++) {
				m1[j] = dd2[i2*n1+j];
			    }
			    sf_wavelet_lop(true,false,n1,n1,d1,m1);
			    for (j=0; j < n1; j++) {
				dd2[i2*n1+j] = d1[j];
			    }
			    sf_wavelet_close();
			} /* Inverse 2-D DWT end */	
			
		    } /* End inner interation */

		    for (i1=0; i1 < n12; i1++) {
			if (known[i1]) {
			    dd3[i1]= dd[i1]+dd3[i1]-dd2[i1]; 
			} else {
			    dd3[i1] = 0.;
			}
		    }
		} /* End outer interation */
		
		for (i1=0; i1 < n12; i1++) {
		    dd[i1] = dd2[i1];
		}
		break;
	} 
	sf_floatwrite (dd,n12,out);
    }
    if (!verb) sf_warning(".");
    exit(0);
}

/* 	$Id$	 */
