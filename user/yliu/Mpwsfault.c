/* Fault detection from plane-wave spray. */
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
#include <rsfpwd.h>

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, i1,i2,i3, is, ns, ns2, ip, fold, niter, nit, order;
    float eps, perc,fact, ***u, **p, *trace, **xk, **yk;
    sf_file inp, out, dip;
    char *type;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(dip,2);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* spray radius */
    ns2 = 2*ns+1;

    if (NULL == (type=sf_getstring("type"))) type="difference";
    /* [difference,sharpen_similarity] calculation type, the default is difference  */


    switch (type[0]) {
	case 'd':
	    break;
	case 's':
	    if (!sf_getint("niter",&niter)) niter=20;
	    /* number of iterations */
	    
	    if (!sf_getfloat("perc",&perc)) perc=98.;
	    /* percentage for sharpen, default is 98*/

	    if (!sf_getfloat("fact",&fact)) fact=0.5;
	    /* factor for sharpen */
	    
	    if (perc < 0. || perc > 100.)  sf_error("Need perc in range [0.0,100.0]"); 
	    sf_sharpen_init(n1*n2,perc,fact);

	    break;
	default:
	    sf_error("Unknown operator \"%s\"",type);
    }

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc3(n1,ns2,n2);
    for (i2=0; i2 < n2; i2++) {
	for (is=0; is < ns2; is++) {
	    for (i1=0; i1 < n1; i1++) {
		u[i2][is][i1] = 0.;
	    }
	}
    }

    p = sf_floatalloc2(n1,n2);
    xk = sf_floatalloc2(n1,n2);
    yk = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	if (verb) fprintf(stderr,"slice %d of %d\n",i3+1,n3);
	sf_floatread(p[0],n1*n2,dip);

	for (i2=0; i2 < n2; i2++) { /* loop over traces */
	    sf_floatread(trace,n1,inp);

	    for (i1=0; i1 < n1; i1++) {
		u[i2][ns][i1] = trace[i1];
	    }

	    /* predict forward */
	    for (is=0; is < ns; is++) {
		ip = i2-is-1;
		if (ip < 0) break;
		predict_step(false,false,trace,p[ip]);
		for (i1=0; i1 < n1; i1++) {
		    u[ip][ns-is-1][i1] = trace[i1];
		}
	    }

	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = u[i2][ns][i1];
	    }

	    /* predict backward */
	    for (is=0; is < ns; is++) {
		ip = i2+is+1;
		if (ip >= n2) break;
		predict_step(false,true,trace,p[ip-1]);
		for (i1=0; i1 < n1; i1++) {
		    u[ip][ns+is+1][i1] = trace[i1];
		}
	    }

	}

	for(is=0; is < n2; is++) {
	    for(i1=0; i1< n1;i1++) {
		p[is][i1]=0.;
	    }
	}
	for(is=0; is < n2; is++) {
	    for(i1=0; i1 < n1; i1++) {
		fold = 0;
		for(ip=0; ip < ns; ip ++) {
		    p[is][i1] +=u[is][ip][i1];
		    if (0!=u[is][ip][i1]) fold++;
		}
		for(ip=ns+1; ip < 2*ns+1; ip ++) {
		    p[is][i1] +=u[is][ip][i1];
		    if (0!=u[is][ip][i1]) fold++;
		}
		p[is][i1] = p[is][i1]/(fold+FLT_EPSILON);
	    }
	}

	switch (type[0]) {
	    case 'd':
		for(is=0; is < n2; is++) {
		    for(i1=0; i1 < n1; i1++) {
			p[is][i1] = (p[is][i1] - u[is][ns][i1])*
			    (p[is][i1] - u[is][ns][i1]);
		    }
		}
		break;
	    case 's':
		for(is=0; is < n2; is++) {
		    for(i1=0; i1 < n1; i1++) {
			xk[is][i1] = 0.;
			yk[is][i1] = 0.;
		    }
		}
		for(nit=0; nit < niter; nit++) {
		    if (verb) sf_warning("Iteration %d of %d",nit+1,niter);

		    for(is=0; is < n2; is++) {
			for(i1=0; i1 < n1; i1++) {
			    xk[is][i1] = ((p[is][i1]-u[is][ns][i1])*u[is][ns][i1] + eps*eps*xk[is][i1])/(u[is][ns][i1]*u[is][ns][i1]+eps*eps);
			    yk[is][i1] = ((u[is][ns][i1]-p[is][i1])*p[is][i1] + eps*eps*yk[is][i1])/(p[is][i1]*p[is][i1]+eps*eps);
			}
		    }
		    sf_sharpen(xk[0]);
		    sf_weight_apply(n1*n2,xk[0]);
		    sf_sharpen(yk[0]);
		    sf_weight_apply(n1*n2,yk[0]);
		}
		for(is=0; is < n2; is++) {
		    for(i1=0; i1 < n1; i1++) {
			
			p[is][i1] = (1+xk[is][i1])*(1+yk[is][i1]);
		    }
		}
		break;
	}

	sf_floatwrite(p[0],n1*n2,out);
    }	    

    exit (0);
}

/* 	$Id$	 */
