/* Structural-oriented smoothing using plane-wave spray and weighted stacking. */
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
    bool verb, bilat, gauss;
    int n1,n2,n3, i1,i2,i3, is, ns, ns2, ip, foldp, foldn, order;
    float eps, ***u, ***w, **p, **norm, *trace, ax, bx, max;
    sf_file inp, out, dip;

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

    if (!sf_getbool("bilat",&bilat)) bilat=false;
    /* if y, bilateral smoothing */

    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, gaussian weight; otherwise, triangle weight */

    if (gauss) {
	if (!sf_getfloat("ax",&ax)) sf_error("Need ax=");
	/* Gaussian weight for the range distance */
    }

    if (bilat) {
	if (!sf_getfloat("bx",&bx)) sf_error("Need bx=");
	/* exponential weight for the domain distance */
    }

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc3(n1,ns2,n2);
    w = sf_floatalloc3(n1,ns2,n2);
    for (i2=0; i2 < n2; i2++) {
	for (is=0; is < ns2; is++) {
	    for (i1=0; i1 < n1; i1++) {
		u[i2][is][i1] = 0.;
		w[i2][is][i1] = 0.;
	    }
	}
    }

    p = sf_floatalloc2(n1,n2);
    norm = sf_floatalloc2(n1,n2);
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

	/* Scaling factor */
	for(is=0; is < n2; is++) {
	    for(i1=0; i1 < n1; i1++) {
		for(ip=0; ip < 2*ns+1; ip ++) {
		    w[is][ip][i1] = u[is][ip][i1]-u[is][ns][i1];
		}
	    }
	}
    
	max=0.;
	for(is=0; is < n2; is++) {
	    for(i1=0; i1 < n1; i1++) {
		for(ip=0; ip < 2*ns+1; ip ++) {
		    if (max < fabsf(w[is][ip][i1])) {
			max = fabsf(w[is][ip][i1]);
		    }
		}
	    }
	}
	
	/* Define weights */
	if (bilat) {
	    if (gauss) {
		for(is=0; is < n2; is++) {
		    for(i1=0; i1 < n1; i1++) {
			for(ip=0; ip < 2*ns+1; ip ++) {
			    w[is][ip][i1] = expf(-0.5*(ip-ns)*(ip-ns)/(ax*ax+FLT_EPSILON)) * expf(-0.5*(u[is][ip][i1]-u[is][ns][i1])*(u[is][ip][i1]-u[is][ns][i1])/(bx*bx*max*max+FLT_EPSILON));
			}
		    }
		}
	    } else {
		for(is=0; is < n2; is++) {
		    for(i1=0; i1 < n1; i1++) {
			for(ip=0; ip < 2*ns+1; ip ++) {
			    w[is][ip][i1] = (1. - fabsf(ip-ns)/(ns+FLT_EPSILON)) * expf(-0.5*(u[is][ip][i1]-u[is][ns][i1])*(u[is][ip][i1]-u[is][ns][i1])/(bx*bx*max*max+FLT_EPSILON));
			}
		    }
		}
	    }
	} else if (gauss) {
	    for(is=0; is < n2; is++) {
		for(i1=0; i1 < n1; i1++) {
		    for(ip=0; ip < 2*ns+1; ip ++) {
			w[is][ip][i1] = expf(-0.5*(ip-ns)*(ip-ns)/(ax*ax+FLT_EPSILON));
		    }
		}
	    }
	} else {
	    for(is=0; is < n2; is++) {
		for(i1=0; i1 < n1; i1++) {
		    for(ip=0; ip < 2*ns+1; ip ++) {
			w[is][ip][i1] = 1. - fabsf(ip-ns)/(ns+FLT_EPSILON);
		    }
		}
	    
	    }
	}

	for(is=0; is < n2; is++) {
	    for(i1=0; i1< n1;i1++) {
		p[is][i1] = 0.;
		norm[is][i1] = 0.;
	    }
	}
	for(is=0; is < n2; is++) {
	    for(i1=0; i1 < n1; i1++) {
		foldp = 0;
		foldn = 0;
		for(ip=0; ip < 2*ns+1; ip ++) {
		    p[is][i1] += u[is][ip][i1]*w[is][ip][i1];
		    if (0 != u[is][ip][i1]*w[is][ip][i1]) foldp++;
		}
		for(ip=0; ip < 2*ns+1; ip ++) {
		    norm[is][i1] +=w[is][ip][i1];
		    if (0 != w[is][ip][i1]) foldn++;
		}
		p[is][i1] = (p[is][i1]*foldn)/(norm[is][i1]*foldp+FLT_EPSILON);
	    }
	}

	sf_floatwrite(p[0],n1*n2,out);
    }	    

    exit (0);
}

/* 	$Id$	 */
