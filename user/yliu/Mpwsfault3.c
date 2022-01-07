/* 3-D fault detection from plane-wave spray. */
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
    int n1,n2,n3, i1,i2,i3, ns2, ns3, ip, np2, np3, np, order; 
    int i4, n4, k2, k3, j2, j3, ud, lr, fold, niter, nit;
    float eps, perc,fact, ****u, ***p1, ***p2, **cost, *trace, ***xk, ***yk;
    sf_file inp, out, dip;
    char *type;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&n3)) sf_error("No n3= in input");
    n4 = sf_leftsize(inp,3);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("ns2",&ns2)) sf_error("Need ns2=");
    if (!sf_getint("ns3",&ns3)) sf_error("Need ns3=");
    /* spray radius */
    np2 = 2*ns2+1;
    np3 = 2*ns3+1;
    np = np2*np3;

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
	    sf_sharpen_init(n1*n2*n3,perc,fact);

	    break;
	default:
	    sf_error("Unknown operator \"%s\"",type);
    }

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    cost = sf_floatalloc2(np2,np3);
    for (i3=0; i3 < np3; i3++) {
	for (i2=0; i2 < np2; i2++) {
	    cost[i3][i2] = 1.;
	}
    }

    dijkstra_init(np2,np3,cost,cost);
    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc4(n1,np,n2,n3);
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (ip=0; ip < np; ip++) {
		for (i1=0; i1 < n1; i1++) {
		    u[i3][i2][ip][i1] = 0.;
		}
	    }
	}
    }

    p1 = sf_floatalloc3(n1,n2,n3);
    p2 = sf_floatalloc3(n1,n2,n3);
    xk = sf_floatalloc3(n1,n2,n3);
    yk = sf_floatalloc3(n1,n2,n3);
    trace = sf_floatalloc(n1);

    for (i4=0; i4 < n4; i4++) {
	if (verb) fprintf(stderr,"slice %d of %d\n",i4+1,n4);
	sf_floatread(p1[0][0],n1*n2*n3,dip);
	sf_floatread(p2[0][0],n1*n2*n3,dip);

	for (i3=0; i3 < n3; i3++) { 
	    for (i2=0; i2 < n2; i2++) { 
		sf_floatread(u[i3][i2][ns3*np2+ns2],n1,inp);
		dijkstra_source(ns2,ns3);
		
		while (dijskstra_step(&k2,&k3,&ud,&lr)) {
		    
		    /* predict k3,k2 from k3-lr,k2-ud */
		    
		    ip = k3*np2+k2;		    
		    j2 = i2+k2-ns2;
		    j3 = i3+k3-ns3;

		    if (j2 < 0 || j2 >= n2 || 
			j3 < 0 || j3 >= n3 ||
			j2-ud < 0 || j2-ud >= n2 || 
			j3-lr < 0 || j3-lr >= n3) continue;

		    trace = u[j3][j2][ip];		    
		    for (i1=0; i1 < n1; i1++) {
			trace[i1] = u[j3-lr][j2-ud][ip-lr*np2-ud][i1];
		    } 

		    if (ud > 0) {
			predict_step(false,true,trace,p1[j3][j2-ud]);
		    } else if (ud < 0) {
			predict_step(false,false,trace,p1[j3][j2]);
		    }
		    
		    if (lr > 0) {
			predict_step(false,true,trace,p2[j3-lr][j2]);
		    } else if (lr < 0) {
			predict_step(false,false,trace,p2[j3][j2]);
		    }
		}
	    }
	}
	
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    p1[i3][i2][i1] = 0.;
		}
	    }
	}
	
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    fold = 0;
		    for(ip=0; ip < np/2; ip ++) {
			p1[i3][i2][i1] +=u[i3][i2][ip][i1];
			if (0!=u[i3][i2][ip][i1]) fold++;
		    }
		    for(ip=np/2+1; ip < np; ip ++) {
			p1[i3][i2][i1] +=u[i3][i2][ip][i1];
			if (0!=u[i3][i2][ip][i1]) fold++;
		    }
		    p1[i3][i2][i1] = p1[i3][i2][i1]/(fold+FLT_EPSILON);
		}
	    }
	}
	switch (type[0]) {
	    case 'd':
		for(i3=0; i3 < n3; i3++) {
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    p1[i3][i2][i1] = (p1[i3][i2][i1] - u[i3][i2][np/2][i1])*
				(p1[i3][i2][i1] - u[i3][i2][np/2][i1]);
			}
		    }
		}
		break;
	    case 's':
		for(i3=0; i3 < n3; i3++) {
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    xk[i3][i2][i1] = 0.;
			    yk[i3][i2][i1] = 0.;
			}
		    }
		}
		for(nit=0; nit < niter; nit++) {
		    if (verb) sf_warning("Iteration %d of %d",nit+1,niter);

		    for(i3=0; i3 < n3; i3++) {
			for (i2=0; i2 < n2; i2++) {
			    for (i1=0; i1 < n1; i1++) {
				xk[i3][i2][i1] = ((p1[i3][i2][i1]-u[i3][i2][np/2][i1])*u[i3][i2][np/2][i1] + 
						  eps*eps*xk[i3][i2][i1])/(u[i3][i2][np/2][i1]*u[i3][i2][np/2][i1]+eps*eps);
				yk[i3][i2][i1] = ((u[i3][i2][np/2][i1]-p1[i3][i2][i1])*p1[i3][i2][i1] + 
						  eps*eps*yk[i3][i2][i1])/(p1[i3][i2][i1]*p1[i3][i2][i1]+eps*eps);
			    }
			}
		    }
		    sf_sharpen(xk[0][0]);
		    sf_weight_apply(n1*n2*n3,xk[0][0]);
		    sf_sharpen(yk[0][0]);
		    sf_weight_apply(n1*n2*n3,yk[0][0]);
		}
		for(i3=0; i3 < n3; i3++) {
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    p1[i3][i2][i1] = (1+xk[i3][i2][i1])*(1+yk[i3][i2][i1]);
			}
		    }
		}
		break;
	}


	sf_floatwrite(p1[0][0],n1*n2*n3,out);
    }

    exit (0);

}

/* 	$Id$	 */
