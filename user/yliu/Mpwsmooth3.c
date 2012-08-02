/* 3-D structural-oriented smoothing using plane-wave spray and weighted stacking. */
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
    int n1,n2,n3, i1,i2,i3, ns2, ns3, ip2, ip3, ip, np2, np3, np; 
    int i4, n4, k2, k3, j2, j3, ud, lr, foldp, foldn, t1, t2, order;
    float eps, ****u, ****w, ***p1, ***p2, ***norm, **cost, *trace;
    float ax, bx, distance, max;
    sf_file inp, out, dip;

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
    /* spray radius (inline) */
    if (!sf_getint("ns3",&ns3)) sf_error("Need ns3=");
    /* spray radius (crossline) */

    np2 = 2*ns2+1;
    np3 = 2*ns3+1;
    np = np2*np3;

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

    cost = sf_floatalloc2(np2,np3);
    for (i3=0; i3 < np3; i3++) {
	for (i2=0; i2 < np2; i2++) {
	    cost[i3][i2] = 1.;
	}
    }

    dijkstra_init(np2,np3,cost,cost);
    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc4(n1,np,n2,n3);
    w = sf_floatalloc4(n1,np,n2,n3);
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (ip=0; ip < np; ip++) {
		for (i1=0; i1 < n1; i1++) {
		    u[i3][i2][ip][i1] = 0.;
		    w[i3][i2][ip][i1] = 0.;
		}
	    }
	}
    }

    p1 = sf_floatalloc3(n1,n2,n3);
    p2 = sf_floatalloc3(n1,n2,n3);
    norm = sf_floatalloc3(n1,n2,n3);
    trace = sf_floatalloc(n1);

    for (i4=0; i4 < n4; i4++) {
	if (verb) fprintf(stderr,"slice %d of %d\n",i4+1,n4);
	sf_floatread(p1[0][0],n1*n2*n3,dip);
	sf_floatread(p2[0][0],n1*n2*n3,dip);

	for (i3=0; i3 < n3; i3++) { 
	    if (verb) sf_warning("[Predictive spraying] crossline slice %d of %d",i3+1,n3);
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

	/* Scaling factor */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    for(ip3=0; ip3 < np3; ip3 ++) {
			for(ip2=0; ip2 < np2; ip2 ++) {
			    w[i3][i2][ip3*np2+ip2][i1] = u[i3][i2][ip3*np2+ip2][i1]-u[i3][i2][ns3*np2+ns2][i1];
			}
		    }
		}
	    }
	}
	max=0.;
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    for(ip=0; ip < np; ip ++) {
			if (max < fabsf(w[i3][i2][ip][i1])) {
			    max = fabsf(w[i3][i2][ip][i1]);
			}
		    }
		}
	    }
	}

	/* Define weights */
	if (bilat) {
	    if (gauss) {
		for (i3=0; i3 < n3; i3++) {
		    if (verb) sf_warning("[Define Gaussian bilateral weights] crossline slice %d of %d",i3+1,n3);
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    for(ip3=0; ip3 < np3; ip3 ++) {
				for(ip2=0; ip2 < np2; ip2 ++) {
				    w[i3][i2][ip3*np2+ip2][i1] = expf(-0.5*((ip2-ns2)*(ip2-ns2)+(ip3-ns3)*(ip3-ns3))/(ax*ax+FLT_EPSILON)) * 
					expf(-0.5*(u[i3][i2][ip3*np2+ip2][i1]-u[i3][i2][ns3*np2+ns2][i1])*(u[i3][i2][ip3*np2+ip2][i1]-u[i3][i2][ns3*np2+ns2][i1])/(bx*bx*max*max+FLT_EPSILON));
				}
			    }
			}
		    }
		}
	    } else {
		for (i3=0; i3 < n3; i3++) {
		    if (verb) sf_warning("[Define triangle bilateral weights] crossline slice %d of %d",i3+1,n3);
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    for(ip3=0; ip3 < np3; ip3 ++) {
				for(ip2=0; ip2 < np2; ip2 ++) {
				    t1 = ip2;
				    t2 = ip3;
				    while (abs(t1-ns2)!=ns2 && abs(t2-ns3)!=ns3) {
					if (0==abs(t1-ns2) && 0==abs(t2-ns3)) {
					    break;
					} else {
					    if ((t1-ns2)>0) {
						t1++;
					    } else {
						if ((t1-ns2)<0) {
						    t1--;
						} 
					    }
					    if ((t2-ns3)>0) {
						t2++;
					    } else {
						if ((t2-ns3)<0) {
						    t2--;
						}
					    }
					}
				    }
				    if (abs(t1-ns2)==ns2) {
					distance = sqrtf(ns2*ns2+(ip3-ns3)*(ip3-ns3));
				    } else {
					distance = sqrtf((t1-ns2)*(t1-ns2)+ns3*ns3);
				    }
				    w[i3][i2][ip3*np2+ip2][i1] = (1. - sqrtf((ip2-ns2)*(ip2-ns2)+(ip3-ns3)*(ip3-ns3))/(distance+FLT_EPSILON)) * 
					expf(-0.5*(u[i3][i2][ip3*np2+ip2][i1]-u[i3][i2][ns3*np2+ns2][i1])*(u[i3][i2][ip3*np2+ip2][i1]-u[i3][i2][ns3*np2+ns2][i1])/(bx*bx*max*max+FLT_EPSILON));
				}
			    }
			}
		    }
		}
	    }
	} else if (gauss) {
	    for (i3=0; i3 < n3; i3++) {
		if (verb) sf_warning("[Define Gaussian weights] crossline slice %d of %d",i3+1,n3);
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			for(ip3=0; ip3 < np3; ip3 ++) {
			    for(ip2=0; ip2 < np2; ip2 ++) {
				w[i3][i2][ip3*np2+ip2][i1] = expf(-0.5*((ip2-ns2)*(ip2-ns2)+(ip3-ns3)*(ip3-ns3))/(ax*ax+FLT_EPSILON));
			    }
			}
		    }
		}
	    }
	} else {
	    for (i3=0; i3 < n3; i3++) {
		if (verb) sf_warning("[Define triangle weights] crossline slice %d of %d",i3+1,n3);
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			for(ip3=0; ip3 < np3; ip3 ++) {
			    for(ip2=0; ip2 < np2; ip2 ++) {
				t1 = ip2;
				t2 = ip3;
				while (abs(t1-ns2)!=ns2 && abs(t2-ns3)!=ns3) {
				    if (0==abs(t1-ns2) && 0==abs(t2-ns3)) {
					break;
				    } else {
					if ((t1-ns2)>0) {
					    t1++;
					} else {
					    if ((t1-ns2)<0) {
						t1--;
					    } 
					}
					if ((t2-ns3)>0) {
					    t2++;
					} else {
					    if ((t2-ns3)<0) {
						t2--;
					    }
					}
				    }
				}
				if (abs(t1-ns2)==ns2) {
				    distance = sqrtf(ns2*ns2+(ip3-ns3)*(ip3-ns3));
				} else {
				    distance = sqrtf((t1-ns2)*(t1-ns2)+ns3*ns3);
				}
				w[i3][i2][ip3*np2+ip2][i1] = (1. - sqrtf((ip2-ns2)*(ip2-ns2)+(ip3-ns3)*(ip3-ns3))/(distance+FLT_EPSILON));
			    }
			}
		    }
		}
	    }
	    
	}
	
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    p1[i3][i2][i1] = 0.;
		    norm[i3][i2][i1] = 0.;
		}
	    }
	}
	
	for (i3=0; i3 < n3; i3++) {
	    if (verb) sf_warning("[Smoothing] crossline slice %d of %d",i3+1,n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    foldp = 0;
		    foldn = 0;
		    for(ip=0; ip < np; ip ++) {
			p1[i3][i2][i1] +=u[i3][i2][ip][i1]*w[i3][i2][ip][i1];
			if (0!=u[i3][i2][ip][i1]*w[i3][i2][ip][i1]) foldp++;
		    }
		    for(ip=0; ip < np; ip ++) {
			norm[i3][i2][i1] +=w[i3][i2][ip][i1];
			if (0!=w[i3][i2][ip][i1]) foldn++;
		    }
		    p1[i3][i2][i1] = (p1[i3][i2][i1]*foldn)/(norm[i3][i2][i1]*foldp+FLT_EPSILON);
		}
	    }
	}
	sf_floatwrite(p1[0][0],n1*n2*n3,out);
    }

    exit (0);

}

/* 	$Id$	 */
