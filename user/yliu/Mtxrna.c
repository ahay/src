/* Causal t-x or t-x-y nonstationary regularized autoregression. */
/*
  Copyright (C) 2014 Jilin University
  
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

#include <math.h>

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM], a[SF_MAX_DIM];
    int mdim, nd, ns, n12, i, j, niter;
    int i1, i2, i3, j1, j2, j3, i4, n4;
    float *d, *f, *g, mean, *ff;
    char key[6];
    sf_file flt, mat, pre;
    bool verb;

    sf_init(argc,argv);

    mat = sf_input("in");
    flt = sf_output("out");

    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    mdim = sf_filedims(mat,m);
    n4 = sf_leftsize(mat,3);

    if (3 < mdim) mdim = 3;

    if (!sf_getints("a",a,mdim)) sf_error("Need a=");

    if (mdim < 3) {
	a[2] = 1;
	n[2] = 1;
	m[2] = 1;
    }
    if (mdim < 2) {
	a[1] = 1;
	n[1] = 1;
	m[1] = 1;
    }

    if (1==m[2]) a[2]=1;

    for (j=0; j < mdim; j++) {
	n[j] = m[j];
    }
    n[mdim] =1;
    for (j=0; j < mdim; j++) {
	n[mdim] *= a[j];
    }
/*    ndim = mdim+1; */

    nd = 1;
    for (j=0; j < mdim; j++) {
	nd *= m[j];
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1;
    }

    ns = n[mdim];
    n12 = nd*ns;

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    d = sf_floatalloc(n12);
    f = sf_floatalloc(n12);
    ff = sf_floatalloc(n12);
    g = sf_floatalloc(nd);

    sf_multidivn_init(ns, mdim, nd, m, rect, d, NULL, verb); 

    if (1==m[2] && 1==m[1]) {
	sf_shiftdim(mat,flt,1);
	sf_putint(flt,"n1",a[0]);
    } else if (1==m[2] && 1!=m[1]) {
	sf_shiftdim2(mat,flt,1);
	sf_putint(flt,"n1",a[0]);
	sf_putint(flt,"n2",a[1]);
    } else {
	sf_putint(flt,"n1",a[0]);
	sf_putint(flt,"n2",a[1]);
	sf_putint(flt,"n3",a[2]);
	sf_putint(flt,"n4",m[0]);
	sf_putint(flt,"n5",m[1]);
	sf_putint(flt,"n6",m[2]);
	sf_putint(flt,"n7",n4);

    }

    for (i4=0; i4 < n4; i4++) {
	sf_floatread(g,nd,mat);
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=-a[0]/2; i1 < (a[0]+1)/2; i1++) {
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
				if (1==m[2]) { /* if 2D data */
				    if ((j1+i1)<0 || 
					(j1+i1)>=m[0] ||
					(j2+i2+1)<0 || 
					(j2+i2+1)>=m[1] ||
					(j3+i3)<0 || 
					(j3+i3)>=m[2]) {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				    } else {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 
					    g[(j3+i3)*m[1]*m[0]+
					      (j2+i2+1)*m[0]+
					      (j1+i1)];
				    }
				} else { /* if 3D data */
				    if ((j1+i1)<0 || 
					(j1+i1)>=m[0] ||
					(j2+i2+1)<0 || 
					(j2+i2+1)>=m[1] ||
					(j3+i3+1)<0 || 
					(j3+i3+1)>=m[2]) {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				    } else {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 
					    g[(j3+i3+1)*m[1]*m[0]+
					      (j2+i2+1)*m[0]+
					      (j1+i1)];
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	

	mean = 0.;
	for(i=0; i < n12; i++) {
	    mean += d[i]*d[i];
	}

	if (mean == 0.) {
	    sf_floatwrite(d,n12,flt);
	    continue;
	}
	
	mean = sqrtf (mean/n12);
	
	for(i=0; i < n12; i++) {
	    d[i] /= mean;
	}
	for(i=0; i < nd; i++) {
	    g[i] /= mean;
	}

	sf_multidivn (g,f,niter);
	
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=0; i1 < a[0]; i1++) {
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
				ff[j3*m[1]*m[0]*a[2]*a[1]*a[0]+
				   j2*m[0]*a[2]*a[1]*a[0]+
				   j1*a[2]*a[1]*a[0]+
				   i3*a[1]*a[0]+
				   i2*a[0]+
				   i1] =
				    f[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
				      i2*a[0]*m[2]*m[1]*m[0]+
				      i1*m[2]*m[1]*m[0]+
				      j3*m[1]*m[0]+
				      j2*m[0]+
				      j1];
			    }
			}
		    }
		}
	    }
	}
	
	sf_floatwrite(ff,n12,flt);	
	
	if (pre) {
	    for(i=0; i < n12; i++) {
		d[i] *= mean;
	    }
	    sf_weight2_lop(false,false,n12,nd,f,g);
	    sf_floatwrite(g,nd,pre);
	}
    }
  
    exit(0);
}

/* 	$Id$	 */
