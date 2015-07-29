/* Estimate adaptive nonstationary PEF on aliased traces. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "mask4apef.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM], a[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j, niter;
    int i1, i2, i3, j1, j2, j3, jump, i4, n4;
    float *d, *f, *g, mean, *ff, *outm;
    char key[6];
    sf_file flt, mat, pre, maskin, maskout;
    bool verb, *mask;

    sf_init(argc,argv);

    mat = sf_input("in");
    flt = sf_output("out");

    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    if (!sf_getint("jump",&jump)) jump=2;
    /* Jump parameter */

    mdim = sf_filedims(mat,m);
    
    if (!sf_getint("dim",&ndim)) ndim=mdim; /* number of dimensions */

    n4 = sf_leftsize(mat,ndim);

    if (!sf_getints("a",a,ndim)) sf_error("Need a=");

    if (ndim < 3) {
	a[2] = 1;
	n[2] = 1;
	m[2] = 1;
    }
    if (ndim < 2) {
	a[1] = 1;
	n[1] = 1;
	m[1] = 1;
    }

    for (j=0; j < ndim; j++) {
	n[j] = m[j];
    }
    n[ndim] =1;
    for (j=0; j < ndim; j++) {
	n[ndim] *= a[j];
    }
/*    ndim = mdim+1; */

    nd = 1;
    for (j=0; j < ndim; j++) {
	nd *= m[j];
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1;
    }

    ns = n[ndim];
    n12 = nd*ns;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL != sf_getstring("maskin")) {
	/* optional input mask file */
	mask = sf_boolalloc(nd);
	maskin = sf_input("maskin");
    } else {
	mask = NULL;
	maskin = NULL;
    }
    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	maskout = sf_output("maskout");
	outm = sf_floatalloc(nd);
    } else {
	outm = NULL;
	maskout = NULL;
    }	

    d = sf_floatalloc(n12);
    f = sf_floatalloc(n12);
    ff = sf_floatalloc(n12);
    g = sf_floatalloc(nd);

    sf_multidivn_init(ns, ndim, nd, m, rect, d, NULL, verb); 

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
    	if (NULL != maskin) {
	    sf_floatread(g,nd,maskin);
	    mask4apef (a, jump, m, g, mask);
	    if (NULL != maskout) {
		for (i3=0; i3 < nd; i3++) {
		    if (!mask[i3]) {
			outm[i3] = 0.;
		    } else {
			outm[i3] = 1.;
		    }
		}
		sf_floatwrite (outm,nd,maskout);
	    }
	}

	sf_floatread(g,nd,mat);
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=-a[0]/2; i1 < (a[0]+1)/2; i1++) {
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
				if (0==i2 && 0==i3 && 0 >=i1) {
				    d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
				      i2*a[0]*m[2]*m[1]*m[0]+
				      (i1+a[0]/2)*m[2]*m[1]*m[0]+
				      j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				} else {
				    if ((j1+i1*jump)<0 || (j1+i1*jump)>=m[0] ||
					(j2+i2*jump)<0 || (j2+i2*jump)>=m[1] ||
					(j3+i3*jump)<0 || (j3+i3*jump)>=m[2]) {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				    } else {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  (i1+a[0]/2)*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 
					    g[(j3+i3*jump)*m[1]*m[0]+
					      (j2+i2*jump)*m[0]+
					      (j1+i1*jump)];
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

    	if (NULL != maskin) {
	    for (j3=0; j3 < m[2]; j3++) {
		for (j2=0; j2 < m[1]; j2++) {
		    for (j1=0; j1 < m[0]; j1++) {
			if (!mask[j3*m[1]*m[0]+j2*m[0]+j1]) {
			    for (i3=0; i3 < a[2]; i3++) {
				for (i2=0; i2 < a[1]; i2++) {
				    for (i1=0; i1 < a[0]; i1++) {
					d[i3*a[1]*a[0]*m[2]*m[1]*m[0]+
					  i2*a[0]*m[2]*m[1]*m[0]+
					  i1*m[2]*m[1]*m[0]+
					  j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
				    }
				}
			    }
			    g[j3*m[1]*m[0]+j2*m[0]+j1] = 0.;
			}
		    }
		}
	    } 
	}
	
	sf_multidivn (g,f,niter);
	
	for (i3=0; i3 < a[2]; i3++) {
	    for (i2=0; i2 < a[1]; i2++) {
		for (i1=0; i1 < a[0]; i1++) {
		    for (j3=0; j3 < m[2]; j3++) {
			for (j2=0; j2 < m[1]; j2++) {
			    for (j1=0; j1 < m[0]; j1++) {
				if (0==i2 && 0==i3 && a[0]/2==i1) {
				    ff[j3*m[1]*m[0]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[2]*a[1]*a[0]+
				       j1*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = 1.0;
				} else if (0==i2 && 0==i3 && a[0]/2>i1) {
				    ff[j3*m[1]*m[0]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[2]*a[1]*a[0]+
				       j1*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = 0.0;
				} else {
				    ff[j3*m[1]*m[0]*a[2]*a[1]*a[0]+
				       j2*m[0]*a[2]*a[1]*a[0]+
				       j1*a[2]*a[1]*a[0]+
				       i3*a[1]*a[0]+
				       i2*a[0]+
				       i1] = -1.*
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

/* 	$Id: Mapef.c 12892 2014-06-26 01:10:04Z sfomel $	 */
