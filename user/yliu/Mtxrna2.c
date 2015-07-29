/* 2D space-noncausal t-x nonstationary regularized autoregression. */
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
    int i1, i2, ii, ix, it, i4, n4;
    float *d, *f, *g, mean;
    char key[6];
    sf_file mat, pre;
    bool verb;

    sf_init(argc,argv);

    mat = sf_input("in");
    pre = sf_output("out");

    mdim = sf_filedims(mat,m);
    n4 = sf_leftsize(mat,2);

    if (2 < mdim) mdim = 2;

    if (!sf_getints("a",a,mdim)) sf_error("Need a=");

    a[1]=a[1]*2;
    if (mdim < 2) sf_error("Need at least two dimension");

    for (j=0; j < mdim; j++) {
	n[j] = m[j];
    }
    n[mdim] = 1;
    for (j=0; j < mdim; j++) {
	n[mdim] *= a[j];
    }

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
    g = sf_floatalloc(nd);

    sf_multidivn_init(ns, mdim, nd, m, rect, d, NULL, verb); 

    for (i4=0; i4 < n4; i4++) {
	sf_warning("slice %d of %d;",i4+1,n4);
	sf_floatread(g,nd,mat);
	ii=0;
	for(ix=-a[1]/2; ix < a[1]/2+1; ix++) {
	    for(it=-a[0]/2; it < (a[0]+1)/2; it++) {
		if(ix!=0) {
		    for (i2=0; i2 < m[1]; i2++) {
			for (i1=0; i1 < m[0]; i1++) { 
			    if ((i2-ix)<0 || 
				(i2-ix)>=m[1] || 
				(i1-it)<0 || 
				(i1-it)>=m[0]) {
				d[ii*m[1]*m[0]+i2*m[0]+i1]=0.;
			    } else {
				d[ii*m[1]*m[0]+i2*m[0]+i1]=
				    g[(i2-ix)*m[0]+(i1-it)];
			    }
			}
		    }
		    ii++;
		}
	    }
	}

	if(ns!=ii) sf_error("Wrong dimension");

	mean = 0.;
	for(i=0; i < n12; i++) {
	    mean += d[i]*d[i];
	}

	mean = sqrtf (mean/n12);
	
	for(i=0; i < n12; i++) {
	    d[i] /= mean;
	}
	for(i=0; i < nd; i++) {
	    g[i] /= mean;
	}

	sf_multidivn (g,f,niter);
	
	
	for(i=0; i < n12; i++) {
	    d[i] *= mean;
	}
	sf_weight2_lop(false,false,n12,nd,f,g);
	sf_floatwrite(g,nd,pre);
	
    }
  
    exit(0);
}

/* 	$Id: Mtxrna2.c 11561 2014-01-05 05:16:51Z yang_liu $	 */
