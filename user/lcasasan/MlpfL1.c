/* Local prediction filter (n-dimensional) in L1 norm. */
/*
  Copyright (C) 2010 Politecnico di Milano
  Copyright (C) 2006 University of Texas at Austin
   
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
    bool verb;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], rect[SF_MAX_DIM];
    int ndim, mdim, nd, ns, n12, i, j, liter, niter, na, ia;
    float *den, *f, *data, mean, a0, perc;
    char key[6], *peffile, *lagfile;
    sf_filter aa;
    sf_file dat, flt, mat, pre, pef, lag;

    sf_init(argc,argv);

    dat = sf_input("in");
    mat = sf_input("match");
    flt = sf_output("out");

    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    ndim = sf_filedims(dat,n);
    mdim = sf_filedims(mat,m);

    if (mdim > ndim) 
	sf_error("Wrong dimensions: %d > %d",mdim,ndim);
 
    nd = 1;
    for (j=0; j < mdim; j++) {
	if (m[j] != n[j]) 
	    sf_error("Size mismatch [n%d]: %d != %d",
		     j+1,m[j],n[j]);
	nd *= m[j];
	snprintf(key,6,"rect%d",j+1);
	if (!sf_getint(key,rect+j)) rect[j]=1;
    }
    for (ns = 1; j < ndim; j++) {
	ns *= n[j];
	if (pre) {
	    snprintf(key,6,"n%d",j+1);
	    sf_putint(pre,key,1);
	}
    }
    n12 = nd*ns;

    if (!sf_getint("liter",&liter)) liter=10;
    /* number of CG iterations */

    if (!sf_getint("niter",&niter)) niter=25;
    /* number of POCS iterations [L1]*/

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening [L1]*/


    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (NULL == (peffile = sf_getstring("pef"))) { 
        /* signal PEF file (optional) */
	aa = NULL;
    } else {
	pef = sf_input(peffile);
	if (!sf_histint(pef,"n1",&na)) 
	    sf_error("No n1= in pef");
	aa = sf_allocatehelix (na);
	
	if (!sf_histfloat(pef,"a0",&a0)) a0=1.;
	sf_floatread (aa->flt,na,pef);
	for( ia=0; ia < na; ia++) {
	    aa->flt[ia] /= a0;
	}
	if (NULL != (lagfile = sf_getstring("lag")) 
	    || 
	    NULL != (lagfile = sf_histstring(pef,"lag"))) {
            /* file with PEF lags (optional) */
	    lag = sf_input(lagfile);
	    sf_intread(aa->lag,na,lag);
	    sf_fileclose(lag);
	} else {
	    for( ia=0; ia < na; ia++) {
		aa->lag[ia] = ia+1;
	    }
	}
    }

    den = sf_floatalloc(n12);
    f = sf_floatalloc(n12); 
    data = sf_floatalloc(nd);

    sf_multidivnL1_init(ns, mdim, nd, m, rect, den, aa, perc, verb);
    
    sf_floatread(den,n12,dat);
    sf_floatread(data,nd,mat);

    mean = 0.;
    for(i=0; i < n12; i++) {
	mean += den[i]*den[i];
    }
    if (mean == 0.) {
	sf_floatwrite(den,n12,flt);
	exit(0);
    }

    mean = sqrtf (mean/n12);
    
    for(i=0; i < n12; i++) {
	den[i] /= mean;
    }
    for(i=0; i < nd; i++) {
	data[i] /= mean;
    }
    
    for (i=0;i<n12;i++) f[i]=0.;
    sf_multidivnL1 (data,f,niter,liter); /*data = num*/


    sf_floatwrite(f,n12,flt);

    if (pre) {
	for(i=0; i < n12; i++) {
	    den[i] *= mean;
	}
	
	sf_weight2_lop(false,false,n12,nd,f,data);
	sf_floatwrite(data,nd,pre);
    }

    sf_multidivnL1_close();

    exit(0);
}
