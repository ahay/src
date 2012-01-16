/* None linear Ricker wavelet spectral fit. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include <float.h>
#include <math.h>
#include <rsf.h>
#include <stdio.h>
#include "gaussel.h"

int main(int argc, char* argv[])
{
    int n2, i2, na, ia, ib, ic, ie, id, niter, iter, n, k, i3, i4, nrow, ncol;
    float eps, d2, f, f0, f2, df, di;
    float *m0, *a, *m, *m2, *m3, *e; /*initial frequency*/
    float *data, **r, **rt, **rs, **rp; /*ricker spectrum, r transpose, spectrum related matrix, partial ricker spectrum*/
    float *rd, *r2, *rpd, **rpr, *rpa, *dd, *ap, **arp, *ra, *gamma, **dm, **aparp, **aparpt;
    bool verb;
    sf_file in, out, ma;    
    
    sf_init (argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    ma  = sf_output("ma");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&na)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&df)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&f0)) sf_error("No o1= in input");
    /*number of terms*/
    if (!sf_getint("n",&n)) n = 3;

    m0 = sf_floatalloc(n);
    for (k=1; k < n; k++) {
/*	gu = 125; */
	/*1 Hz per step, up to na steps, equally distributed k steps*/
	if (!sf_getfloat("m",&m0[k])) m0[k] = f0+0.008*(na-1)*df;
    }

    if (!sf_getint("niter",&niter)) niter=100;
    if (!sf_getbool("verb",&verb)) verb=false;

    sf_putint(ma,"n1",2);
    sf_putint(ma,"nf",na);
    sf_putfloat(ma,"df",df);
    sf_putfloat(ma,"f0",f0);
    sf_fileflush(ma,in);

    data = sf_floatalloc(na);

    r = sf_floatalloc2(n,na);
    rt = sf_floatalloc2(na,n);    
    rp = sf_floatalloc2(n,na);  
    rs = sf_floatalloc2(n,n);
    eps = 10.*FLT_EPSILON;
    eps *= eps;

    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);
	sf_floatread(data,na,in);

	d2 = 0.;
	for (ia=0; ia < na; ia++) {
	    di = data[ia];
	    if (di < 0.) data[ia] = -di;
	    d2 += di*di;
	}

	m = m0;
	m2 = sf_floatalloc(n);
	m3 = sf_floatalloc(n);
	e = sf_floatalloc(n);
	a = sf_floatalloc(n);

	for (iter = 0; iter < niter; iter++) {
	    for (i3 = 0; i3 < n; i3++) {
		m[i3] = m0[i3];
		m2[i3] = m[i3]*m[i3];
		m3[i3] = m2[i3]*m[i3];

		rd = sf_floatalloc(n);
		r2 = sf_floatalloc(n);
		rpd = sf_floatalloc(n);
		rpr = sf_floatalloc2(n,n);

	      for (ia = 0; ia < na; ia++) {
		  f = f0 + ia*df;
		  f2 = f*f;
		  e[i3] = exp(-f2/m2[i3]);  
		  r[i3][ia] = e[i3]*f2/m2[i3];
		  rp[i3][ia] = 2.*e[i3]*f2*(f2-m2[i3])/(m3[i3]*m2[i3]);
		  rd[i3] += r[i3][ia]*data[ia];
		  rpd[i3] += rp[i3][ia]*data[ia];
	       }    
	    }
	    /*ricker spectrum matrix transpose*/
	    for (ib = 0; ib < na; ib++) {
		for (i4 = 0; i4 < n; i4++) {
		    rt[ib][i4] = r[i3][ia];
		}
	    }
	    /*ricker spectrum, partial ricker spectrum multiplication */
		       for (ic = 1; ic < n; ic++) {
			   for (id = 1; id < n; id++) {
			       rs[ic][id] = 0;
			       rpr[ic][id] = 0;
			       for (ie = 1; ie < na; ie++) {
				   rs[ic][id] += r[ic][ie]*rt[ie][id];
				   rpr[ic][id] += rp[ic][ie]*rt[ie][id];
			       }
			   }
		       }
		       /*partial ricker spectrum*/
		       for (nrow = 1; nrow < n; nrow++) {
			   for (ncol = 1; ncol < n; ncol++) {
			       if (nrow == ncol)
				     {
					 rpr[nrow][ncol] = 2*rpr[nrow][ncol];
				     }
			   }
		       }
	}
	gaussel_init(n);
	/*getting amplitude*/
	gaussel_solve(rs, rd, a);
    }
    rpa = sf_floatalloc(n);
    for (nrow=1; nrow<n; nrow++) {
	rpa[nrow] = 0;
	for (ncol=1; ncol<n; ncol++) {
	    rpa[ncol] += rpr[nrow][ncol]*a[nrow];
	    }
    }
    ap = sf_floatalloc(n);
    dd = sf_floatalloc(n);
    for (nrow=1; nrow<n; nrow++) {
	dd[nrow] = rpd[nrow]-rpa[nrow]; 
    }
	gaussel_init(n);
	/*solving for partial amplitude*/
	gaussel_solve(rs, dd, ap);

	arp = sf_floatalloc2(n,na);
	/*column is each frequency, dot product*/
	for (ncol=1; ncol<na; ncol++) {
	    for (nrow=1; nrow<n; nrow++) {
		arp[nrow][ncol] = a[nrow]*rp[nrow][ncol];
	    }
	}

	gamma = sf_floatalloc(na);
	ra = sf_floatalloc(na);
	for (nrow=1; nrow<n; nrow++) {
	    for (ia=0; ia<na; ia++) {
		ra[ia] = 0;
		ra[ia] += rt[ia][nrow]*a[nrow];
	    }
	}

	for (ia=0; ia<na; ia++) {
	    gamma[ia] = 0;
	    gamma[ia] = data[ia]-ra[ia];
	}

	aparp = sf_floatalloc2(n,na);
	for (nrow=1; nrow<n; nrow++) {
	    for (ncol=1; ncol<n; ncol++) {
	    aparp[nrow][ncol] = arp[nrow][ncol]+ap[nrow];
	    }
	}
	
	aparpt = sf_floatalloc2(na,n);
	/*aparp transpose*/
	for (ncol=0; ncol<na; ncol++) {
	    for (nrow=1; nrow<n; nrow++) {
		aparpt[ncol][nrow] = aparp[nrow][ncol];
	    }
	}
/*
	dm = sf_floatalloc2(n,na);
	for (ia=0; ia<na; ia++) {
	    gaussel_init(na);
	    gaussel_solve(aparpt[ia], gamma[ia], dm[ia]);
	}
*/

    sf_warning(".");
    exit (0);
}
