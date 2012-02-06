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
    int n2, na, ia, i, j, niter, k, iter, n, l, ib;
    float eps, f, f0, f2, df;
    float *m0=NULL, *a, *m, *m2, *m3, *e, *ap; /*initial frequency*/
    float *data, **rt, **r, **rs, **rp, **rpt; /*ricker spectrum, r transpose, spectrum related matrix, partial ricker spectrum*/
    float *rtd, **rptr, **rtrp, *dd, *ra, *gamma, **rk, *rka, *rptd, *rkd; 
    float **rpa, **rap, **raprpa, **raprpat, **mt, *gt, *dm, d, est, r2;
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
    if (sf_getint("n",&n)) {
	m0 = sf_floatalloc(n);
	if (n==1 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.25*(na-1)*df;
	    }
	} else if (n==2 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.16*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==3 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.13*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==4 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.12*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==5 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.10*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==6 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.09*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==7 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.075*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==8 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.07*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==9 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.06*(i+1)*(na-1)*df;
	    }    
	}
	else if (n==10 && !sf_getfloats("m",m0,n)) {
	    for (i=0; i<n; i++) {
		m0[i] = f0+0.05*(i+1)*(na-1)*df;
	    }    
	}
    } 
    else
    {
	sf_error("n is not specified.");
    }
	    for (i=0; i<n; i++) {
		sf_warning("i=%d m0=%g f0=%g", i, m0[i], f0);
	    }

	    if (!sf_getint("niter",&niter)) niter=100;
	    if (!sf_getbool("verb",&verb)) verb=false;

	    sf_putint(ma,"n1",2);
	    sf_putint(ma,"nf",na);
	    sf_putfloat(ma,"df",df);
	    sf_putfloat(ma,"f0",f0);
	    sf_fileflush(ma,in);

	    data = sf_floatalloc(na);
	    eps = 10.*FLT_EPSILON;
	    eps *= eps;

	    for (i=0; i < n2; i++) {
		sf_warning("slice %d of %d;",i+1,n2);
		sf_floatread(data,na,in);

		m = sf_floatalloc(n);

		for (j = 0; j < n; j++) {
		    m[j] = m0[j];
		}

		m2 = sf_floatalloc(n);
		m3 = sf_floatalloc(n);
		e = sf_floatalloc(n);
		a = sf_floatalloc(n);

		r = sf_floatalloc2(n,na);
		rp = sf_floatalloc2(n,na);
		rt = sf_floatalloc2(na,n);
		rpt = sf_floatalloc2(na,n);

		rtd = sf_floatalloc(n);
		rptd = sf_floatalloc(n);

		for (k=0;k<n;k++) {
		    m[k] = m0[k];
		}

		for (iter = 0; iter < niter; iter++) {
		    for (k = 0; k < n; k++) {
			m2[k] = m[k]*m[k];
			m3[k] = m2[k]*m[k];
			rtd[k] = 0.;
			rptd[k] = 0.;
			for (ia = 0; ia < na; ia++) {
			    f = f0 + ia*df;
			    f2 = f*f;
			    e[k] = exp(-f2/m2[k]);
			    rt[k][ia] = e[k]*f2/m2[k];
			    rpt[k][ia] = 2.*e[k]*f2*(f2-m2[k])/(m3[k]*m2[k]);
			    rtd[k] += rt[k][ia]*data[ia];
			    rptd[k] += rpt[k][ia]*data[ia];
			}
		    }

		    for (i=0;i<n;i++) {
			sf_warning ("rtd of %d =%g",i,rtd[i]);
		    }

		    for (ib = 0; ib < na; ib++) {	    
			for (l = 0; l < n; l++) {
			    r[ib][l] = rt[l][ib];
			    rp[ib][l] = rpt[l][ib];
			}
		    }

		    rs = sf_floatalloc2(n,n);
		    rptr = sf_floatalloc2(n,n);
		    rtrp = sf_floatalloc2(n,n);
		    /*ricker square*/
		    for (k = 0; k < n; k++) {
			for (l = 0; l < n; l++) {
			    rs[k][l] = 0.;
			    for (ib = 0; ib < na; ib++) {
				rs[k][l] += rt[k][ib]*r[ib][l];
			    }
			}
		    }

/* partial ricker spectrum multiplication */

		    for (l = 0; l < n; l++) {
			for (k = 0; k < n; k++) {
			    rptr[l][k] = 0.;
			    rtrp[l][k] = 0.;
			    for (ib = 0; ib < na; ib++) {

				rptr[l][k] += rpt[l][ib]*r[ib][k];
				rtrp[l][k] += rt[l][ib]*rp[ib][k];

			    }
			}
		    }

		    for (k = 0; k < n; k++) {
			for (l = 0; l < n; l++) {
			    if (k == l) {
				rs[k][l] = rs[k][l]+eps;
//			sf_warning ("rs=%g", rs[k][l]);
			    }
			}
		    }

		    gaussel_init(n);

		    gaussel_solve(rs, rtd, a);

		    rk = sf_floatalloc2(n,n);
		    for (k = 0; k < n; k++) {
			for (l = 1; l < n; l++) {
			    rk[k][l] = rptr[k][l]+rtrp[k][l];
			}
		    }

		    rka = sf_floatalloc(n);
		    for (k = 0; k < n; k++) {
			rka[k] = 0;
			for (l = 0; l < n; l++) {
			    rka[k] += rk[k][l]*a[l];
			}
		    }

		    rptd = sf_floatalloc(n);
		    for (k = 0; k < n ; k++) {
			rptd[k] = 0;
			for (ib = 0; ib < na; ib++) {
			    rptd[k] += rpt[k][ib]*data[ib];
			}
		    }

		    rkd = sf_floatalloc(n);
		    for (k = 0; k < n; k++) {
			rkd[k] = rptd[k]-rka[k];
		    }

		    ap = sf_floatalloc(n);

		    gaussel_init(n);
		    gaussel_solve(rs, rkd, ap);

/*DM=inverse((X transpose X)) X transpose Y */
/*aprarp is X; gamma is Y*/
		    dd = sf_floatalloc(n);
		    gamma = sf_floatalloc(na);
		    ra = sf_floatalloc(na);
		    for (ib = 0; ib < na; ib++) {
			for (k = 0; k < n; k++) {
			    ra[ib] = 0;
			    ra[ib] += r[ib][k]*a[k];
			}
		    }

		    for (ib = 0; ib < na; ib++) {
			gamma[ib] = data[ib]-ra[ib];
		    }

		    raprpa = sf_floatalloc2(n,na);
		    rpa = sf_floatalloc2(n,na);
		    rap = sf_floatalloc2(n,na);

		    for (ib = 0; ib < na; ib++) {
			for (k = 0; k < n; k++) {
			    rpa[ib][k] = rp[ib][k]*a[k];
			}
		    }

		    for (ib = 0; ib < na; ib++) {
			for (k = 0; k < n; k++) {
			    rap[ib][k] = r[ib][k]*ap[k];
			}
		    }

		    for (ib = 0; ib < na; ib++) {
			for (k = 0; k < n; k++) {
			    raprpa[ib][k] = rap[ib][k] + rpa[ib][k];
			}
		    }

		    raprpat = sf_floatalloc2(na,n);

		    for (k = 0; k < n; k++) {
			for (ib = 0; ib < na ; ib++) {
			    raprpat[k][ib] = raprpa[ib][k];
			}
		    }

// least squares for delta m.	    
		    mt = sf_floatalloc2(n,n);
		    for (k = 0; k < n; k++) {
			for (l = 0; l < n; l++) {
			    mt[k][l] = 0;
			    for (ib = 0; ib < na; ib++) {
				mt[k][l] += raprpat[k][ib]*raprpa[ib][l];
			    }
			}
		    }

		    gt = sf_floatalloc(n);
		    for (k = 0; k < n; k++) {
			gt[k] = 0;
			for (ib = 0; ib < na; ib++) {
			    gt[k] += raprpat[k][ib]*gamma[ib];
			}
		    }
		    dm = sf_floatalloc(n);
		    gaussel_init(n);
		    /*soliving for amplitude partial*/
		    gaussel_solve(mt, gt, dm);
//estimation, residual
		    est = 0;
		    for (k = 0; k < n; k++) {
			est +=ra[k] ;
		    }
		    d = 0;
		    for (ib = 0; ib < na; ib++) {
			d += data[ib];
		    }
		    r2 = (d - est)*(d - est); 

		    for (k = 0; k < n; k++) {
			m[k] += dm[k];
		    }

		    for (k = 0; k < n; k++) {
			sf_warning("m=%g",m[k]);
		    }

		    if (verb && 5000 > n2) sf_warning("iter=%d r2=%g", iter,r2);

		    if (r2 < eps) break;
		}

		for (k = 0; k < n; k++) {
		    m[k] = fabsf(m[k]);
		    sf_warning("m = %g", m[k]);
		    m2[k] = m[k]*m[k];
		}
      
		sf_floatwrite(m2,n,ma);
		sf_floatwrite(a,n,ma);

		data = sf_floatalloc(na);
		for (ib = 0; ib < na; ib++) {
		    f = f0 + ib*df;
		    f2 = f*f;
		    data[ib] = 0;
		    for (k = 0; k < n; k++) {
			data[ib] += a[k]*exp(-f2/m2[k])*f2/m2[k];
		    }
		}
        
		sf_floatwrite(data,na,out);

	    }
	    sf_warning(".");

	    exit (0);
	}
