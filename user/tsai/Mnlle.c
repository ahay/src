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

/*
remove m, m0, m2, m3, r, rp later.
m1f, m2f: peak frequencies to estimate; r1s, r2s: ricker spectrum 
*/

int main(int argc, char* argv[])
{
    int n2, i2, na, ia, niter, iter;
    float f0, df, f, f2, m, m0, m2, m3, *data, *r, *rp, m1f0, m1f, m1f2, m1f3, m2f0, m2f, m2f2, m2f3;
    float rd, r2, rpd, rp2, rpr, ap, num, den, a, e, dm, eps, di, d2, a1, a2, *r1s, r1s2, *r2s, r2s2;
    float r1sd, r2sd, r1spd, r2spd, r1sp2, r1spr1s, r2sp2, r2spr2s, e1, e2, *r1sp, *r2sp, r1sr2s, r1sr2s2;
    float denpm1, denpm2, r1spr2s, r1sr2sp, r1s2p, r2s2p, pa1m1, pa2m2, pa1m2, pa2m1;
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

/*the following line removed later*/
    if (!sf_getfloat("m",&m0)) m0 = f0+0.25*(na-1)*df;
/*the following two lines added*/
    if (!sf_getfloat("m",&m1f0)) m1f0 = f0+0.25*(na-1)*df;
    if (!sf_getfloat("m",&m2f0)) m2f0 = f0+0.25*(na-1)*df;
    /* initial frequency */
    if (!sf_getint("niter",&niter)) niter = 100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(ma,"n1",2);
    sf_putint(ma,"nf",na);
    sf_putfloat(ma,"df",df);
    sf_putfloat(ma,"f0",f0);
    sf_fileflush(ma,in);

    data = sf_floatalloc(na);
/*the following two lines removed later*/
    r = sf_floatalloc(na);
    rp = sf_floatalloc(na);
/*added*/
    r1s = sf_floatalloc(na);
    r1sp = sf_floatalloc(na);
    r2s = sf_floatalloc(na);
    r2sp = sf_floatalloc(na);

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

/*removed later*/
	m = m0;

/*added*/
	m1f = m1f0;
	m2f = m2f0;

	for (iter = 0; iter < niter; iter++) {
/*removed later*/
	    m2 = m*m;
	    m3 = m2*m;
/*added*/	    
	    m1f2 = m1f*m1f;
	    m1f3 = m1f2*m1f;

	    m2f2 = m2f*m2f;
	    m2f3 = m2f2*m2f;

/*removed later*/	    
	    rd = r2 = rpd = rp2 = rpr = 0.;
/*added*/
	    r1sd = r1s2 = r2sd = r2s2 = r1spd = r2spd = r1sp2 = r1spr1s = r2sp2 = r2spr2s = 0.;

	    for (ia = 0; ia < na; ia++) {
		f = f0 + ia*df;
		f2 = f*f;
/*removed later*/
		e = exp(-f2/m2);
/*added*/
		e1 = exp(-f2/m1f2);
		e2 = exp(-f2/m2f2);
		/*removed later*/
		r[ia] = e*f2/m2;
		rp[ia] = 2.*e*f2*(f2-m2)/(m3*m2);
		/*added*/
		r1s[ia] = e1*f2/m1f2;
		r2s[ia] = e2*f2/m2f2;
		r1sp[ia] = 2.*e1*f2*(f2-m1f2)/(m1f3*m1f2);
		r2sp[ia] = 2.*e2*f2*(f2-m2f2)/(m2f3*m2f2);
		/*removed later*/
		rd += r[ia]*data[ia];
		r2 += r[ia]*r[ia];
		rpd += rp[ia]*data[ia];
		rp2 += rp[ia]*rp[ia];
		rpr += rp[ia]*r[ia];
		/*added*/
		r1sd += r1s[ia]*data[ia];
		r1s2 += r1s[ia]*r1s[ia];
		r2sd += r2s[ia]*data[ia];
		r2s2 += r2s[ia]*r2s[ia];
		r1sr2s += r1s[ia]*r2s[ia];
		r1spd += r1sp[ia]*data[ia];
		r2spd += r2sp[ia]*data[ia];
		r1sp2 += r1sp[ia]*r1sp[ia];
		r2sp2 += r2sp[ia]*r2sp[ia];
		r1spr2s += r1sp[ia]*r2s[ia];
		r1sr2sp += r1s[ia]*r2sp[ia];
		r1s2p += 2*r1s[ia]*r1sp[ia];
		r1spr2s += r1sp[ia]*r2s[ia];
		r2s2p += 2*r2s[ia]*r2sp[ia];
	    }
            r1sr2s2 = r1sr2s*r1sr2s;  
/*removed later*/
	    a = rd/(r2 + eps);
/*added*/
	    a1 = (r1sd*(r2s2+eps)-r1sr2s*r2sd)/((r1s2+eps)*(r2s2+eps)-r1sr2s2);
	    a2 = (r2sd*(r1s2+eps)-r1sr2s*r1sd)/((r1s2+eps)*(r2s2+eps)-r1sr2s2);
/*pa over pm*/
	    denpm1 = denpm2 = ((r1s2+eps)*(r2s2+eps)-r1sr2s2)*((r1s2+eps)*(r2s2+eps)-r1sr2s2);

/*pa1 numerator*/
	    pa1m1 = ((r2s2+eps)*r1spd-r2sd*r1spr2s)*((r1s2+eps)*(r2s2+eps)-r1sr2s2)-(r1sd*(r2s2+eps)-r2sd*r1sr2s)*(r1s2p*(r2s2+eps)-2*r1sr2s*r1spr2s);
/*pa2 numerator*/
	    pa2m2 = ((r1s2+eps)*r2spd-r1sd*r1sr2sp)*((r1s2+eps)*(r2s2+eps)-r1sr2s2)-(r2sd*(r1s2+eps)-r1sd*r1sr2s)*(r2s2p*(r1s2+eps)-2*r1sr2s*r1sr2sp);
/*pa1m2 numerator*/
	    pa1m2 = (r2sp2*r1sd-r2sd*r1sr2sp)*((r1s2+eps)*(r2s2+eps)-r1sr2s2)-(r1sd*(r1s2+eps)-r2sd*r1sr2s)*((r1s2+eps)*r2s2p-2*r1sr2s*r1sr2sp);
/*pa2m1 numerator*/
	    pa2m1 = (r1s2p*r2sd-r1spd*r1spr2s)*(r1s2*r2s2-r1sr2s2)-(r2sd*r1s2-r1sd*r1sr2s)*(r1s2p*r2s2-2*r1sr2s*r1spr2s);

	    ap = (rpd-2.*rpr*a)/(r2 + eps);
	    num =  a*(rpd-rpr*a)+ap*(rd-r2*a);     
	    den = a*a*rp2 + 2.*a*ap*rpr + ap*ap*(r2+eps) + eps;
	    dm = num/den;
        
	    r2 = d2 - 2.*rd*a + r2*a*a;
	    rp2 = dm*dm;

	    if (verb && 5000 > n2) sf_warning("iter=%d r2=%g rp2=%g m=%g a=%g",
					      iter,r2,rp2,m,a);

	    m += dm;
	    if (r2 < eps || rp2 < eps) break;
	}
        
	m = fabsf(m);
	m2 = m*m;
        
	sf_floatwrite(&m2,1,ma);
	sf_floatwrite(&a,1,ma);
        
	for (ia = 0; ia < na; ia++) {
	    f = f0 + ia*df;
	    f2 = f*f;
	    data[ia] = a*exp(-f2/m2)*f2/m2;
	}
        
	if (verb) sf_warning("m=%g a=%g",m,a*m*sqrtf(SF_PI)*0.5);
	if (verb) sf_warning ("%d of %d, %d iterations", i2+1, n2, iter);
        
	sf_floatwrite (data,na,out);
    }

    sf_warning(".");

    exit (0);
}

/* 	tsai $	 */
