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

/* m1f, m2f, m3f: peak frequencies to estimate; r1s, r2s, r3s: ricker spectrum */

int main(int argc, char* argv[])
{
    int n2, i2, na, ia, niter, iter;
    float f0, df, f, f2, *data, m1f0, m1f, m1f2, m1f3, m2f0, m2f, m2f2, m2f3, m3f0, m3f, m3f2, m3f3, dm1f, dm2f, dm3f;
    float eps, r2, di, d2, a1, a2, a3, a1num, a1den, a2num, a2den; 
    float *r1s, r1s2, *r2s, r2s2, *r3s, r3s2;
    float r1sd, r2sd, r3sd, e1, e2, e3, r1sr2s, r1sr2s2, r1sr3s, r1sr3s2, r2sr3s, r2sr3s2;
    float r1spd, r2spd, r1sp2, r1spr1s, r2sp2, r2spr2s;
    float r1spr2s, r1sr2sp, r1s2p, r2s2p;
    float *r1sp, *r2sp, *r3sp;
    /*  float a1m1, a2m2, a1m2, a2m1;*/
    float numm1, denm1, numm2, denm2;
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


    if (!sf_getfloat("m",&m1f0)) m1f0 = f0+0.05*(na-1)*df;
    if (!sf_getfloat("m",&m2f0)) m2f0 = f0+0.25*(na-1)*df;
    if (!sf_getfloat("m",&m3f0)) m3f0 = f0+0.40*(na-1)*df;
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

    r1s = sf_floatalloc(na);
    r1sp = sf_floatalloc(na);
    r2s = sf_floatalloc(na);
    r2sp = sf_floatalloc(na);
    r3s = sf_floatalloc(na);
    r3sp = sf_floatalloc(na);

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

	m1f = m1f0;
	m2f = m2f0;
	m3f = m3f0;
	for (iter = 0; iter < niter; iter++) {
	    
	    m1f2 = m1f*m1f;
	    m1f3 = m1f2*m1f;
	    m2f2 = m2f*m2f;
	    m2f3 = m2f2*m2f;
	    m3f2 = m3f*m3f;
	    m3f3 = m3f2*m3f;

	    r1sd = r1s2 = r2sd = r2s2 = r3sd = r3s2 = 0.;
	    r1sr2s = r1sr3s = r2sr3s = r1spr2s = r1sr2sp = r1s2p = r2s2p = 0.;
            r1spd = r2spd = r1sp2 = r1spr1s = r2sp2 = r2spr2s = 0.;
	    denm1 = numm1 = denm2 = numm2 = 0.;
	    dm1f = dm2f = dm3f = r2 = 0.;
	    for (ia = 0; ia < na; ia++) {
		f = f0 + ia*df;
		f2 = f*f;

		e1 = exp(-f2/m1f2);
		e2 = exp(-f2/m2f2);
		e3 = exp(-f2/m3f2);

		r1s[ia] = e1*f2/m1f2;
		r2s[ia] = e2*f2/m2f2;
		r3s[ia] = e3*f2/m3f2;
		r1sp[ia] = 2.*e1*f2*(f2-m1f2)/(m1f3*m1f2);
		r2sp[ia] = 2.*e2*f2*(f2-m2f2)/(m2f3*m2f2);
		r3sp[ia] = 2.*e3*f2*(f2-m3f2)/(m3f3*m3f2);

		r1sd += r1s[ia]*data[ia];
		r1s2 += r1s[ia]*r1s[ia];
		r2sd += r2s[ia]*data[ia];
		r2s2 += r2s[ia]*r2s[ia];
		r3sd += r3s[ia]*data[ia];
		r3s2 += r3s[ia]*r3s[ia];	

		r1sr2s += r1s[ia]*r2s[ia];
		r2sr3s += r2s[ia]*r3s[ia];
		r1sr3s += r1s[ia]*r3s[ia];
	    }
            r1sr2s2 = r1sr2s*r1sr2s;
            r1sr3s2 = r1sr3s*r1sr3s;
            r2sr3s2 = r2sr3s*r2sr3s;


	    a1num = (r1sd*(r3s2+eps)-r3sd*r1sr3s)*((r2s2+eps)*(r3s2+eps)-r2sr3s2)-(r1sr2s*(r3s2+eps)-r2sr3s*r1sr3s)*(r2sd*(r3s2+eps)-r3sd*r2sr3s);
	    a1den = ((r1s2+eps)*(r3s2+eps)-r1sr3s2)*((r2s2+eps)*(r3s2+eps)-r2sr3s2)-(r1sr2s*(r3s2+eps)-r1sr3s*r2sr3s)*(r1sr2s*(r3s2+eps)-r2sr3s*r1sr3s);
	    a1 = a1num/a1den;

	    a2num = (r2sd*(r3s2+eps)-r3sd*r2sr3s)-a1*(r1sr2s*(r3s2+eps)-r1sr3s*r2sr3s);
	    a2den = (r2s2+eps)*(r3s2+eps)-r2sr3s2;
	    a2 = a2num/a2den;

	    a3 = (r3sd-a2*r2sr3s-a1*r1sr3s)/(r3s2+eps);

	    if (verb && 5000 > n2) sf_warning("iter=%d r2=%g numm1=%g denm1=%g dm1=%g numm2=%g denm2=%g dm2=%g m1f=%g m2f=%g a1=%g a2=%g",
					      iter,r2,numm1,denm1,dm1f,numm2,denm2,dm2f,m1f,m2f,a1,a2);
	    m1f += dm1f;
            m2f += dm2f;
	    m3f += dm3f;
	    if (r1s2 < eps || r2s2 < eps || r1sp2 < eps || r2sp2 < eps || r3s2 < eps) break;
	}     

	m1f = 6;
	m2f = 28;
	m3f = 60;

	m1f = fabsf(m1f);
	m1f2 = m1f*m1f;
	m2f = fabsf(m2f);
	m2f2 = m2f*m2f;
	m3f = fabsf(m3f);
	m3f2 = m3f*m3f;

	sf_floatwrite(&m1f2,1,ma);
	sf_floatwrite(&a1,1,ma);

	sf_floatwrite(&m2f2,1,ma);
	sf_floatwrite(&a2,1,ma);
	
	sf_floatwrite(&m3f2,1,ma);
	sf_floatwrite(&a3,1,ma);
       
	for (ia = 0; ia < na; ia++) {
	    f = f0 + ia*df;
	    f2 = f*f;

	    data[ia] = a1*exp(-f2/m1f2)*f2/m1f2 + a2*exp(-f2/m2f2)*f2/m2f2 + a3*exp(-f2/m3f2)*f2/m3f2;
	}
        
	if (verb) sf_warning("m1f=%g m2f=%g a1=%g a2=%g",m1f,m2f,a1*m1f*sqrtf(SF_PI)*0.5,a2*m2f*sqrtf(SF_PI)*0.5);
	if (verb) sf_warning ("%d of %d, %d iterations", i2+1, n2, iter);
        
	sf_floatwrite (data,na,out);
    }

    sf_warning(".");

    exit (0);
}

/*tsai $*/
