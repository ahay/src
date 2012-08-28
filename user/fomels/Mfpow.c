/* Time/frequency power estimation */
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

static float newton(float bi, int imin, int imax, const float* spec, const float *w, float *db);

int main(int argc, char* argv[])
{
    bool verb, time;
    int n1, n2, i1, i2, imin, imax, iter, niter, ib, nb, i;
    float fmin, fmax, o1, d1, *spec, bi, *b, *w, semb, db, deltab; 
    float omega, rho, tol, smax, b0, bmin, bmax, oldsemb;
    sf_file inp, out, beta;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    beta = sf_output("beta");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"o1",&o1)) o1=0.0; 
    if (!sf_histfloat(inp,"d1",&d1)) d1=1.0; 
    n2 = sf_leftsize(inp,1);

    spec = sf_floatalloc(n1);
    w = sf_floatalloc(n1);

    sf_unshiftdim(inp,beta,1);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of Newton iterations */
    if (!sf_getfloat("fmin",&fmin)) fmin=o1;
    /* minimum frequency */
    if (!sf_getfloat("fmax",&fmax)) fmax=o1+(n1-1)*d1;
    /* maximum frequency */
    if (!sf_getfloat("bmin",&bmin)) bmin=-1.0;
    /* minimum value of beta */
    if (!sf_getfloat("bmax",&bmax)) bmax=1.0;
    /* maximum value of beta */
    if (!sf_getint("nb",&nb)) nb=10;
    db = (bmax-bmin)/(nb-1);

    if (!sf_getfloat("tol",&tol)) tol=SF_EPS;
    /* accuracy tolerance for beta */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getbool("time",&time)) time=false;
    /* time axis */

    imin = SF_MAX(floorf(0.5+fmin/d1),0);
    imax = SF_MIN(floorf(1.5+fmax/d1),n1);

    if (imax <= imin) sf_error("Wrong fmin or fmax");

    rho = 1.0-1.0/n1;

    for (i1=0; i1 < n1; i1++) {
	if (time) {
	    w[i1] = o1+i1*d1;
	} else {
	    omega = SF_PI*i1/(n1-1);
	    w[i1] = hypotf(1.-rho*cosf(omega),rho*sinf(omega));
	}
    }

    b = sf_floatalloc(n2);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(spec,n1,inp);

	/* rough search for the best beta */

	smax = 0.;
	b0 = bmin;

	for (ib=0; ib < nb; ib++) {
	    bi = bmin + ib*db;

	    semb = newton(bi, imin, imax, spec, w, NULL);
	    if (semb > smax) {
		smax = semb;
		b0 = bi;
	    }
	}

	/* Newton optimization */

	bi=b0;
	deltab=0.;
	oldsemb=0.;

	for (iter=0; iter < niter; iter++) {
	    semb = newton(bi+deltab, imin, imax, spec, w, &deltab);

	    if (semb < oldsemb) {
		/* try decresing or reversing step */
		for (i=0; i < 10 && semb < oldsemb; i++) {
		    if (i%2) {
			deltab = -deltab;
		    } else {
			deltab *= 0.5;
		    }
		    semb = newton(bi+deltab, imin, imax, spec, w, NULL);
		}
	    }

	    bi += deltab;
	    oldsemb = semb;

	    if (verb) sf_warning("beta=%g semb=%g",bi,semb);
	    if (fabsf(deltab) < tol) break;
	}

	b[i2] = bi;
	
	for (i1=0; i1 < n1; i1++) {
	    spec[i1] = powf(w[i1],bi);
	}

	sf_floatwrite(spec,n1,out);
    }

    sf_floatwrite(b,n2,beta);

    exit(0);
}

static float newton(float bi, int imin, int imax, const float* spec, const float *w, float *db) 
/* Newton step */
{
    int i1;
    float s, s2, sl, s2l, sl2, s2l2, f, l, semb;

    s=s2=sl=s2l=sl2=s2l2=0.0;
    
    for (i1=imin; i1 < imax; i1++) {
	f = spec[i1]*powf(w[i1],bi);
	l = logf(w[i1]);
	s += f;
	s2 += f*f;
	sl += f*l;
	s2l += f*f*l;
	sl2 += f*l*l;
	s2l2 += f*f*l*l;
    }

    semb = s*s/(s2*(imax-imin));

    if (NULL != db) 
	*db = (s*s2*(s*s2l - s2*sl))/
	    (s*s*(4*s2l*s2l - 2*s2*s2l2) + s2*s2*sl*sl + s*s2*(s2*sl2-4*s2l*sl));

    return semb;
}
