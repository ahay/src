/* Jacobi-like method for general complex matrix */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

static sf_complex *ak;
static bool verb;

void jacobi2_init(int n /* matrix size */,
		  bool verb1 /* verbosity */)
/*< initialize >*/
{
    ak = sf_complexalloc(n);
    verb = verb1;
}

void jacobi2_close(void)
/*< free allocated storage >*/
{
    free(ak);
}

float jacobi2(sf_complex** a /* matrix to rotate */, 
	      int n, int k, int l   /* indeces */) 
/*< Jacobi-like transformation >*/
/* The algorithm is modified from 

   J.P. Killingbeck, A.Grosjean, and G. Jolicard, 2004, A simple
   method for complex eigenvalues: Journal of Physics A: Mathematical
   and General, v. 37, L567-L572. */
{
    int i, j1,j2;
    sf_complex t, s, akl;

    if (k==l || cabsf(a[k][l]) < FLT_EPSILON) {
	if (verb) sf_warning("k==l or abs(a(k,l) = 0");
	return 0.0f;
    }

#ifdef SF_HAS_COMPLEX_H
    t = 0.5*(a[l][l]-a[k][k]);
    s = csqrtf(t*t+a[l][k]*a[k][l]);
#else
    t = sf_crmul(sf_cadd(a[l][l],sf_cneg(a[k][k])),0.5);
    s = csqrtf(sf_cadd(sf_cmul(t,t),sf_cmul(a[l][k],a[k][l])));
#endif

    if (cabsf(t) < FLT_EPSILON && cabsf(s) < FLT_EPSILON) {
	if (verb) sf_warning("abs(t) and abs(s) = 0");
	return 0.0f;
    }
	
#ifdef SF_HAS_COMPLEX_H
    if (cabsf(t+s) > cabsf(t-s)) {
	t = a[k][l]/(t+s);
    } else {
	t = a[k][l]/(t-s);
    }
#else
    if (cabsf(sf_cadd(t,s)) > cabsf(sf_cadd(t,sf_cneg(s)))) {
	t = sf_cdiv(a[k][l],sf_cadd(t,s));
    } else {
	t = sf_cdiv(a[k][l],sf_cadd(t,sf_cneg(s)));
    }
#endif

#ifdef SF_HAS_COMPLEX_H
    akl = a[k][l] + t*(a[k][k]-a[l][l])-t*t*a[l][k];
#else
    akl = sf_cadd(a[k][l],
		  sf_cadd(sf_cmul(t,
				  sf_cadd(a[k][k],sf_cneg(a[l][l]))),
			  sf_cneg(sf_cmul(t,sf_cmul(t,a[l][k])))));
#endif
    for (i=0; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H
	ak[i] = a[i][l] + t*a[i][k];
#else
	ak[i] = sf_cadd(a[i][l],sf_cmul(t,a[i][k]));
#endif
    }
    for (i=0; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H
	a[k][i] -= t*a[l][i];
#else
	a[k][i] = sf_cadd(a[k][i],sf_cmul(t,a[l][i]));
#endif
    }
    for (i=0; i < n; i++) {
	a[i][l] = ak[i];
    }
    a[k][l] = akl;
	
    if (verb) {
	sf_warning("t=%f+i%f (after),   |t|=%f,   akl=%f + i%f\n",
		   crealf(t),cimagf(t), cabsf(t),crealf(akl),cimagf(akl)); 
	for (j1=0; j1<n; j1++) {
	    for (j2=0; j2<n; j2++) {
		if (cabsf(a[j2][j1])>FLT_EPSILON)
		    sf_warning("a[%d,%d] = %f + i%f",j1,j2,
			       crealf(a[j2][j1]),cimagf(a[j2][j1]));
	    }
	}
    }

    return cabsf(t);
}
