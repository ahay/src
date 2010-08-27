/* Wilson-Burg spectral factorization */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>
#include <float.h>

#include <rsf.h>
/*^*/

#include "wilson.h"

static int n, n2;
static float *au, *bb, *cc, *b, *c;    

void wilson_init( int nmax /* maximum data size */) 
/*< initialize >*/
{
    n  = nmax;
    n2 = 2*n-1;
    au = sf_floatalloc (n2);
    bb = sf_floatalloc (n2);
    cc = sf_floatalloc (n2);
    b  = sf_floatalloc (n);
    c  = sf_floatalloc (n);
}

float wilson_factor(int niter    /* number of iterations */, 
		    float s0     /* zero-lag auto-correlation */, 
		    sf_filter ss /* input auto-correlation */, 
		    sf_filter aa /* output factor */, 
		    bool verb    /* verbosity flag */,
		    float tol    /* tolerance */)
/*< Factor >*/ 
{
    float eps;
    int i, iter;
    
    for(i=0; i < n2; i++) au[i] = 0.;
    au[n-1] = s0;
    b[0] = 1.;       /* initialize */
    
    for(i=0; i < ss->nh; i++) {  /* symmetrize input auto */
	au[n-1+ss->lag[i]] =  ss->flt[i];        
	au[n-1-ss->lag[i]] =  ss->flt[i];
    }                   
    
    sf_helicon_init( aa);           /* multiply polynoms */
    sf_polydiv_init( n2, aa);       /* divide   polynoms */
    for(iter=0; iter < niter; iter++) {
	sf_polydiv_lop(false,false, n2, n2, au, bb);  /* S/A */
	sf_polydiv_lop(true, false, n2, n2, cc, bb);  /* S/(AA') */
	eps = 0.;
	for(i=1; i < n; i++) {   /* b = plusside(1+cc) */
	    b[i] = 0.5*(cc[n-1+i] + cc[n-1-i]) / cc[n-1]; 
	    if (fabs(b[i]) > eps) eps = fabs(b[i]);
	}
	if (verb) sf_warning("wilson %d %f",iter,eps);
	if (eps < tol) break;
	
	sf_helicon_lop( false, false, n, n, b, c);   /* c = A b */
	
        /* put on helix */
	for(i=0; i < aa->nh; i++) aa->flt[i] = c[aa->lag[i]];              
    } 
    return sqrtf(cc[n-1]);
}

void wilson_close( void) 
/*< free allocated storage >*/
{
    sf_polydiv_close();
    free (au);
    free (bb);
    free (cc);
    free (b);
    free (c);
}
