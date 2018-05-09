/* N-dimensional smooth divisioni with OPENMP parallelization */
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
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "conjgradp.h"
#include "weightp.h"
#include <time.h>
#include "trianglenp.h"
/*^*/

static int niter, n;
static float *p;

void sf_divnp_init(int ndim   /* number of dimensions */, 
		  int nd     /* data size */, 
		  int *ndat  /* data dimensions [ndim] */, 
		  int *nbox  /* smoothing radius [ndim] */, 
		  int niter1 /* number of iterations */,
		  bool verb  /* verbosity */) 
/*< initialize >*/
{
    niter = niter1;
    n = nd;

    sf_trianglenp_init(ndim, nbox, ndat);
    sf_conjgradp_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nd);
}

void sf_divnp_close (void)
/*< free allocated storage >*/
{
    sf_trianglenp_close();
    sf_conjgradp_close();
    free (p);
}

void sf_divnp (float* num, float* den,  float* rat)
/*< smoothly divide rat=num/den >*/
{
    sf_weightp_init(den);
    sf_conjgradp( NULL, sf_weightp_lop,sf_trianglenp_lop,p,rat,num,niter); 
}

void sf_divnep (float* num, float* den,  float* rat, float eps)
/*< smoothly divide rat=num/den with preconditioning >*/
{
    int i;
    double norm;
    double t_start, t_end, t_count;
    double td_start, td_end, td_count;
    if (eps > 0.0f) {

#ifdef _OPENMP
#pragma omp parallel for private(norm)
#endif
	for (i=0; i < n; i++) {
	    norm = 1.0/hypot(den[i],eps);
	    num[i] *= norm;
	    den[i] *= norm;
	}
    } 
#ifdef _OPENMP
      td_start = omp_get_wtime();
#endif

      norm =  pblas_dsdot(n,den,1,den,1);

#ifdef _OPENMP
      td_end = omp_get_wtime();
      td_count = (td_end - td_start);
//      sf_warning("\nSingle Dot Product Time (sec) : %g\n",td_count);
#endif

    if (norm == 0.0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i < n; i++) {
	    rat[i] = 0.0;
	}
	return;
    }
    norm = sqrt(n/norm);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < n; i++) {
	num[i] *= norm;
	den[i] *= norm;
    }   

    sf_weightp_init(den);
#ifdef _OPENMP
    t_start = omp_get_wtime();
#endif
    sf_conjgradp(NULL, sf_weightp_lop,sf_trianglenp_lop,p,rat,num,niter); 
#ifdef _OPENMP
    t_end = omp_get_wtime();
    t_count = (t_end-t_start);
sf_warning("\nSingle dot product time (sec) : %g\n",td_count);
sf_warning("\nConj Grad Elapsed Time (sec): %g\n",t_count);
#endif
}


