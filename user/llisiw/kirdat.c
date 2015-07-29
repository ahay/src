/* Utilities for Kirchhoff redatuming. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

static float dt;
static int nsam, mts;
static float** filt;

void filt_init(float dt0 /* time sampling */,
	       float length /* filter length */)
/*< initialize filter >*/
{
    int mts;

    dt = dt0;
    nsam = (int)(length/dt)+2;

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    filt = sf_floatalloc2(nsam,mts);
}

void filt_close()
/*< close >*/
{
    int its;

    for (its=0; its < mts; its++) free(filt[its]);

    free(filt);
}

void filt_set(float tau /* time delay */)
/*< set up filter >*/
{
    int its, isam;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    /* value */
    filt[its][0] = 0.;

    for (isam=1; isam < nsam; isam++) {
	filt[its][isam] = sqrtf(powf((tau+isam*dt)/tau,2.)-1.);
    }

    /* first derivative */
    for (isam=0; isam < nsam-1; isam++) {
	filt[its][isam] = filt[its][isam+1]-filt[its][isam];
    }

    /* second derivative */
    for (isam=nsam-2; isam > 0; isam--) {
	filt[its][isam] = filt[its][isam]-filt[its][isam-1];
    }
}

float pick(float delta /* sample position */,
	   float* trace /* input trace */,
	   int shift /* sample shift */)
/*< filter input trace for one output sample >*/
{
    float value=0.;
    int its, isam;    

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    for (isam=0; (isam < nsam-1) && (shift-isam >= 0); isam++) {
	value += ((1.-delta)*trace[shift-isam]+delta*trace[shift-isam+1])
	    *filt[its][isam]/dt;
    }

    return value;
}
