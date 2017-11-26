/* Complex-valued compensation (between two wavefields) */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

void stable_div(int n, float eps, float *num, float *den, float *ratio)
/*< stable division of real numbers >*/
{
    int i;
    float r,vnum,vden;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,vnum,vden,r)
#endif
    for (i = 0; i < n; i++) {
        vnum = num[i];
        vden = den[i];

        if (vden == 0.f || fabsf(vden) >= fabsf(vnum))
            r = 1.f;
        else
            r = vnum*vden/(vden*vden + eps);

        ratio[i] = r;
    }
}

void stable_cdiv(int n, float eps, sf_complex *num, sf_complex *den, sf_complex *ratio)
/*< stable division of complex numbers >*/
{
    int i;
    sf_complex r,vnum,vden;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,vnum,vden,r)
#endif
    for (i = 0; i < n; i++) {
        vnum = num[i];
        vden = den[i];

        if (cabsf(vden) == 0.f || cabsf(vden) >= cabsf(vnum))
            r = sf_cmplx(1.f,0.f);
        else {
#ifdef SF_HAS_COMPLEX_H
            r = vnum*conjf(vden)/(vden*conjf(vden) + eps);
#else
            r = sf_cdiv(sf_cmul(vnum,sf_conjf(vden)),sf_cadd(sf_cmul(vden,sf_conjf(vden)) + sf_cmplx(eps,0.f)));
#endif
        }

        ratio[i] = r;
    }
}

void stable_cdiv_f(int n, float eps, sf_complex *num, sf_complex *den, float *ratio)
/*< stable division of complex numbers using envelope, therefore outputs real numbers >*/
{
    int i;
    float r,vnum,vden;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,vnum,vden,r)
#endif
    for (i = 0; i < n; i++) {
        vnum = cabsf(num[i]);
        vden = cabsf(den[i]);

        if (vden == 0.f || vden >= vnum)
            r = 1.f;
        else
            r = vnum*vden/(vden*vden + eps);

        ratio[i] = r;
    }
}

float find_max(int n, float *vals)
/*< find maximum abs val >*/
{
    int i;
    float v,maxval;

    maxval = 0.f;
    for (i = 0; i < n; i++) {
        v = fabsf(vals[i]);
        if (v>maxval) maxval=v;
    }

    return maxval;
}

float find_cmax(int n, sf_complex *vals)
/*< find maximum abs val >*/
{
    int i;
    float v,maxval;

    maxval = 0.f;
    for (i = 0; i < n; i++) {
        v = cabsf(vals[i]);
        if (v>maxval) maxval=v;
    }

    return maxval;
}

void invert(int n, float *vals, float *ivals)
/*< calculate stable inversion >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < n; i++) {
        if (vals[i] < 1) ivals[i] = 1.f;
        else ivals[i] = 1.f/vals[i];
    }

}

void cinvert(int n, sf_complex *vals, sf_complex *ivals)
/*< calculate stable inversion >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < n; i++) {
        if (cabsf(vals[i]) < 1) 
            ivals[i] = sf_cmplx(1.f,0.f);
        else 
#ifdef SF_HAS_COMPLEX_H
            ivals[i] = 1/vals[i];
#else
            ivals[i] = sf_cdiv(sf_cmplx(1.f,0.f),vals[i]);
#endif
    }
}
