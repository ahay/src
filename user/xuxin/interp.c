/*
  Copyright (C) 2012 KAUST

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

#include <stdlib.h>

typedef struct {
    int nz; float z0,dz,*z; /* mesh */
    int nx; float x0,dx,*x;
} Int2;
/*^*/

Int2 *int2_init(const int   *n, /* {nz,nx} */
                const float *o, /* {z0,x0} */
                const float *d  /* {dz,dx} */)
/*< >*/
{
    int i,nz,nx;
    float *z,*x,z0,x0,dz,dx;
    Int2 *I2;

    nz = n[0]; nx = n[1];
    z0 = o[0]; x0 = o[1];
    dz = d[0]; dx = d[1];

    z = (float *)calloc(nz,sizeof(float));
    x = (float *)calloc(nx,sizeof(float));

    for (i=0; i < nz; i++) z[i] = z0 + i * dz;
    for (i=0; i < nx; i++) x[i] = x0 + i * dx;

    I2 = malloc(sizeof(Int2));
    I2->nz = nz; I2->nx = nx;
    I2->z0 = z0; I2->x0 = x0;
    I2->dz = dz; I2->dx = dx;
    I2->z  = z;  I2->x  = x;
    return I2;
}

void int2_apply(float *F,       /* [m] */
                const float *f, /* [nx][nz] */
                int m,
                const float *c, /* [m][2] coord */
                const Int2 *I2)
/*< evaluate f at c, return as F. >*/
{
	int   i,iz,ix;
	float zo,xo,wx,wz;

    int   nz = I2->nz; int   nx = I2->nx;
    float dz = I2->dz; float dx = I2->dx;
    float z0 = I2->z0; float x0 = I2->x0;
    float *z = I2->z;  float *x = I2->x;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i,zo,xo,iz,ix,wz,wx)
#endif
    for (i=0; i < m; i++) {
        zo = c[i*2  ]; iz = (int)((zo - z0) / dz);
        xo = c[i*2+1]; ix = (int)((xo - x0) / dx);

        iz = (iz < 0   ) ? 0    : iz;
        iz = (iz > nz-2) ? nz-2 : iz;
        ix = (ix < 0   ) ? 0    : ix;
        ix = (ix > nx-2) ? nx-2 : ix;

        wx = (xo - x[ix]) / dx;
        wz = (zo - z[iz]) / dz;
	
        F[i] = 
            (1.-wx)*(1.-wz) * f[ ix   *nz+ iz   ] +			\
            wx     *(1.-wz) * f[(ix+1)*nz+ iz   ] +			\
            (1.-wx)*    wz  * f[ ix   *nz+(iz+1)] +			\
            wx     *    wz  * f[(ix+1)*nz+(iz+1)];
    }
}

void int2_inject(const float *F, /* [m] */
                 float *f,       /* [nx][nz] */
                 int m,
                 const float *c, /* [m][2] coord */
                 const Int2 *I2)
/*< update f using value F at c >*/
{
	int   i,iz,ix;
	float zo,xo,wx,wz;

    int   nz = I2->nz; int   nx = I2->nx;
    float dz = I2->dz; float dx = I2->dx;
    float z0 = I2->z0; float x0 = I2->x0;
    float *z = I2->z;  float *x = I2->x;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i,zo,xo,iz,ix,wz,wx)
#endif
    for (i=0; i < m; i++) {
        zo = c[i*2  ]; iz = (int)((zo - z0) / dz);
        xo = c[i*2+1]; ix = (int)((xo - x0) / dx);

        iz = (iz < 0   ) ? 0    : iz;
        iz = (iz > nz-2) ? nz-2 : iz;
        ix = (ix < 0   ) ? 0    : ix;
        ix = (ix > nx-2) ? nx-2 : ix;

        wx = (xo - x[ix]) / dx;
        wz = (zo - z[iz]) / dz;

        f[ ix   *nz+ iz   ] += F[i] * (1.-wx)*(1.-wz);
        f[(ix+1)*nz+ iz   ] += F[i] * wx     *(1.-wz);
        f[ ix   *nz+(iz+1)] += F[i] * (1.-wx)*    wz ;
        f[(ix+1)*nz+(iz+1)] += F[i] * wx     *    wz ;
    }
}

void int2_inject2(const float *F, /* [m] */
                  float *f,       /* [nx][nz] */
                  int m,
                  const float *c, /* [m][2] coord */
                  const Int2 *I2)
/*< set f using value F at c >*/
{
	int   i,iz,ix;
	float zo,xo,wx,wz;

    int   nz = I2->nz; int   nx = I2->nx;
    float dz = I2->dz; float dx = I2->dx;
    float z0 = I2->z0; float x0 = I2->x0;
    float *z = I2->z;  float *x = I2->x;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i,zo,xo,iz,ix,wz,wx)
#endif
    for (i=0; i < m; i++) {
        zo = c[i*2  ]; iz = (int)((zo - z0) / dz);
        xo = c[i*2+1]; ix = (int)((xo - x0) / dx);

        iz = (iz < 0   ) ? 0    : iz;
        iz = (iz > nz-2) ? nz-2 : iz;
        ix = (ix < 0   ) ? 0    : ix;
        ix = (ix > nx-2) ? nx-2 : ix;

        wx = (xo - x[ix]) / dx;
        wz = (zo - z[iz]) / dz;

        f[ ix   *nz+ iz   ] = F[i] * (1.-wx)*(1.-wz);
        f[(ix+1)*nz+ iz   ] = F[i] * wx     *(1.-wz);
        f[ ix   *nz+(iz+1)] = F[i] * (1.-wx)*    wz ;
        f[(ix+1)*nz+(iz+1)] = F[i] * wx     *    wz ;
    }
}

void int2_free(Int2 *I2)
/*< >*/ 
{
    free(I2->z); free(I2->x); free(I2);
}
