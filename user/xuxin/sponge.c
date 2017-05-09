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
#include <rsf.h>

typedef struct {
    int bzl; float *wzl;
    int bzh; float *wzh;
    int bxl; float *wxl;
    int bxh; float *wxh;
} Sponge;
/*^*/

Sponge *sponge_init(const int *b,   /* {bzl,bzh,bxl,bxh} */
					const float *c  /* {czl,czh,cxl,cxh} */)
					
/*< >*/
{
    int i;
    Sponge *S;

    int bzl = b[0]; float czl = c[0];
    int bzh = b[1]; float czh = c[1];
    int bxl = b[2]; float cxl = c[2];
    int bxh = b[3]; float cxh = c[3];

    float *wzl = (float *)calloc(bzl,sizeof(float));
    float *wzh = (float *)calloc(bzh,sizeof(float));
    float *wxl = (float *)calloc(bxl,sizeof(float));
    float *wxh = (float *)calloc(bxh,sizeof(float));

    for (i=0; i < bzl; i++) wzl[i] = expf(-czl * i*i);
    for (i=0; i < bzh; i++) wzh[i] = expf(-czh * i*i);
    for (i=0; i < bxl; i++) wxl[i] = expf(-cxl * i*i);
    for (i=0; i < bxh; i++) wxh[i] = expf(-cxh * i*i);

    S = malloc(sizeof(Sponge));
    S->bzl = bzl; S->wzl = wzl;
    S->bzh = bzh; S->wzh = wzh;
    S->bxl = bxl; S->wxl = wxl;
    S->bxh = bxh; S->wxh = wxh;
    return S;
}

void sponge_apply(float *u,     /* [Nz][Nx] */
                  const int *n, /* {nz,nx} */
                  const Sponge *S)
/*< >*/
{
    int i,iz,ix;
    
    int nz = n[0];
    int nx = n[1];
    int bzl = S->bzl; float *wzl = S->wzl;
    int bzh = S->bzh; float *wzh = S->wzh;
    int bxl = S->bxl; float *wxl = S->wxl;
    int bxh = S->bxh; float *wxh = S->wxh;
    int Nz = nz + bzl + bzh;
    int Nx = nx + bxl + bxh;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(ix,i)
#endif
    for (ix=0; ix < Nx; ix++) {
        for (i=1; i <= bzl; i++)
            u[ix*Nz + bzl-i] *= wzl[i-1];
        for (i=1; i <= bzh; i++)
            u[ix*Nz + bzl+i+nz-1] *= wzh[i-1];
		}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,i)
#endif
    for (iz=0; iz < Nz; iz++) {
        for (i=1; i <= bxl; i++)
            u[(bxl-i)*Nz + iz] *= wxl[i-1];
        for (i=1; i <= bxh; i++)
            u[(bxl+i+nx-1)*Nz + iz] *= wxh[i-1];
    }
}

void sponge_free(Sponge *S)
/*< >*/
{
    free(S->wzl); free(S->wzh);
    free(S->wxl); free(S->wxh);
    free(S);
}
