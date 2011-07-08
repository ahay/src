/* acoustic modeling */

/*
  Copyright (C) 2011 KAUST
  
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
#include "wave.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _wave_h

typedef struct wave2d *wave2dp;
/*^*/

typedef struct abcone2d *abcone2dp;
/*^*/

typedef struct sponge2d *sponge2dp;
/*^*/

typedef struct fftdiff2d *fftdiff2dp;
/*^*/

struct wave2d{
    sf_axis ax,az; /* original */
    sf_axis px,pz; /* expanded */
    sf_axis at;
    int nbx,nbz;   /* num of boundary */
    bool fsrf;     /* free surface */
    /* iz at surface = W->fsrf ? 0 : W->nbz */
};
/*^*/

struct abcone2d{
    float *bxl, *bxh; /* [nz] */
    float *bzl, *bzh; /* [nx] */
};
/*^*/

struct sponge2d{
    float *wx; /* [nx] */
    float *wz; /* [nz] */
};
/*^*/

struct fftdiff2d{
    int Nx,Nz; /* fft size */
    int nkx,nkz;
    float *kx; /* [nkx] = 0 .. 1/2 */
    float *kz; /* [nkz] */
    kiss_fftr_cfg cfgx,icfgx,cfgxc,icfgxc;
    kiss_fftr_cfg cfgz,icfgz,cfgzc,icfgzc;
    float *wx; /* [nkx] low-pass filter */
    float *wz; /* [nkz] */
    sf_complex *mx;  /* [nkx] first deriv */
    sf_complex *mz;  /* [nkz] */
    sf_complex *mmx; /* [nkx] second deriv */
    sf_complex *mmz; /* [nkz] */
    int ompchunk;
    sf_complex **ccx; /* [nz][nkx] temporary */
    sf_complex **ccz; /* [nx][nkz] */
    float      **rrz; /* [nx][Nz]  */
};
/*^*/

#endif

wave2dp wave2d_init(sf_axis x, sf_axis z, sf_axis t,
		    int nbx, int nbz, bool fsrf)
/*< initialize grid >*/
{
    int nx,nz,nt,npx,npz;
    float ox,oz,ot,dx,dz,dt,opx,opz;
    wave2dp W;
    
    W = (wave2dp)sf_alloc(1,sizeof(*W));
    nx = sf_n(x); ox = sf_o(x); dx = sf_d(x);
    nz = sf_n(z); oz = sf_o(z); dz = sf_d(z);
    nt = sf_n(t); ot = sf_o(t); dt = sf_d(t);

    npx = nx + 2*nbx;
    opx = ox - nbx*dx;
    npz = nz + (fsrf ? 1 : 2)*nbz;
    opz = oz - (fsrf ? 0 : 1)*nbz*dz;

    W->ax = sf_maxa(nx,ox,dx);
    W->az = sf_maxa(nz,oz,dz);
    W->at = sf_maxa(nt,ot,dt);
    W->px = sf_maxa(npx,opx,dx);
    W->pz = sf_maxa(npz,opz,dz);

    W->fsrf = fsrf;
    W->nbx = nbx;
    W->nbz = nbz;

    return W;
}

void wave2d_info(wave2dp W)
/*< debug >*/
{
    sf_raxa(W->ax);
    sf_raxa(W->az);
    sf_raxa(W->at);
    sf_raxa(W->px);
    sf_raxa(W->pz);
    sf_warning("nbx=%d,nbz=%d,fsrf=%d",W->nbx,W->nbz,W->fsrf);
}

void expand2(float **a,      /* input */
	     float **b,      /* output */
	     wave2dp W)
/*< expand 2D >*/
{
    int ix,iz,nx,nz,nbx,nbz,npx,npz,srf;
    nx = sf_n(W->ax);
    nz = sf_n(W->az);
    npx= sf_n(W->px);
    npz= sf_n(W->pz);
    srf = W->fsrf ? 0 : W->nbz;
    nbx = W->nbx;
    nbz = W->nbz;

    for     (iz=0; iz < nz; iz++) {
	for (ix=0; ix < nx; ix++) {
	    b[ix+nbx][iz+srf] = a[ix][iz];
	}
    }
    for     (iz=0; iz < nz; iz++) {
	for (ix=0; ix < nbx; ix++) {
	    b[ix      ][iz+srf]=a[0][iz];
	    b[npx-1-ix][iz+srf]=a[nx-1][iz];
	}
    }
    for     (ix=0; ix < npx; ix++) {
	for (iz=0; iz < nbz; iz++) {
	    if (!W->fsrf)
		b[ix][iz] = b[ix][nbz];
	    b[ix][npz-1-iz] = b[ix][npz-nbz-1];
	}
    }
}

void cut2(float **a,      /* input */
	  float **b,      /* output */
	  wave2dp W)
/*< cut 2D >*/
{
    int ix,iz,nz,nx,srf,nbx;
    nx = sf_n(W->ax);
    nz = sf_n(W->az);
    srf= W->fsrf ? 0 : W->nbz;
    nbx= W->nbx;

    for     (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    b[ix][iz] = a[ix+nbx][iz+srf];
	}
    }
}

abcone2dp abcone2d_init(float **vel, /* [npz][npx] */
			wave2dp W)
/*< initialize absorbing BC >*/
{
    abcone2dp abc;
    float nu,dt,dx,dz;
    int ix,iz,x,z,nx,nz,srf,nbx,nbz,npx,npz;

    abc=(abcone2dp)sf_alloc(1,sizeof(*abc));

    nx = sf_n(W->ax); npx = sf_n(W->px); nbx= W->nbx;
    nz = sf_n(W->az); npz = sf_n(W->pz); nbz= W->nbz;
    dx = sf_d(W->ax);
    dz = sf_d(W->az);
    dt = sf_d(W->at);
    srf = W->fsrf ? 0 : W->nbz;

    abc->bxl = sf_floatalloc(nz);
    abc->bxh = sf_floatalloc(nz);
    abc->bzl = sf_floatalloc(nx);
    abc->bzh = sf_floatalloc(nx);

    for (iz=0; iz < nz; iz++) {
	z = iz + srf;	

	x = nbx;       /* left */
	nu= vel[x][z]*dt/dx;
	abc->bxl[iz] = (1.0f-nu)/(1.0f+nu);

	x = npx-nbx-1; /* right */
	nu= vel[x][z]*dt/dx;
	abc->bxh[iz] = (1.0f-nu)/(1.0f+nu);
    }
    for (ix=0; ix < nx; ix++) {
	x = ix + nbx; 

	if (W->fsrf) 
	    abc->bzl[ix] = 0.0f;
	else {
	    z = nbz;   /* top */
	    nu= vel[x][z]*dt/dz;
	    abc->bzl[ix] = (1.0f-nu)/(1.0f+nu);
	}

	z = npz-nbz-1; /* bottom */
	nu= vel[x][z]*dt/dz;
	abc->bzh[ix] = (1.0f-nu)/(1.0f+nu);
    }

    return abc;
}


abcone2dp abcone2d_tau_init(float **vel, /* [npz][npx] */
			    float v0,
			    wave2dp W)
/*< initialize absorbing BC >*/
{
    abcone2dp abc;
    float nu,dt,dx,dz;
    int ix,iz,x,z,nx,nz,srf,nbx,nbz,npx,npz;

    abc=(abcone2dp)sf_alloc(1,sizeof(*abc));

    nx = sf_n(W->ax); npx = sf_n(W->px); nbx= W->nbx;
    nz = sf_n(W->az); npz = sf_n(W->pz); nbz= W->nbz;
    dx = sf_d(W->ax);
    dz = sf_d(W->az);
    dt = sf_d(W->at);
    srf = W->fsrf ? 0 : W->nbz;

    abc->bxl = sf_floatalloc(nz);
    abc->bxh = sf_floatalloc(nz);
    abc->bzl = sf_floatalloc(nx);
    abc->bzh = sf_floatalloc(nx);

    for (iz=0; iz < nz; iz++) {
	z = iz + srf;	

	x = nbx;       /* left */
	nu= vel[x][z]*dt/dx;
	abc->bxl[iz] = (1.0f-nu)/(1.0f+nu);

	x = npx-nbx-1; /* right */
	nu= vel[x][z]*dt/dx;
	abc->bxh[iz] = (1.0f-nu)/(1.0f+nu);
    }
    for (ix=0; ix < nx; ix++) {
	x = ix + nbx; 

	if (W->fsrf) 
	    abc->bzl[ix] = 0.0f;
	else {
	    z = nbz;   /* top */
	    nu= v0*dt/dz;
	    abc->bzl[ix] = (1.0f-nu)/(1.0f+nu);
	}

	z = npz-nbz-1; /* bottom */
	nu= v0*dt/dz;
	abc->bzh[ix] = (1.0f-nu)/(1.0f+nu);
    }

    return abc;
}

void abcone2d_absorb(float **up, /* [npz][npx] = output */
                     float **uo, /* [npz][npx] = input */
                     wave2dp W,
                     abcone2dp abc)
/*< apply absorbing BC >*/
{
    int ix,iz,x,z,nx,nz,nbx,nbz,npx,npz,srf;

    nx = sf_n(W->ax); nbx = W->nbx; npx = sf_n(W->px);
    nz = sf_n(W->az); nbz = W->nbz; npz = sf_n(W->pz);
    srf = W->fsrf ? 0 : W->nbz;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ix,iz,x,z)				\
    shared(up,uo,nx,nz,nbx,nbz,npx,npz,abc,srf)
#endif
    for (iz=0; iz < nz; iz++) {
	z = iz + srf;
	for (ix=0; ix <= nbx; ix++) {
	    x = nbx-ix;        /* left */
	    up[x][z] = uo[x+1][z] + abc->bxl[iz]*(uo[x][z] - up[x+1][z]);
	    x = npx-nbx-1+ix;  /* right */
	    up[x][z] = uo[x-1][z] + abc->bxh[iz]*(uo[x][z] - up[x-1][z]);
	}
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ix,iz,x,z)				\
    shared(up,uo,nx,nz,nbx,nbz,npx,npz,abc,srf,W)
#endif
    for (ix=0; ix < nx; ix++) {
	x = ix + nbx;
	for (iz=0; iz <= nbz; iz++) {
	    if (!W->fsrf) {
		z = nbz-iz;   /* top */
		up[x][z] = uo[x][z+1] + abc->bzl[ix]*(uo[x][z] - up[x][z+1]);
	    }
	    z = npz-nbz-1+iz; /* bottom */
	    up[x][z] = uo[x][z-1] + abc->bzh[ix]*(uo[x][z] - up[x][z-1]); 
	}
    }
} 

sponge2dp sponge2d_init(float eta,
			wave2dp W)
/*< initialze sponge >*/
/* adapted from GEOPHYSICS VOL 50 1985 Cerjan. eta=0.0015^2 is adviced for nb=20 */
{
    sponge2dp spo;
    int ix,iz,nx,nbx,nz,nbz;

    nx = sf_n(W->ax); nbx = W->nbx; 
    nz = sf_n(W->az); nbz = W->nbz; 

    spo=(sponge2dp)sf_alloc(1,sizeof(*spo));
    spo->wx = sf_floatalloc(nbx);
    spo->wz = sf_floatalloc(nbz);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(ix) shared(nbx,spo,eta)
#endif
    for (ix=1; ix <= nbx; ix++)
	spo->wx[ix-1] = exp(-eta*ix*ix);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(iz) shared(nbz,spo,eta)
#endif
    for (iz=1; iz <= nbz; iz++)
	spo->wz[iz-1] = exp(-eta*iz*iz);

    return spo;
}

void sponge2d_sponge(float **u,
		     wave2dp W,
		     sponge2dp spo)
/*< apply sponge >*/
{
    int ix,iz,x,z,nx,nz,nbx,nbz,npx,npz,srf;

    nx = sf_n(W->ax); nbx = W->nbx; npx = sf_n(W->px);
    nz = sf_n(W->az); nbz = W->nbz; npz = sf_n(W->pz);
    srf = W->fsrf ? 0 : W->nbz;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(ix,iz,x,z) shared(nbx,npx,srf,spo,u)
#endif
    for (iz=0; iz < nz; iz++) {
	z = iz + srf;
	for (ix=1; ix <= nbx; ix++) {
	    x = nbx-ix;
	    u[x][z] *= spo->wx[ix];
	    x = npx-nbx-1+ix;
	    u[x][z] *= spo->wx[ix];
	}
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(ix,iz,x,z) shared(nbz,npz,srf,W,spo,u)
#endif
    for (ix=0; ix < nx; ix++) {
	x = ix + nbx;
	for (iz=1; iz <= nbz; iz++) {
	    if (!W->fsrf) {
		z = nbz-iz;
		u[x][z] *= spo->wz[iz];
	    }
	    z = npz-nbz-1+iz;
	    u[x][z] *= spo->wz[iz];
	}
    }
}


void diffx2d(float **a, /* [nz][nx] = input */
	     float **b, /* [nz][nx] = output */
	     int nx, int nz,
	     float dx)
/*< differencing in x >*/
{
    int ix,iz;
    float C1= 0.666666666666666666667; /*  2/3 */
    float C2=-0.083333333333333333333; /* -1/12 */
    float F0=-1.500000000000000000000;
    float F1= 2.000000000000000000000;
    float F2=-0.500000000000000000000;
    float idx=1.0f/dx;
    
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<2; ix++)
	    b[iz][ix]=(F0*a[iz][ix]+F1*a[iz][ix+1]+F2*a[iz][ix+2])*idx;
	for (ix=2; ix<nx-2; ix++)
	    b[iz][ix]=(C1*(a[iz][ix+1]-a[iz][ix-1])
		      +C2*(a[iz][ix+2]-a[iz][ix-2]))*idx;
	for (ix=nx-2; ix<nx; ix++)
	    b[iz][ix]=-(F0*a[iz][ix]+F1*a[iz][ix-1]+F2*a[iz][ix-2])*idx;
    }
}

void diffz2d(float **a, /* [nz][nx] = input */
	     float **b, /* [nz][nx] = output */
	     int nx, int nz,
	     float dz)
/*< differencing in z >*/
{
    int ix,iz;
    float C1= 0.666666666666666666667; /*  2/3 */
    float C2=-0.083333333333333333333; /* -1/12 */
    float F0=-1.500000000000000000000;
    float F1= 2.000000000000000000000;
    float F2=-0.500000000000000000000;
    float idz=1.0f/dz;
    
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<2; iz++)
	    b[iz][ix]=(F0*a[iz][ix]+F1*a[iz+1][ix]+F2*a[iz+2][ix])*idz;
	for (iz=2; iz<nz-2; iz++)
	    b[iz][ix]=(C1*(a[iz+1][ix]-a[iz-1][ix])
		      +C2*(a[iz+2][ix]-a[iz-2][ix]))*idz;
	for (iz=nz-2; iz<nz; iz++)
	    b[iz][ix]=-(F0*a[iz][ix]+F1*a[iz-1][ix]+F2*a[iz-2][ix])*idz;
    }
}
