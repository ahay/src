/* 3-D  Fourier transform */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static int nkx,nky,nkz;
static kiss_fft_cfg cfgx,cfgxi,cfgy,cfgyi,cfgz,cfgzi;
static kiss_fft_cpx *ctracex,*ctracey,*ctracez;

int fft_nk(int n, int opt)
/*< wavenumber for FT >*/
{
	int nk;

    nk = opt? kiss_fft_next_fast_size(n): n;
	return nk;
}

void fft2d_init(int nz, int nx, int opt)
/*< 2D fft init >*/
{
    nkx = opt? kiss_fft_next_fast_size(nx): nx;
	nkz = opt? kiss_fft_next_fast_size(nz): nz;

    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);

    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);
}

void fft2d_close()
/*< 2D fft close >*/
{
	free(ctracex);
	free(ctracez);
	free(cfgx);
	free(cfgxi);
	free(cfgz);
	free(cfgzi);
}

void fft2d_JYAN(bool inv,       /* forward or inverse */ 
	    int nz, int nx, /* dimensions */
	    kiss_fft_cpx **data,    /* data [nkx][nkz] */
		float **datrl /* if inv==1 * datrl[nx][nz] */)
/*< 2-D transform >*/
{
    int ix, iz, ikx, ikz;

    if (inv) {
		/* Inverse FFT*/
		for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kx to x */
             kiss_fft_stride(cfgxi,(kiss_fft_cpx *)data[0]+ikz,(kiss_fft_cpx *)ctracex,nkz); 
             for (ikx=0; ikx < nkx; ikx++) data[ikx][ikz] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0); 
        }
        for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
            kiss_fft_stride(cfgzi,(kiss_fft_cpx *)data[ikx],(kiss_fft_cpx *)ctracez,1); 
            for (ikz=0; ikz < nkz; ikz++) data[ikx][ikz] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0); 
        }

        for (ix=0; ix < nx; ix++){
            for (iz=0; iz < nz; iz++){ 
                datrl[ix][iz] = sf_crealf(data[ix][iz]); 
                datrl[ix][iz] /= (nkx*nkz); 
            }
		}
	} else {
	/* compute u(kx,kz) */
        for (ix=0; ix < nx; ix++){
            /* Fourier transform z to kz */
            for (iz=1; iz < nz; iz+=2){
                data[ix][iz] = sf_cneg(data[ix][iz]);
            }
            kiss_fft_stride(cfgz,data[ix],ctracez,1); 
            for (ikz=0; ikz<nkz; ikz++) data[ix][ikz] = ctracez[ikz]; 
        }

        for (ikz=0; ikz < nkz; ikz++){
             /* Fourier transform x to kx */
            for (ikx=1; ikx<nkx; ikx+=2){
                data[ikx][ikz] = sf_cneg(data[ikx][ikz]);
            }
            kiss_fft_stride(cfgx,data[0]+ikz,ctracex,nkz); 
            for (ikx=0; ikx<nkx; ikx++) data[ikx][ikz] = ctracex[ikx];
            }
	}
}

void fft3d_init(int nz, int nx, int ny, int opt)
/*< 3D fft init >*/
{
    nkx = opt? kiss_fft_next_fast_size(nx): nx;
    nky = opt? kiss_fft_next_fast_size(ny): ny;
	nkz = opt? kiss_fft_next_fast_size(nz): nz;

    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgy = kiss_fft_alloc(nky,0,NULL,NULL);
    cfgyi = kiss_fft_alloc(nky,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);

    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracey = (kiss_fft_cpx *) sf_complexalloc(nky);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);
}

void fft3d_close()
/*< 3D fft close >*/
{
	free(ctracex);
	free(ctracey);
	free(ctracez);
	free(cfgx);
	free(cfgxi);
	free(cfgy);
	free(cfgyi);
	free(cfgz);
	free(cfgzi);
}

#ifdef _sgsdfgsh
void fft3d_JYAN(bool inv,       /* forward or inverse */ 
	    int nz, int nx, int ny,  /* dimensions */
	    kiss_fft_cpx ***data,    /* data [nky][nkx][nkz] */
		float ***datrl /* if inv==1 * datrl[ny][nx][nz] */)
/*< 3-D transform >*/
{
    int ix, iy, iz, ikx, iky, ikz;

    if (inv) {
		/* Inverse FFT*/
		for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kx to x */
             kiss_fft_stride(cfgxi,(kiss_fft_cpx *)data[0]+ikz,(kiss_fft_cpx *)ctracex,nkz); 
             for (ikx=0; ikx < nkx; ikx++) data[ikx][ikz] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0); 
        }
        for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
            kiss_fft_stride(cfgzi,(kiss_fft_cpx *)data[ikx],(kiss_fft_cpx *)ctracez,1); 
            for (ikz=0; ikz < nkz; ikz++) data[ikx][ikz] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0); 
        }

        for (ix=0; ix < nx; ix++){
            for (iz=0; iz < nz; iz++){ 
                datrl[ix][iz] = sf_crealf(data[ix][iz]); 
                datrl[ix][iz] /= (nkx*nkz); 
            }
		}
	} else {
	/* compute u(kx,kz) */
		for (iy=0; iy < ny; iy++) {
        for (ix=0; ix < nx; ix++) {
            /* Fourier transform z to kz */
            for (iz=1; iz < nz; iz+=2){
                data[iy][ix][iz] = sf_cneg(data[iy][ix][iz]);
            }
            kiss_fft_stride(cfgz,data[iy][ix],ctracez,1); 
            for (ikz=0; ikz<nkz; ikz++) data[iy][ix][ikz] = ctracez[ikz]; 
        }

        for (ikz=0; ikz < nkz; ikz++){
             /* Fourier transform x to kx */
            for (ikx=1; ikx<nkx; ikx+=2){
                data[ikx][ikz] = sf_cneg(data[ikx][ikz]);
            }
            kiss_fft_stride(cfgx,data[0]+ikz,ctracex,nkz); 
            for (ikx=0; ikx<nkx; ikx++) data[ikx][ikz] = ctracex[ikx];
            }
	}
}

#endif
