/* Modeling of pure acoustic wave in 3-D transversely isotropic meida using optimized pseudo-Laplacian operator */
/*
  Copyright (C)2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
					2009 University of Texas at Austin
  
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
#include "opostep3.h"
#include "fft3.h"
//#include "fft3d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex *cwave, *cwavem;
static float	**Cxxl, **Cxxr, **Cyyl, **Cyyr, **Czzl, **Czzr, **Cxyl, **Cxyr,
		 		  	**Cxzl, **Cxzr, **Cyzl, **Cyzr, **Cxxxxl, **Cxxxxr, **Cyyyyl,
		   		**Cyyyyr, **Czzzzl, **Czzzzr, **Cxxxyl, **Cxxxyr, **Cxxxzl,
		   		**Cxxxzr, **Cxyyyl, **Cxyyyr, **Cyyyzl, **Cyyyzr, **Cxzzzl,
		   		**Cxzzzr, **Cyzzzl, **Cyzzzr, **Cxxyyl, **Cxxyyr, **Cxxzzl,
		   		**Cxxzzr, **Cyyzzl, **Cyyzzr, **Cxxyzl, **Cxxyzr, **Cxyyzl,
		   		**Cxyyzr, **Cxyzzl, **Cxyzzr;

static float ***uuxx, ***uuyy, ***uuzz, ***uuxy, ***uuxz, ***uuyz, 
				 ***uxxxx, ***uyyyy, ***uzzzz, ***uxxxy, ***uxxxz, ***uxyyy, 
				 ***uyyyz, ***uxzzz, ***uyzzz, ***uxxyy, ***uxxzz, ***uyyzz,
				 ***uxxyz, ***uxyyz, ***uxyzz, ***curtmp, **wave;
static int nx, ny, nz, nk, nx2, ny2, nz2, m2, nxyzb2, opt;


void lowrank_init3(int nzb, int nxb, int nyb, int nkxyz, int nkzz, int nkxx, int nkyy, int m, int nxyzb, float **cxxl, float **cxxr, float **cyyl, float **cyyr, float **czzl, float **czzr, float **cxyl, float **cxyr, float **cxzl, float **cxzr, float **cyzl, float **cyzr, float **cxxxxl, float **cxxxxr, float **cyyyyl, float **cyyyyr, float **czzzzl, float **czzzzr, float **cxxxyl, float **cxxxyr, float **cxxxzl, float **cxxxzr, float **cxyyyl, float **cxyyyr, float **cyyyzl, float **cyyyzr, float **cxzzzl, float **cxzzzr, float **cyzzzl, float **cyzzzr, float **cxxyyl, float **cxxyyr, float **cxxzzl, float **cxxzzr, float **cyyzzl, float **cyyzzr, float **cxxyzl, float **cxxyzr, float **cxyyzl, float **cxyyzr, float **cxyzzl, float **cxyzzr)
/*< 3D lowrank initiazation >*/
{
	int ix, iy, iz, ik;

    nz = nzb;
	ny = nyb;
	nx = nxb;
	nk = nkxyz;

	nz2 = nkzz;
	ny2 = nkyy;
	nx2 = nkxx;

	m2 = m;
	nxyzb2 = nxyzb;

	Cxxl = cxxl;
	Cxxr = cxxr;
	Cyyl = cyyl;
	Cyyr = cyyr;
	Czzl = czzl;
	Czzr = czzr;

	Cxyl = cxyl;
	Cxyr = cxyr;
	Cxzl = cxzl;
	Cxzr = cxzr;
	Cyzl = cyzl;
	Cyzr = cyzr;

	Cxxxxl = cxxxxl;
	Cxxxxr = cxxxxr;
	Cyyyyl = cyyyyl;
	Cyyyyr = cyyyyr;
	Czzzzl = czzzzl;
	Czzzzr = czzzzr;

	Cxxxyl = cxxxyl;
	Cxxxyr = cxxxyr;
	Cxxxzl = cxxxzl;
	Cxxxzr = cxxxzr;
	Cxyyyl = cxyyyl;
	Cxyyyr = cxyyyr;
	Cyyyzl = cyyyzl;
	Cyyyzr = cyyyzr;
	Cxzzzl = cxzzzl;
	Cyzzzr = cyzzzr;
	Cyzzzl = cyzzzl;
	Cxzzzr = cxzzzr;

	Cxxyyl = cxxyyl;
	Cxxyyr = cxxyyr;
	Cxxzzl = cxxzzl;
	Cxxzzr = cxxzzr;
	Cyyzzl = cyyzzl;
	Cyyzzr = cyyzzr;

	Cxxyzl = cxxyzl;
	Cxxyzr = cxxyzr;
	Cxyyzl = cxyyzl;
	Cxyyzr = cxyyzr;
	Cxyzzl = cxyzzl;
	Cxyzzr = cxyzzr;

	uuxx = sf_floatalloc3(nz,nx,ny);
    uuyy = sf_floatalloc3(nz,nx,ny);
    uuzz = sf_floatalloc3(nz,nx,ny);
    uuxy = sf_floatalloc3(nz,nx,ny);
    uuxz = sf_floatalloc3(nz,nx,ny);
    uuyz = sf_floatalloc3(nz,nx,ny);
    curtmp = sf_floatalloc3(nz2,nx2,ny2);

    uxxxx = sf_floatalloc3(nz,nx,ny);
    uyyyy = sf_floatalloc3(nz,nx,ny);
    uzzzz = sf_floatalloc3(nz,nx,ny);
    uxxxy = sf_floatalloc3(nz,nx,ny);
    uxxxz = sf_floatalloc3(nz,nx,ny);
    uxyyy = sf_floatalloc3(nz,nx,ny);
    uyyyz = sf_floatalloc3(nz,nx,ny);
    uxzzz = sf_floatalloc3(nz,nx,ny);
    uyzzz = sf_floatalloc3(nz,nx,ny);
    uxxyy = sf_floatalloc3(nz,nx,ny);
    uxxzz = sf_floatalloc3(nz,nx,ny);
    uyyzz = sf_floatalloc3(nz,nx,ny);
    uxxyz = sf_floatalloc3(nz,nx,ny);
    uxyyz = sf_floatalloc3(nz,nx,ny);
    uxyzz = sf_floatalloc3(nz,nx,ny);
 


	cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nxyzb2,m2);

	for (ik = 0; ik < nk; ik++) {
		cwave[ik] = sf_cmplx(0.0, 0.0);
		cwavem[ik] = sf_cmplx(0.0, 0.0);
	}
	for (ix = 0; ix < m2; ix++) {
	    for (iz = 0; iz < nxyzb2; iz++) {
			wave[ix][iz] = 0.0;
		}
	}

}

void lowrank_close3()
/*< free the work space for opo >*/
{

}



void lowrank_comp3(float ***uu, sf_complex *cwave, sf_complex *cwavem, float **lt, float **rt)
/*< K space mix-matrix approximated by lowrank >*/
{
	int i, j, ix, iy, iz, im, ik;
   
	ifft3_allocate(cwavem);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_crmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    ifft3(wave[im],cwavem);
	}

	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz = 0; iz < nz; iz++) {
		    i = iz+nz *(ix+nx *iy); /* original grid */
		    j = iz+nz2*(ix+nx2*iy); /* padded grid */

		    for (im = 0; im < m2; im++) {
			uu[iy][ix][iz]= lt[im][i]*wave[im][j];
		    }
		}
		}
	}
}


void opostep3(float ***old /*previous step*/,
             float ***cur /*current step*/,
             int nz, int nx, int ny /*model size*/,
			 float dz, float dx, float dy /*model grid */,
             float v0 /*reference vel*/,
             float ***v /*reference vel*/,
             float ***sigma /*reference vel*/,
             float ***delta /*reference vel*/,
             float ***seta /*reference vel*/,
             float ***phi /*reference vel*/,
             float dt /*time step size*/)
/*< Optimized pseudo-Laplacian operator>*/
{
    int ix, ikx, iy, iky, ikz, iz;

    for (iky=0; iky < ny; iky++){ 
    for (ikx=0; ikx < nx; ikx++){ 
        for (ikz=0; ikz < nz; ikz++){ 
			uuxx[iky][ikx][ikz] = 0.;
			uuyy[iky][ikx][ikz] = 0.;
			uuzz[iky][ikx][ikz] = 0.;
			uuxy[iky][ikx][ikz] = 0.;
			uuxz[iky][ikx][ikz] = 0.;
			uuyz[iky][ikx][ikz] = 0.;
			
			uxxxx[iky][ikx][ikz] = 0.;
			uyyyy[iky][ikx][ikz] = 0.;
			uzzzz[iky][ikx][ikz] = 0.;
			uxxxy[iky][ikx][ikz] = 0.;
			uxxxz[iky][ikx][ikz] = 0.;
			uxyyy[iky][ikx][ikz] = 0.;
			uyyyz[iky][ikx][ikz] = 0.;
			uxzzz[iky][ikx][ikz] = 0.;
			uyzzz[iky][ikx][ikz] = 0.;
			uxxyy[iky][ikx][ikz] = 0.;
			uxxzz[iky][ikx][ikz] = 0.;
			uyyzz[iky][ikx][ikz] = 0.;
			uxxyz[iky][ikx][ikz] = 0.;
			uxyyz[iky][ikx][ikz] = 0.;
			uxyzz[iky][ikx][ikz] = 0.;
         }  
	}
	}

    for (iy=0; iy < ny2; iy++){ 
    for (ix=0; ix < nx2; ix++){ 
        for (iz=0; iz < nz2; iz++){ 
		    	curtmp[iy][ix][iz] = 0.;
         }  
	}
	}
	
    for (iy=0; iy < ny; iy++){ 
    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
		    curtmp[iy][ix][iz] = cur[iy][ix][iz];	
         }  
	}
	}

	/* Just do one FFT and severn IFFT in one time step */

	/* Transform the input data into wavenumber space domain */
	fft3(curtmp[0][0], cwave);

	lowrank_comp3(uuxx, cwave, cwavem, Cxxl, Cxxr);
	lowrank_comp3(uuyy, cwave, cwavem, Cyyl, Cyyr);
	lowrank_comp3(uuzz, cwave, cwavem, Czzl, Czzr);
	lowrank_comp3(uuxy, cwave, cwavem, Cxyl, Cxyr);
	lowrank_comp3(uuxz, cwave, cwavem, Cxzl, Cxzr);
	lowrank_comp3(uuyz, cwave, cwavem, Cyzl, Cyzr);
	lowrank_comp3(uxxxx, cwave, cwavem, Cxxxxl, Cxxxxr);
	lowrank_comp3(uyyyy, cwave, cwavem, Cyyyyl, Cyyyyr);
	lowrank_comp3(uzzzz, cwave, cwavem, Czzzzl, Czzzzr);
	
	lowrank_comp3(uxxxy, cwave, cwavem, Cxxxyl, Cxxxyr);
	lowrank_comp3(uxxxz, cwave, cwavem, Cxxxzl, Cxxxzr);
	lowrank_comp3(uxyyy, cwave, cwavem, Cxyyyl, Cxyyyr);
	lowrank_comp3(uyyyz, cwave, cwavem, Cyyyzl, Cyyyzr);
	lowrank_comp3(uxzzz, cwave, cwavem, Cxzzzl, Cxzzzr);
	lowrank_comp3(uyzzz, cwave, cwavem, Cyzzzl, Cyzzzr);

	lowrank_comp3(uxxyy, cwave, cwavem, Cxxyyl, Cxxyyr);
	lowrank_comp3(uxxzz, cwave, cwavem, Cxxzzl, Cxxzzr);
	lowrank_comp3(uyyzz, cwave, cwavem, Cyyzzl, Cyyzzr);

	lowrank_comp3(uxxyz, cwave, cwavem, Cxxyzl, Cxxyzr);
	lowrank_comp3(uxyyz, cwave, cwavem, Cxyyzl, Cxyyzr);
	lowrank_comp3(uxyzz, cwave, cwavem, Cxyzzl, Cxyzzr);

#ifdef _OPENMP
#pragma omp parallel for private(ix, iy, iz)
#endif
    for (iy=0; iy < ny; iy++) {  
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            old[iy][ix][iz]  = (-1.0*dt*dt*v[iy][ix][iz]*v[iy][ix][iz])* (
			   	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uuxx[iy][ix][iz]
			+	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uuyy[iy][ix][iz]
			+	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz]))*uuzz[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uuxy[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uuxz[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uuyz[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxx[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyyy[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uzzzz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxy[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uxyyy[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyyz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxzzz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyzzz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uxxyy[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxzz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyzz[iy][ix][iz]
			+	( 12.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxyz[iy][ix][iz]
			+	( 12.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uxyyz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uxyzz[iy][ix][iz]
			) + 2.0*cur[iy][ix][iz]-old[iy][ix][iz];

		}
	}  
	}

/*

fprintf(stderr, "\n Begin Old \n");
	for (iky=0; iky<ny;iky++) {
	for (ikx=0; ikx<nx;ikx++) {
		for (ikz=0;ikz < nz;ikz++){
			fprintf(stderr, "%f ", old[iky][ikx][ikz]);
			if (ikz==nkz-1)
				fprintf(stderr, "\n");
		}
	}
	}
fprintf(stderr, "\n End Old \n");

*/

}

