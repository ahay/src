/* Correct time-shift gathers with two-step lowrank propagator*/
/*
 Copyright (C) 2014 University of Texas at Austin
 
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
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fftwave.h"
#include "fft2w.h"
#include "absorb.h"

void zero2( float **data, int n1, int n2)
{
  int i1, i2;
#ifdef _OPENMP
#pragma omp parallel for private(i2, i1)
#endif
  for(i2=0; i2<n2; i2++){
    for(i1=0; i1<n1; i1++){
      data[i2][i1]=0.0;
    }
  }
}

void zero1( float *data, int n1)
{
  int i1;
#ifdef _OPENMP
#pragma omp parallel for private(i1)
#endif
  for(i1=0; i1<n1; i1++){
    data[i1]=0.0;
  }
}

int main(int argc, char *argv[])
{
  bool cmplx, abc, taper;
  int ix, iz, it, itau, n2, m2;
  int nx, nz, ntau, nt, pad1, nb, nx2, nz2, nzx2, fnx, fnz, fnzx, nk;
  float dt, dtau, par, dz, dx, thresh;
  float tau0, tau;
	
  float **lt, **rt;
  float *curr, *prev;
    
  float **dd, **mm, **v0, ***der;

  sf_axis ax, az, atau;
  sf_file tgather, cgather, left, right;
  FILE *swap;
    
  int cpuid, numprocs, nth;
  MPI_Comm comm=MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &cpuid);
  MPI_Comm_size(comm, &numprocs);
    
  sf_init(argc, argv);

#ifdef _OPENMP
#pragma omp parallel
  {
    nth=omp_get_num_threads();
  }
  sf_warning(">>> Using %d threads <<<", nth);
#endif
    
  tgather=sf_input("--input");
  cgather=sf_output("--output");
  left=sf_input("left");
  right=sf_input("right");
    
  az=sf_iaxa(tgather, 1);
  ax=sf_iaxa(tgather, 2);
  atau=sf_iaxa(tgather, 3);
    
  nz=sf_n(az); dz = sf_d(az); 
  nx=sf_n(ax); dx = sf_d(ax);
  ntau=sf_n(atau); dtau=sf_d(atau); tau0=sf_o(atau);
    
    if(cpuid==0){
      sf_oaxa(cgather, az, 1);
      sf_oaxa(cgather, ax, 2);
      sf_oaxa(cgather, atau, 3);
    }

    if(!sf_getfloat("dt", &dt)) dt=0.001;
    if(!sf_getint("nb", &nb)) nb=60;
    if(!sf_getfloat("par", &par)) par=0.01;
    if(!sf_getbool("cmplx", &cmplx)) cmplx=false; /* use complex FFT */
    if(!sf_getbool("abc", &abc)) abc=true; /* use complex FFT */
    if(!sf_getint("pad1", &pad1)) pad1=1; /* padding factor on the first axis */
    if(!sf_getbool("taper", &taper)) taper=true; /* tapering */
    if(!sf_getfloat("thresh", &thresh)) thresh=0.92; /* threshold */
    
    nx2=nx+2*nb;
    nz2=nz+2*nb;
    
    nk=fft2_init(cmplx,pad1,nz2,nx2,&fnz,&fnx);
    
    nzx2=nz2*nx2;
    fnzx=fnz*fnx;

    if (!sf_histint(left,"n1",&n2) || n2 != nzx2) sf_error("Need n1=%d in left",nzx2);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_floatalloc2(nzx2,m2);
    rt = sf_floatalloc2(m2,nk);

    sf_floatread(lt[0],nzx2*m2,left);
    sf_floatread(rt[0],m2*nk,right);
   
    dd=sf_floatalloc2(nz, nx);
    mm=sf_floatalloc2(nz, nx);
    der=sf_floatalloc3(nz, nx, 2);
    v0=sf_floatalloc2(nz2, nx2);
    curr=sf_floatalloc(fnzx);
    prev=sf_floatalloc(fnzx);
    
    /* initialize functions */
    lrinit(nx2, nz2, m2, cmplx, pad1, nb, par, abc, lt, rt, taper, thresh, dz, dx);
    
    swap=fopen("temswap.bin", "wb+");
    /* tau loop */
    for(itau=cpuid; itau<ntau; itau+=numprocs){
      sf_warning("itau=%d/%d", itau+1, ntau);
        
      tau=tau0+itau*dtau;

      sf_seek(tgather, itau*nz*nx*sizeof(float), SEEK_SET);
      sf_floatread(dd[0], nz*nx, tgather);

      // tau=0
      if(itau==(ntau-1)/2){
	for(ix=0; ix<nx; ix++)
	  for(iz=0; iz<nz; iz++)
	    mm[ix][iz]=dd[ix][iz];
	
	fseeko(swap, itau*nx*nz*sizeof(float), SEEK_SET);
    fwrite(mm[0], sizeof(float)*nz, nx, swap);
	continue;
      }
        
      // calculate v0 (derivative with respect to tau)
      zero2(v0, nz2, nx2);
      if(itau==0){
	sf_seek(tgather, itau*nz*nx*sizeof(float), SEEK_SET);
	sf_floatread(der[0][0], 2*nz*nx, tgather);

	for(ix=0; ix<nx; ix++){
	  for(iz=0; iz<nz; iz++){
	    v0[ix+nb][iz+nb]=(der[1][ix][iz]-der[0][ix][iz])/dtau;
	  }
	}
      } else if (itau==ntau-1){
	sf_seek(tgather, (itau-1)*nz*nx*sizeof(float), SEEK_SET);
	sf_floatread(der[0][0], 2*nz*nx, tgather);

	for(ix=0; ix<nx; ix++){
	  for(iz=0; iz<nz; iz++){
	    v0[ix+nb][iz+nb]=(der[1][ix][iz]-der[0][ix][iz])/dtau;
	  }
	}
      } else {
	sf_seek(tgather, (itau-1)*nz*nx*sizeof(float), SEEK_SET);
	sf_floatread(der[0][0], nz*nx, tgather);
	sf_seek(tgather, (itau+1)*nz*nx*sizeof(float), SEEK_SET);
	sf_floatread(der[1][0], nz*nx, tgather);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for(ix=0; ix<nx; ix++){
	  for(iz=0; iz<nz; iz++){
	    v0[ix+nb][iz+nb]=(der[1][ix][iz]-der[0][ix][iz])/dtau/2.0;
	  }
	}
      }
        
      // calculate u1
      zero1(curr, fnzx);
      zero1(prev, fnzx);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
      for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	  curr[(ix+nb)*fnz+iz+nb]=dd[ix][iz];
        
      // tau>0
      if(itau>(ntau-1)/2){
	// calculate u(t-lt)
	for (ix=0; ix<nx2; ix++)
	  for (iz=0; iz<nz2; iz++)
	    prev[ix*fnz+iz]=2.*v0[ix][iz]*dt;

	lrupdate(curr,prev);

	for (ix=0; ix<fnzx; ix++)
	  curr[ix]=curr[ix]/2.;
	// it loop
	nt=tau/dt+0.5;
	for(it=1; it<nt; it++){
	  sf_warning("it=%d/%d;", it+1, nt);

	  lrupdate(curr,prev);
	} //end of it
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for(ix=0; ix<nx; ix++){
	  for(iz=0; iz<nz; iz++){
	    mm[ix][iz]=curr[(ix+nb)*fnz+iz+nb];
	  }
	}
      } // end of positive tau
        
	// tau<0
      if(itau<(ntau-1)/2){
	//calculate u(t+lt)
	for (ix=0; ix<nx2; ix++)
	  for(iz=0; iz<nz2; iz++)
	    prev[ix*fnz+iz]=-2.*v0[ix][iz]*dt;

	lrupdate(curr, prev);

	for (ix=0; ix<fnzx; ix++)
	  curr[ix]=curr[ix]/2.;

	// it loop
	nt=-tau/dt+0.5;
	for(it=1; it<nt; it++){
	  sf_warning("it=%d/%d;", it+1, nt);

	  lrupdate(curr, prev);
	}//end of it
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for(ix=0; ix<nx; ix++){
	  for(iz=0; iz<nz; iz++){
	    mm[ix][iz]=curr[(ix+nb)*fnz+iz+nb];
	  }
	}
      }// end of negative tau

      // output the corrected time-shift gather
      fseeko(swap, itau*nx*nz*sizeof(float), SEEK_SET);
      fwrite(mm[0], sizeof(float)*nz, nx, swap);
    } //end of itau
    MPI_Barrier(comm);

    if(cpuid==0){
      for (itau=0; itau<ntau; itau++){
	fseeko(swap, itau*nz*nx*sizeof(float), SEEK_SET);
	fread(mm[0], sizeof(float)*nz, nx, swap);
	sf_floatwrite(mm[0], nz*nx, cgather);
      }
      fclose(swap);
    }
    MPI_Barrier(comm);
    
    //lrclose();
}
