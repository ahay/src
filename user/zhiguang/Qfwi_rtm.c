/* Acoustic/Visco-acoustic RTM */
/*
 Copyright (C) 2016 The University of Texas at Austin
 
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
#include "qfwi_commons.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/*^*/

void rtm_a(sf_file Fdat, sf_file Fimg, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, bool verb)
/*< acoustic rtm >*/
{
	int ix, iz, is, ir, it, wit;
	int sx, rx, sz, rz, rectx, rectz;
	int nz, nx, nzx, padnz, padnx, padnzx, nt, nr, nb, wnt;

	float dx2, dz2, dt2, dt;
	float **vv, **dd, **mm;
	float **p0, **p1, **p2, **term, **tmparray, *rr, ***wave;
	float *sendbuf, *recvbuf;

	MPI_Comm comm=MPI_COMM_WORLD;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;
	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nr=acpar->nr;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	rectx=soupar->rectx;
	rectz=soupar->rectz;

	nt=acpar->nt;
	wnt=(nt-1)/acpar->interval+1;

	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;
	dt2=acpar->dt*acpar->dt;
	dt=acpar->dt;

	/* memory allocation */
	vv = sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);
	mm=sf_floatalloc2(nz, nx);

	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wave=sf_floatalloc3(nz, nx, wnt);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);

	memset(mm[0], 0., nzx*sizeof(float));

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* save wavefield */
			if(it%acpar->interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(wave,p1,wit,nb,nx,nz)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[wit][ix][iz]=p1[ix+nb][iz+nb];
				wit++;
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* load source */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(term,rr,padnx,padnz,it)
#endif
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					term[ix][iz] += rr[ix*padnz+iz]*array->ww[it];
				}
			}

			/* update */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p0,p1,p2,vv,term,padnx,padnz,dt2)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, acpar->bc, padnx, padnz, nb);
			apply_sponge(p1, acpar->bc, padnx, padnz, nb);
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;

		/* read data */
		sf_seek(Fdat, is*nr*nt*sizeof(float), SEEK_SET);
		sf_floatread(dd[0], nr*nt, Fdat);

		/* initialization */
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* load data */
			for(ir=0; ir<acpar->nr2[is]; ir++){
				rx=acpar->r0_v[is]+ir*acpar->dr_v;
				term[rx][rz] += dd[acpar->r02[is]+ir][it];
			}

			/* update */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p0,p1,p2,vv,term,padnx,padnz,dt2)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}

			/* calculate image */
			if(it%acpar->interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p1,wave,nx,nz,wit,mm,nb)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						mm[ix][iz] += wave[wit][ix][iz]*p1[ix+nb][iz+nb];
				wit--;
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, acpar->bc, padnx, padnz, nb);
			apply_sponge(p1, acpar->bc, padnx, padnz, nb);
		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);

	if(mpipar->cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=mm[0];
	}else{
		sendbuf=mm[0];
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);

	if(mpipar->cpuid==0) sf_floatwrite(mm[0], nzx, Fimg);
	MPI_Barrier(comm);
}

void rtm(sf_file Fdat, sf_file Fimg, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, bool verb)
/*< visco-acoustic rtm >*/
{
	int ix, iz, is, ir, it, wit;
	int sx, rx, sz, rz, rectx, rectz;
	int nz, nx, nzx, padnz, padnx, padnzx, nt, nr, nb, wnt;

	float dx2, dz2, dt2, dt;
	float **vv, **tau, **taus, **dd, **mm;
	float **p0, **p1, **p2, **r1, **r2, **term, **tmparray, *rr, ***wave;
	float *sendbuf, *recvbuf;

	MPI_Comm comm=MPI_COMM_WORLD;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;
	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nr=acpar->nr;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	rectx=soupar->rectx;
	rectz=soupar->rectz;

	nt=acpar->nt;
	wnt=(nt-1)/acpar->interval+1;

	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;
	dt2=acpar->dt*acpar->dt;
	dt=acpar->dt;

	/* memory allocation */
	vv = sf_floatalloc2(padnz, padnx);
	tau= sf_floatalloc2(padnz, padnx);
	taus=sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);
	mm=sf_floatalloc2(nz, nx);

	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	r1=sf_floatalloc2(padnz, padnx);
	r2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wave=sf_floatalloc3(nz, nx, wnt);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

	memset(mm[0], 0., nzx*sizeof(float));

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* save wavefield */
			if(it%acpar->interval==0){
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[wit][ix][iz]=p1[ix+nb][iz+nb];
				wit++;
			}

			/* load source */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					p1[ix][iz] += rr[ix*padnz+iz]*array->ww[it];
				}
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* update */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					r2[ix][iz]=
						(-tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (1./dt-0.5/taus[ix][iz])*r1[ix][iz])
						/(1./dt+0.5/taus[ix][iz]);
					term[ix][iz]=term[ix][iz]*(1.+tau[ix][iz])+(r2[ix][iz]+r1[ix][iz])*0.5;
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;
			tmparray=r1; r1=r2; r2=tmparray;

			/* boundary condition */
			apply_sponge(p0, acpar->bc, padnx, padnz, nb);
			apply_sponge(p1, acpar->bc, padnx, padnz, nb);
			apply_sponge(r1, acpar->bc, padnx, padnz, nb);
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		/* read data */
		sf_seek(Fdat, is*nr*nt*sizeof(float), SEEK_SET);
		sf_floatread(dd[0], nr*nt, Fdat);
		/* initialization */
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);

			/* load data */
			for(ir=0; ir<acpar->nr2[is]; ir++){
				rx=acpar->r0_v[is]+ir*acpar->dr_v;
				p1[rx][rz]=dd[acpar->r02[is]+ir][it];
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* update */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					r2[ix][iz]=
						(-tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (1./dt-0.5/taus[ix][iz])*r1[ix][iz])
						/(1./dt+0.5/taus[ix][iz]);
					term[ix][iz]=term[ix][iz]*(1.+tau[ix][iz])+(r2[ix][iz]+r1[ix][iz])*0.5;
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}

			/* calculate image */
			if(it%acpar->interval==0){
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						mm[ix][iz] += wave[wit-1][ix][iz]*p2[ix+nb][iz+nb];
				wit--;
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;
			tmparray=r1; r1=r2; r2=tmparray;

			/* boundary condition */
			apply_sponge(p0, acpar->bc, padnx, padnz, nb);
			apply_sponge(p1, acpar->bc, padnx, padnz, nb);
			apply_sponge(r1, acpar->bc, padnx, padnz, nb);
		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);

	if(mpipar->cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=mm[0];
	}else{
		sendbuf=mm[0];
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);

	if(mpipar->cpuid==0) sf_floatwrite(mm[0], nzx, Fimg);
	MPI_Barrier(comm);
}
