/* Visco-acoustic gradient */
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
#include "Qfwi_commons.h"
#include "triutil.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <stdlib.h>
/*^*/

static bool verb, first;
static int cpuid, numprocs, nturn;
static int nz, nx, nzx, padnz, padnx, padnzx, nb, nt;
static int ns, ds_v, s0_v, sz, nr, dr_v, rz;
static int *nr2, *r02, *r0_v;
static int rectx, rectz, grectx, grectz, interval, wnt;
static int waterz, waterzb, wtn1, wtn2, woffn1, woffn2;

static float dt, idt, dt2, dx2, dz2, wdt, wdt2, scaling;
static float wt1, wt2, woff1, woff2;
static float ***dd, **vv, **tau, **taus, *ww, *bc, **weight;

MPI_Comm comm=MPI_COMM_WORLD;

void gradient_init(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, bool verb1)
/*< initialize >*/
{
	int iturn, is;

	verb=verb1;
	first=true; // only at the first iteration, need to calculate the gradient scaling parameter

	cpuid=mpipar->cpuid;
	numprocs=mpipar->numprocs;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;
	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nb=acpar->nb;
	nt=acpar->nt;

	ns=acpar->ns;
	ds_v=acpar->ds_v;
	s0_v=acpar->s0_v;
	sz=acpar->sz;
	nr=acpar->nr;
	dr_v=acpar->dr_v;
	nr2=acpar->nr2;
	r02=acpar->r02;
	r0_v=acpar->r0_v;
	rz=acpar->rz;

	rectx=soupar->rectx;
	rectz=soupar->rectz;
	grectx=fwipar->rectx;
	grectz=fwipar->rectz;
	interval=acpar->interval;
	wnt=(nt-1)/interval+1;

	dt=acpar->dt;
	idt=1./dt;
	dt2=dt*dt;
	wdt=dt*interval;
	wdt2=wdt*wdt;
	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;

	wt1=fwipar->wt1;
	wt2=fwipar->wt2;
	woff1=fwipar->woff1;
	woff2=fwipar->woff2;
	waterz=fwipar->waterz;

	ww=array->ww;
	bc=acpar->bc;
	
	/* read data */
	if(ns%numprocs==0) nturn=ns/numprocs;
	else nturn=ns/numprocs+1;
	dd=sf_floatalloc3(nt, nr, nturn);
	memset(dd[0][0], 0., nt*nr*nturn*sizeof(float));
	for(iturn=0; iturn<nturn; iturn++){
		is=iturn*numprocs+cpuid;
		if(is<ns){
			sf_seek(Fdat, is*nr*nt*sizeof(float), SEEK_SET);
			sf_floatread(dd[iturn][0], nr*nt, Fdat);
		}
	}
	
	/* data residual weights */
	wtn1=(wt1-acpar->t0)/dt+0.5;
	wtn2=(wt2-acpar->t0)/dt+0.5;
	woffn1=(woff1-acpar->r0)/acpar->dr+0.5;
	woffn2=(woff2-acpar->r0)/acpar->dr+0.5;
	weight=sf_floatalloc2(nt, nr);
	residual_weighting(weight, nt, nr, wtn1, wtn2, woffn1, woffn2, fwipar->oreo);

	/* padding and convert vector to 2-d array */
	vv = sf_floatalloc2(padnz, padnx);
	tau= sf_floatalloc2(padnz, padnx);
	taus=sf_floatalloc2(padnz, padnx);
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

	return;
}

void gradient_av(float *x, float *fcost, float *grad)
/*< acoustic velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **p0, **p1, **p2, **term, **tmparray, *rr, ***wave, **pp;
	float *sendbuf, *recvbuf;

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	/* memory allocation */
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wave=sf_floatalloc3(nz, nx, wnt);
	pp=sf_floatalloc2(nt, nr);

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(pp[0], 0., nr*nt*sizeof(float));
		
		sx=s0_v+is*ds_v;
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* output predicted data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p1[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
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
			shared(term,rr,padnx,padnz,ww,it)
#endif
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					term[ix][iz] += rr[ix*padnz+iz]*ww[it];
				}
			}

			/* update */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p0,p1,p2,vv,term,padnx,padnz,dt2)
#endif
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=dd[iturn][ir][it]-pp[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		
		/* window the data residual */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it] *= weight[ir][it];
			}
		}
		iturn++;

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
			
			/* load data residual*/
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

			/* update */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p0,p1,p2,vv,term,padnx,padnz,dt2)
#endif
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) \
			shared(nx,nz,vv,wave,p1,wit,wdt2,grad)
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=vv[ix+nb][iz+nb];
							temp=temp*temp*temp;
							temp=-2./temp;
							grad[ix*nz+iz] += (wave[wit+1][ix][iz]-2.*wave[wit][ix][iz]+wave[wit-1][ix][iz])/wdt2*p1[ix+nb][iz+nb]*temp;
						}
					}
				}
				wit--;
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=grad;
	}else{
		sendbuf=grad;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(grad, nzx, MPI_FLOAT, 0, comm);

	/* scaling gradient */
	if(first){
		dmax=0.;
		for(ix=0; ix<nzx; ix++)
			if(fabsf(grad[ix])>dmax)
				dmax=fabsf(grad[ix]);
		scaling=0.1/dmax;
		first=false;
	}

	/* smooth gradient */
	gradient_smooth2(grectx, grectz, nx, nz, waterz, scaling, grad);

	/* free allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*pp); free(pp);
	free(**wave); free(*wave); free(wave);
	free(rr); free(*term); free(term);
}

void gradient_v(float *x, float *fcost, float *grad)
/*< velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **p0, **p1, **p2, **r1, **r2, **term, **tmp, **tmparray, *rr, ***wave, **pp;
	float *sendbuf, *recvbuf;

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	/* memory allocation */
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	r1=sf_floatalloc2(padnz, padnx);
	r2=sf_floatalloc2(padnz, padnx);
	tmp=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wave=sf_floatalloc3(nz, nx, wnt);
	pp=sf_floatalloc2(nt, nr);

	
	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		memset(tmp[0], 0., padnzx*sizeof(float));
		
		sx=s0_v+is*ds_v;
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* output predicted data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p1[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(nz,nx,wit,nb,p1,wave)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[wit][ix][iz]=p1[ix+nb][iz+nb];
				wit++;
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* calculate r, load source and update wavefield */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(r2,taus,tau,term,p0,p1,p2)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					r2[ix][iz]=
						(tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (idt-0.5/taus[ix][iz])*r1[ix][iz])
						/(idt+0.5/taus[ix][iz]);
					term[ix][iz]=term[ix][iz]*(1.+tau[ix][iz]) - (r2[ix][iz]+r1[ix][iz])*0.5 + rr[ix*padnz+iz]*ww[it];
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;
			tmparray=r1; r1=r2; r2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
			apply_sponge(r1, bc, padnx, padnz, nb);
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=dd[iturn][ir][it]-pp[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
				pp[ir][it] *= weight[ir][it];
			}
		}
		iturn++;

		/* initialization */
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		memset(tmp[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);

			/* calculate and load r term */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(tau,taus,p1,r1,tmp)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					r2[ix][iz]=
						(-tau[ix][iz]/taus[ix][iz]*p1[ix][iz]
						 + (-idt+0.5/taus[ix][iz])*r1[ix][iz])
						/(-idt-0.5/taus[ix][iz]);
					tmp[ix][iz]=p1[ix][iz]*(1.+tau[ix][iz]) - 0.5*(r2[ix][iz]+r1[ix][iz]);
				}
			}

			/* laplacian operator */
			laplace(tmp, term, padnx, padnz, dx2, dz2);
			
			/* load data residual*/
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

			/* update */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(p0,p1,term)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) \
			shared(vv,nz,nx,nb,wave,p1)
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=vv[ix+nb][iz+nb];
							temp=temp*temp*temp;
							temp=-2./temp;
							grad[ix*nz+iz] += (wave[wit+1][ix][iz]-2.*wave[wit][ix][iz]+wave[wit-1][ix][iz])/wdt2*p1[ix+nb][iz+nb]*temp;
						}
					}
				}
				wit--;
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;
			tmparray=r1; r1=r2; r2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
			apply_sponge(r1, bc, padnx, padnz, nb);
		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=grad;
	}else{
		sendbuf=grad;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(grad, nzx, MPI_FLOAT, 0, comm);

	/* scaling gradient */
	if(first){
		dmax=0.;
		for(ix=0; ix<nzx; ix++)
			if(fabsf(grad[ix])>dmax)
				dmax=fabsf(grad[ix]);
		scaling=0.1/dmax;
		first=false;
	}

	/* smooth gradient */
	gradient_smooth2(grectx, grectz, nx, nz, waterz, scaling, grad);

	/* free allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*pp); free(pp);
	free(*r1); free(r1); free(*r2); free(r2);
	free(**wave); free(*wave); free(wave);
	free(rr); free(*term); free(term);
	free(*tmp); free(tmp);
}

void lstri_op(float **dd, float **dwt, float ***ww, float ***mwt, sf_acqui acpar, sf_vec array, sf_pas paspar, bool verb)
/*< ls TRI operator >*/
{
    float **vv1;
    int ix,iz,it;

    if (paspar->inv) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
        for         (it=0; it<acpar->nt; it++)
            for     (ix=0; ix<acpar->nx; ix++)
                for (iz=0; iz<acpar->nz; iz++)
                    ww[it][ix][iz] = 0.0f;
        if (NULL!=mwt) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
            for         (it=0; it<acpar->nt; it++)
                for     (ix=0; ix<acpar->nx; ix++)
                    for (iz=0; iz<acpar->nz; iz++)
                        mwt[it][ix][iz] = 1.0f;
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,it)
#endif
        for     (ix=0; ix<acpar->nx; ix++)
            for (it=0; it<acpar->nt; it++)
                dd[ix][it] = 0.0f;
    }

    /* map 1d to 2d */
    vv1 = (float**) sf_alloc (acpar->nx,sizeof(float*)); 
    vv1[0] = array->vv;
    for (ix=1; ix<acpar->nx; ix++) vv1[ix] = vv1[0]+ix*acpar->nz;

    timerev_init(verb, true, acpar->nt, acpar->nx, acpar->nz, acpar->nb, acpar->rz-acpar->nb, acpar->dt, acpar->dx, acpar->dz, acpar->coef, vv1);

    /* calculate model weighting using correlative imaging condition */
    if (paspar->inv && paspar->prec) { 
        if (paspar->ctr) {
            ctimerev(paspar->ngrp,mwt,dd);
            absval(acpar->nz*acpar->nx*acpar->nt,mwt[0][0]);
        } else {
            timerev_lop(true, false, acpar->nz*acpar->nx*acpar->nt, acpar->nt*acpar->nx, mwt[0][0], dd[0]);
            autopow(acpar->nz*acpar->nx*acpar->nt,(float)paspar->ngrp,mwt[0][0]);
        }
        /* smoothing */
        smooth(acpar->nz, acpar->nx, acpar->nt, paspar->rectz, paspar->rectx, paspar->rectt, paspar->repeat, mwt[0][0]);
        /* local normalizaiton */
        swnorm(verb, paspar->sw, acpar->nz, acpar->nx, acpar->nt, paspar->size, paspar->perc, mwt[0][0]);
        /* hard thresholding */
        if (paspar->hard>0) threshold(false, acpar->nz*acpar->nx*acpar->nt, paspar->hard, mwt[0][0]);
    }

    /* apply time-reversal imaging linear operator */
    if (paspar->inv) {
        if (NULL!=dwt) sf_solver(timerev_lop,sf_cgstep,acpar->nz*acpar->nx*acpar->nt,acpar->nt*acpar->nx,ww[0][0],dd[0],paspar->niter,"mwt",mwt[0][0],"wt",dwt[0],"verb",verb,"end");
        else sf_solver(timerev_lop,sf_cgstep,acpar->nz*acpar->nx*acpar->nt,acpar->nt*acpar->nx,ww[0][0],dd[0],paspar->niter,"mwt",mwt[0][0],"verb",verb,"end");
    } else {
        timerev_lop(false, false, acpar->nz*acpar->nx*acpar->nt, acpar->nt*acpar->nx, ww[0][0], dd[0]);
    }
    
    /* close */
    timerev_close();
    free(vv1);

}

static float ****ww3;
static float ****gwt;
static sf_butter blo=NULL, bhi=NULL;

/* for passive source and fwi */
void gradient_pas_init(sf_file Fdat, sf_file Fsrc, sf_file Fmwt, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_pas paspar, bool verb1)
/*< initialize >*/
{
        float **dwt=NULL,***mwt,***wwt;
        int it,ix,iz,iturn,is,rdn;
        char filename[20]="tempbin",srdn[10];
        FILE *temp;

	verb=verb1;
	first=true; // only at the first iteration, need to calculate the gradient scaling parameter

	cpuid=mpipar->cpuid;
	numprocs=mpipar->numprocs;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;
	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nb=acpar->nb;
	nt=acpar->nt;

	ns=acpar->ns;
	nr=acpar->nr;
	dr_v=acpar->dr_v;
	r0_v=acpar->r0_v;
	rz=acpar->rz;

	grectx=fwipar->rectx;
	grectz=fwipar->rectz;
	interval=acpar->interval;
	wnt=(nt-1)/interval+1;

	dt=acpar->dt;
	idt=1./dt;
	dt2=dt*dt;
	wdt=dt*interval;
	wdt2=wdt*wdt;
	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;

	wt1=fwipar->wt1;
	wt2=fwipar->wt2;
	woff1=fwipar->woff1;
	woff2=fwipar->woff2;
	waterz=fwipar->waterz;
	waterzb=fwipar->waterzb;

	bc=acpar->bc;

        if (cpuid==0) {
            srand(time(NULL));
            rdn = rand()%1000000000;
            sprintf(srdn,"%d",rdn);
            strcat(filename,srdn);
        }
	MPI_Bcast(filename, 20, MPI_CHAR, 0, comm);
        if(verb && cpuid==0) sf_warning("filename=%s",filename);

        temp=fopen(filename, "wb+");

	if(ns%numprocs==0) nturn=ns/numprocs;
	else nturn=ns/numprocs+1;

        /* allocate data/source/weight */
        dd  = sf_floatalloc3(nt, nx, nturn);
        ww3 = sf_floatalloc4(nz, nx, nt, nturn);
        gwt = sf_floatalloc4(nz, nx, nt, nturn);
        wwt = sf_floatalloc3(nz, nx, nt); /* temporary output var */
        if (!paspar->onlyvel) {
            mwt = sf_floatalloc3(nz, nx, nt); /* src model weight */
            /*
            dwt = sf_floatalloc2(acpar->nt, acpar->nx);

            wtn1=(fwipar->wt1-acpar->t0)/acpar->dt+0.5;
            wtn2=(fwipar->wt2-acpar->t0)/acpar->dt+0.5;
            woffn1=(fwipar->woff1-acpar->r0)/acpar->dr+0.5;
            woffn2=(fwipar->woff2-acpar->r0)/acpar->dr+0.5;
            residual_weighting(dwt, acpar->nt, acpar->nx, wtn1, wtn2, woffn1, woffn2, !fwipar->oreo);
            */
        } else {
            mwt=NULL;
            dwt=NULL;
        }

        /* read data/source */
        for(iturn=0; iturn<nturn; iturn++){
            is=iturn*numprocs+cpuid;
            if(is<ns){
                /* read data */
                sf_seek(Fdat, is*nt*nx*sizeof(float), SEEK_SET);
                sf_floatread(dd[iturn][0], nt*nx, Fdat);
                if (paspar->onlyvel) {
                    /* read source */
                    sf_seek(Fsrc, is*nz*nx*nt*sizeof(float), SEEK_SET);
                    sf_floatread(ww3[iturn][0][0], nz*nx*nt, Fsrc);
                } else {
                    /* linear inversion of source */
                    lstri_op(dd[iturn], dwt, ww3[iturn], mwt, acpar, array, paspar, verb);
                    /* write source */
                    fseeko(temp, is*nz*nx*nt*sizeof(float), SEEK_SET);
                    fwrite(ww3[iturn][0][0], sizeof(float), nz*nx*nt, temp);
                    if (NULL!=Fmwt && is==0) sf_floatwrite(mwt[0][0], nz*nx*nt, Fmwt);
                }

                /* calculate gradient mask */
                if (!paspar->onlyvel && paspar->prec && paspar->hidesrc) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
                    for         (it=0; it<nt; it++)
                        for     (ix=0; ix<nx; ix++)
                            for (iz=0; iz<nz; iz++)
                                gwt[iturn][it][ix][iz] = mwt[it][ix][iz];
                    threshold(true, nz*nx*nt, paspar->hard, gwt[iturn][0][0]);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
                    for         (it=0; it<nt; it++)
                        for     (ix=0; ix<nx; ix++)
                            for (iz=0; iz<nz; iz++)
                                gwt[iturn][it][ix][iz] = 1.-gwt[iturn][it][ix][iz];
                } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
                    for         (it=0; it<nt; it++)
                        for     (ix=0; ix<nx; ix++)
                            for (iz=0; iz<nz; iz++)
                                gwt[iturn][it][ix][iz] = 1.;
                }
            } /* if is<ns */
        }
        fclose(temp);
        MPI_Barrier(comm);

        if(!paspar->onlyvel && cpuid==0) {
            temp=fopen(filename, "rb");
            for(is=0; is<ns; is++){
                fseeko(temp, is*nz*nx*nt*sizeof(float), SEEK_SET);
                fread(wwt[0][0], sizeof(float), nz*nx*nt, temp);
                sf_floatwrite(wwt[0][0], nz*nx*nt, Fsrc);
            }
            fclose(temp);
            remove(filename);
        }
        MPI_Barrier(comm);

	/* data residual weights */
	wtn1=(wt1-acpar->t0)/dt+0.5;
	wtn2=(wt2-acpar->t0)/dt+0.5;
	woffn1=(woff1-acpar->r0)/acpar->dr+0.5;
	woffn2=(woff2-acpar->r0)/acpar->dr+0.5;
	weight=sf_floatalloc2(nt, nr);
	residual_weighting(weight, nt, nr, wtn1, wtn2, woffn1, woffn2, fwipar->oreo);

	/* padding and convert vector to 2-d array */
	vv = sf_floatalloc2(padnz, padnx);
	tau= sf_floatalloc2(padnz, padnx);
	taus=sf_floatalloc2(padnz, padnx);
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

        /* multiscale gradient */
	if(soupar->flo > 0.0001) blo=sf_butter_init(false, soupar->flo, 3);
	if(soupar->fhi < 0.5-0.0001) bhi=sf_butter_init(true, soupar->fhi, 3);

        free(**wwt); free(*wwt); free(wwt);
        if (NULL!=mwt) { free(**mwt); free(*mwt); free(mwt); }
	return;
}

//JS
//static int counter=0;
void gradient_pas_av(float *x, float *fcost, float *grad)
/*< acoustic velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int rx;

	float temp, dmax;
	float **p0, **p1, **p2, **term, **tmparray, ***wave, **pp;
	float *sendbuf, *recvbuf;

        /*
        //JS
        sf_file Fwfl1, Fwfl2, Fres;
        counter++;
        if (cpuid==0 && counter==3) {
            Fwfl1=sf_output("Fwfl1");
            Fwfl2=sf_output("Fwfl2");
            Fres=sf_output("Fres");
            sf_putint(Fwfl1,"n1",padnz);
            sf_putint(Fwfl1,"n2",padnx);
            sf_putint(Fwfl1,"n3",(nt-1)/50+1);
            sf_putint(Fwfl2,"n1",padnz);
            sf_putint(Fwfl2,"n2",padnx);
            sf_putint(Fwfl2,"n3",(nt-1)/50+1);
            sf_putint(Fres,"n1",nt);
            sf_putint(Fres,"n2",nr);
        }
        */

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	/* memory allocation */
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	wave=sf_floatalloc3(nz, nx, wnt);
	pp=sf_floatalloc2(nt, nr);

        iturn=0;
        for(is=cpuid; is<ns; is+=numprocs){
            if(cpuid==0) sf_warning("###### is=%d ######", is+1);

            memset(p0[0], 0., padnzx*sizeof(float));
            memset(p1[0], 0., padnzx*sizeof(float));
            memset(p2[0], 0., padnzx*sizeof(float));
            memset(pp[0], 0., nr*nt*sizeof(float));

            wit=0;
            /* forward propagation */
            for(it=0; it<nt; it++){
                if(verb) sf_warning("Forward propagation it=%d;", it);

                /* output predicted data */
                for(ir=0; ir<nr; ir++){
                    rx=r0_v[0]+ir*dr_v;
                    pp[ir][it]=p1[rx][rz];
                }

                /* save wavefield */
                if(it%interval==0){
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

                /*
                //JS
                if(is==0 && counter==3 && it%50==0) sf_floatwrite(p1[0],padnzx,Fwfl1);
                */

                /* laplacian operator */
                laplace(p1, term, padnx, padnz, dx2, dz2);

                /* load source */
#ifdef _OPENMP 
#pragma omp parallel for \
                private(ix,iz) \
                shared(term,nb,ww3,it)
#endif
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        term[ix+nb][iz+nb] += ww3[iturn][it][ix][iz];
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
                apply_sponge(p0, bc, padnx, padnz, nb);
                apply_sponge(p1, bc, padnx, padnz, nb);
            } // end of time loop

            /* check */
            if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
            wit--;

            /* calculate data residual and data misfit */
            for(ir=0; ir<nr; ir++){
                for(it=0; it<nt; it++){
                    pp[ir][it]=dd[iturn][ir][it]-pp[ir][it];
                    *fcost += 0.5*pp[ir][it]*pp[ir][it];
                }
            }

            /* window the data residual */
            for(ir=0; ir<nr; ir++){
                /* multiscale */
                if(NULL != blo){
                    sf_butter_apply(blo, nt, pp[ir]);
                    sf_reverse(nt, pp[ir]);
                    sf_butter_apply(blo, nt, pp[ir]);
                    sf_reverse(nt, pp[ir]);
                }
                if(NULL != bhi){
                    sf_butter_apply(bhi, nt, pp[ir]);
                    sf_reverse(nt, pp[ir]);
                    sf_butter_apply(bhi, nt, pp[ir]);
                    sf_reverse(nt, pp[ir]);
                }
                for(it=0; it<nt; it++){
                    pp[ir][it] *= weight[ir][it];
                }
            }

            /*
            // JS
            if(is==0 && counter==3) sf_floatwrite(pp[0], nr*nt, Fres);
            */

            /* initialization */
            memset(p0[0], 0., padnzx*sizeof(float));
            memset(p1[0], 0., padnzx*sizeof(float));
            memset(p2[0], 0., padnzx*sizeof(float));
            memset(term[0], 0., padnzx*sizeof(float));

            /* backward propagation */
            for(it=nt-1; it>=0; it--){
                if(verb) sf_warning("Backward propagation it=%d;", it);

                /* laplacian operator */
                laplace(p1, term, padnx, padnz, dx2, dz2);

                /* load data residual*/
                for(ir=0; ir<nr; ir++){
                    rx=r0_v[0]+ir*dr_v;
                    term[rx][rz] += pp[ir][it];
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

                /*
                // JS
                if(is==0 && counter==3 && it%50==0) sf_floatwrite(p1[0],padnzx,Fwfl2);
                */

                /* calculate gradient  */
                if(it%interval==0){
                    if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
                        private(ix,iz,temp) \
                        shared(nx,nz,vv,wave,p1,wit,wdt2,grad)
#endif
                        for(ix=0; ix<nx; ix++){
                            for(iz=0; iz<nz; iz++){
                                temp=vv[ix+nb][iz+nb];
                                temp=temp*temp*temp;
                                temp=-2./temp;
                                grad[ix*nz+iz] += gwt[iturn][it][ix][iz]*(wave[wit+1][ix][iz]-2.*wave[wit][ix][iz]+wave[wit-1][ix][iz])/wdt2*p1[ix+nb][iz+nb]*temp;
                            }
                        }
                    }
                    wit--;
                }

                /* swap wavefield pointer of different time steps */
                tmparray=p0; p0=p1; p1=p2; p2=tmparray;

                /* boundary condition */
                apply_sponge(p0, bc, padnx, padnz, nb);
                apply_sponge(p1, bc, padnx, padnz, nb);
            } // end of time loop

            iturn++;

            /*
            // JS
            if(is==0 && counter==3){
            sf_fileclose(Fwfl1);
            sf_fileclose(Fwfl2);
            sf_fileclose(Fres);
            }
            sf_warning("---counter=%d fcost=%3.3e---",counter, *fcost);
            */
        } // end of shot loop
	MPI_Barrier(comm);

	/* misfit reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=grad;
	}else{
		sendbuf=grad;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(grad, nzx, MPI_FLOAT, 0, comm);

        /* scaling gradient */
        if(first){
            dmax=0.;
            for(ix=0; ix<nzx; ix++)
                if(fabsf(grad[ix])>dmax)
                    dmax=fabsf(grad[ix]);
            scaling=0.1/dmax;
            first=false;
        }

	/* smooth gradient */
	gradient_smooth2b(grectx, grectz, nx, nz, waterz, waterzb, scaling, grad);

	/* free allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*pp); free(pp);
	free(**wave); free(*wave); free(wave);
	free(*term); free(term);
}

