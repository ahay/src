/* Visco-acoustic FWI */
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
#include "fwi_commons.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/*^*/

bool verb, first;
int cpuid, numprocs, nturn; // mpi-related
int nz, nx, nzx, padnz, padnx, padnzx, nb, nt; // dimension
int ns, ds_v, s0_v, sz, nr, dr_v, rz, *nr2, *r02, *r0_v; // acquisition
int frectx, frectz, interval, wnt; // wavefield
int waterz, wtn1, wtn2, woffn1, woffn2, grectx, grectz; // gradient

float dt, idt, dt2, dx2, dz2, wdt, wdt2; // wavefield
float wt1, wt2, woff1, woff2, gain, scaling, scaling2; // gradient
float ***dd, **vv, **tau, **taus, *ww, *bc, **weight; // arrays

MPI_Comm comm=MPI_COMM_WORLD;

void gradient_init(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_q array, sf_fwi_q fwipar, bool verb1)
/*< initialize >*/
{
	int iturn, is;

	verb=verb1;
	first=true; // only at the first iteration (for calculating the gradient scaling parameter)

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
	rz=acpar->rz;
	nr2=acpar->nr2;
	r02=acpar->r02;
	r0_v=acpar->r0_v;

	frectx=soupar->frectx;
	frectz=soupar->frectz;
	interval=acpar->interval;
	wnt=(nt-1)/interval+1;
	
	dt=acpar->dt;
	idt=1./dt;
	dt2=dt*dt;
	wdt=dt*interval;
	wdt2=wdt*wdt;
	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;

	/* gradient preconditioning */
	waterz=fwipar->waterz;
	grectx=fwipar->grectx;
	grectz=fwipar->grectz;
	
	/* data residual weighting */
	gain=fwipar->gain;
	wt1=fwipar->wt1;
	wt2=fwipar->wt2;
	woff1=fwipar->woff1;
	woff2=fwipar->woff2;
	wtn1=(wt1-acpar->t0)/dt+0.5;
	wtn2=(wt2-acpar->t0)/dt+0.5;
	woffn1=(woff1-acpar->r0)/acpar->dr+0.5;
	woffn2=(woff2-acpar->r0)/acpar->dr+0.5;
	weight=sf_floatalloc2(nt, nr);
	residual_weighting(weight, nt, nr, wtn1, wtn2, woffn1, woffn2, gain);

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
	
	/* padding and convert vector to 2-d array */
	vv = sf_floatalloc2(padnz, padnx);
	tau= sf_floatalloc2(padnz, padnx);
	taus=sf_floatalloc2(padnz, padnx);
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

	return;
}

void gradient_v(float *x, float *fcost, float *grad)
/*< velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **p0, **p1, **p2, **r1, **r2, **term, **tmp, **tmparray, *rr, ***wave, **pp;
	float *sendbuf, *recvbuf;
	//sf_file Fwfl1, Fwfl2, Fres, Fv, Ft;
	//Fwfl1=sf_output("Fwfl1");
	//Fwfl2=sf_output("Fwfl2");
	//Fres=sf_output("Fres");
	//Fv=sf_output("Fv");
	//Ft=sf_output("Ft");
	//counter++;

	/*
	sf_putint(Fwfl1,"n1",padnz);
	sf_putint(Fwfl1,"n2",padnx);
	sf_putint(Fwfl1,"n3",(nt-1)/50+1);
	sf_putint(Fwfl2,"n1",padnz);
	sf_putint(Fwfl2,"n2",padnx);
	sf_putint(Fwfl2,"n3",(nt-1)/50+1);
	sf_putint(Fres,"n1",nt);
	sf_putint(Fres,"n2",nr);
	sf_putint(Fv,"n1",padnz);
	sf_putint(Fv,"n2",padnx);
	sf_putint(Ft,"n1",padnz);
	sf_putint(Ft,"n2",padnx);
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
		
		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

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

			//if(is==0 && counter==2  && it%50==0) sf_floatwrite(p1[0],padnzx,Fwfl1);

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
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

		//if(is==0 && counter==2) sf_floatwrite(pp[0], nr*nt, Fres);
		//if(is==0 && counter==2) sf_floatwrite(vv[0], padnzx, Fv);
		//if(is==0 && counter==2) sf_floatwrite(tau[0], padnzx, Ft);

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
			
			//if(cpuid==0 && counter==2 && it%50==0) sf_floatwrite(p1[0],padnzx,Fwfl2);

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
	
	/*
	if(counter==2){
		sf_fileclose(Fwfl1);
		sf_fileclose(Fwfl2);
		sf_fileclose(Fres);
		sf_fileclose(Fv);
		sf_fileclose(Ft);
	}
	*/

	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
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
	//MPI_ALLreduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, comm);
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

void gradient_q(float *x, float *fcost, float *grad)
/*< tau gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float dmax;
	float **p0, **p1, **p2, **r1, **r2, **term, **tmp, **tmparray, *rr, ***wave, **pp;
	float *sendbuf, *recvbuf;

	/* initialize fcost */
	*fcost=0.;
	/* update tau */
	pad2d(x, tau, nz, nx, nb);
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
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* output predicted data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p1[rx][rz];
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* calculate r, load source and update wavefield */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(r2,taus,tau,term,tmp,p0,p1,p2)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					r2[ix][iz]=
						(tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (idt-0.5/taus[ix][iz])*r1[ix][iz])
						/(idt+0.5/taus[ix][iz]);
					tmp[ix][iz]=term[ix][iz]*(1.+tau[ix][iz]) - (r2[ix][iz]+r1[ix][iz])*0.5 + rr[ix*padnz+iz]*ww[it];
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*tmp[ix][iz];
				}
			}
			
			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(nz,nx,wit,nb,term,r1,r2,tau,wave)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[wit][ix][iz]= - (term[ix+nb][iz+nb] - (r1[ix+nb][iz+nb]+r2[ix+nb][iz+nb])*0.5/tau[ix+nb][iz+nb]);
				wit++;
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
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
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
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(nz,nx,wit,nb,wave,p1,grad)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						grad[ix*nz+iz] += wave[wit][ix][iz]*p1[ix+nb][iz+nb];
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
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
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
	//MPI_ALLreduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, comm);
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
		scaling=0.02/dmax;
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

void gradient_vq(float *x, float *fcost, float *grad)
/*< velocity and tau gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx, nm;

	float **p0, **p1, **p2, **r1, **r2, **term, **tmp, **tmparray, *rr, ***wave1, ***wave2, **pp;
	float *sendbuf, *recvbuf;
	float temp, dmax, *x1, *x2;

	/* initialize fcost */
	*fcost=0.;
	/* update velocity and tau */
	x1=x;
	x2=x+nzx;
	pad2d(x1, vv, nz, nx, nb);
	pad2d(x2, tau, nz, nx, nb);
	/* model dimension */
	nm=2*nzx;
	/* initialize gradient */
	memset(grad, 0., nm*sizeof(float));

	/* memory allocation */
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	r1=sf_floatalloc2(padnz, padnx);
	r2=sf_floatalloc2(padnz, padnx);
	tmp=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wave1=sf_floatalloc3(nz, nx, wnt);
	wave2=sf_floatalloc3(nz, nx, wnt);
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
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* output predicted data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p1[rx][rz];
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* calculate r, load source and update wavefield */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(r2,taus,tau,term,tmp,p0,p1,p2)
#endif
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					r2[ix][iz]=
						(tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (idt-0.5/taus[ix][iz])*r1[ix][iz])
						/(idt+0.5/taus[ix][iz]);
					tmp[ix][iz]=term[ix][iz]*(1.+tau[ix][iz]) - (r2[ix][iz]+r1[ix][iz])*0.5 + rr[ix*padnz+iz]*ww[it];
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*tmp[ix][iz];
				}
			}
			
			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(nz,nx,wit,nb,term,r1,r2,tau,p1,wave1,wave2)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++){
						wave1[wit][ix][iz]= p1[ix+nb][iz+nb];
						wave2[wit][ix][iz]= - (term[ix+nb][iz+nb] - (r1[ix+nb][iz+nb]+r2[ix+nb][iz+nb])*0.5/tau[ix+nb][iz+nb]);
					}
				wit++;
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
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
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
			private(ix,iz) \
			shared(nz,nx,wit,nb,wave1,wave2,p1,grad)
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=vv[ix+nb][iz+nb];
							temp=temp*temp*temp;
							temp=-2./temp;
							grad[ix*nz+iz] += (wave1[wit+1][ix][iz]-2.*wave1[wit][ix][iz]+wave1[wit-1][ix][iz])/wdt2*p1[ix+nb][iz+nb]*temp;
							grad[ix*nz+iz+nzx] += wave2[wit][ix][iz]*p1[ix+nb][iz+nb];
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
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
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
	//MPI_ALLreduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, comm);
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=grad;
	}else{
		sendbuf=grad;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nm, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(grad, nm, MPI_FLOAT, 0, comm);

	/* scaling gradient */
	if(first){
		dmax=0.;
		for(ix=0; ix<nzx; ix++)
			if(fabsf(grad[ix])>dmax)
				dmax=fabsf(grad[ix]);
		scaling=0.1/dmax;
		
		dmax=0.;
		for(ix=nzx; ix<nm; ix++)
			if(fabsf(grad[ix])>dmax)
				dmax=fabsf(grad[ix]);
		scaling2=0.02/dmax;
		first=false;
	}

	/* smooth gradient */
	gradient_smooth2(grectx, grectz, nx, nz, waterz, scaling, grad);
	gradient_smooth2(grectx, grectz, nx, nz, waterz, scaling2, grad+nzx);

	/* free allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*pp); free(pp);
	free(*r1); free(r1); free(*r2); free(r2);
	free(**wave1); free(*wave1); free(wave1);
	free(**wave2); free(*wave2); free(wave2);
	free(rr); free(*term); free(term);
	free(*tmp); free(tmp);
}

void fwi(sf_file Fdat, sf_file Finv, sf_file Ferr, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_q array, sf_fwi_q fwipar, sf_optim optpar, bool verb1)
/*< fwi >*/
{
	int i, iter=0, flag, nm=0, type=0;
	float fcost, threshold[4];
	float *x=NULL, *direction, *grad;
	sf_gradient gradient=NULL;
	FILE *fp=NULL;

	/* initialize */
	gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb1);

	/* gradient type */
	if(fwipar->grad_type==1){
		gradient=gradient_v;
		nm=nzx;
		x=array->vv;
		threshold[0]=fwipar->v1;
		threshold[1]=fwipar->v2;
		type=1;
	}else if(fwipar->grad_type==2){
		gradient=gradient_q;
		nm=nzx;
		x=array->tau;
		threshold[0]=fwipar->tau1;
		threshold[1]=fwipar->tau2;
		type=1;
	}else if(fwipar->grad_type==3){
		gradient=gradient_vq;
		nm=2*nzx;
		x=sf_floatalloc(nm);;
		for(i=0; i<nzx; i++){
			x[i]=array->vv[i];
			x[nzx+i]=array->tau[i];
		}
		threshold[0]=fwipar->v1;
		threshold[1]=fwipar->v2;
		threshold[2]=fwipar->tau1;
		threshold[3]=fwipar->tau2;
		type=2;
	}


	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	/* if onlygrad=y, program terminates */
	if(fwipar->onlygrad) return; 

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");

	direction=sf_floatalloc(nm);
	optpar->sk=sf_floatalloc2(nm, optpar->npair);
	optpar->yk=sf_floatalloc2(nm, optpar->npair);

	optpar->igrad=1;
	optpar->ipair=0;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
	optpar->alpha=1.;
	/* initialize data error vector */
	for(iter=0; iter<optpar->niter+1; iter++){
		optpar->err[iter]=0.;
	}
	optpar->err[0]=fcost;

	iter=0;
	if(mpipar->cpuid==0){
		l2norm(nm, grad, &optpar->gk_norm);
		print_iteration(fp, iter, optpar);
	}

	/* optimization loop */
	for(iter=0; iter<optpar->niter; iter++){
		if(mpipar->cpuid==0) sf_warning("-------------------iter=%d----------------------", iter+1);
		
		optpar->ils=0;
		if(iter%optpar->repeat==0) optpar->alpha=1.;

		/* search direction */
		if(iter==0){
			reverse(nm, grad, direction);
		}else{
			lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
			lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
		}

		/* line search */
		lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
		line_search(nm, x, grad, direction, gradient, optpar, threshold, &flag, mpipar->cpuid, type);
		optpar->err[iter+1]=optpar->fk;
		
		if(mpipar->cpuid==0){
			l2norm(nm, grad, &optpar->gk_norm);
			print_iteration(fp, iter+1, optpar);
			fclose(fp); /* get written to disk right away */
			fp=fopen("iterate.txt","a");
		}

		if(mpipar->cpuid==0 && flag==2){
			fprintf(fp, "Line Search Failed\n");
			break;
		}

		if(mpipar->cpuid==0 && optpar->fk/optpar->f0 < optpar->conv_error){
			fprintf(fp, "Convergence Criterion Reached\n");
			break;
		}
	} // end of iter

	if(mpipar->cpuid==0 && iter==optpar->niter){
		fprintf(fp, "Maximum Iteration Number Reached\n");
	}

	/* output vel & misfit */
	if(mpipar->cpuid==0) sf_floatwrite(x, nm, Finv);
	if(mpipar->cpuid==0) sf_floatwrite(optpar->err, optpar->niter+1, Ferr);
	if(mpipar->cpuid==0) fclose(fp);

	return;
}
