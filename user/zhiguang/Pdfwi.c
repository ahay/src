/* Acoustic FWI */
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

float dt, idt, dx, dz, wdt; // wavefield
float wt1, wt2, woff1, woff2, gain, v0, t0, scaling; // gradient
float ***dd, **vv, **den, *ww, *bcxp, *bczp, *bcxv, *bczv, **weight; // arrays

// reflection fwi gradient
int rfwi;
float *den00, **den0;
sf_file Fd0=NULL;

MPI_Comm comm=MPI_COMM_WORLD;

void gradient_init(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_d array, sf_fwi_d fwipar, bool verb1)
/*< initialize >*/
{
	int iturn, is, ir, it, sx, rx;
	float tdis;

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
	wdt=dt*interval*2.;
	dx=acpar->dx;
	dz=acpar->dz;

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
	v0=fwipar->v0;
	t0=fwipar->t0;
	wtn1=(wt1-acpar->t0)/dt+0.5;
	wtn2=(wt2-acpar->t0)/dt+0.5;
	woffn1=(woff1-acpar->r0)/acpar->dr+0.5;
	woffn2=(woff2-acpar->r0)/acpar->dr+0.5;
	weight=sf_floatalloc2(nt, nr);

	for(ir=0; ir<nr; ir++)
		for(it=0; it<nt; it++)
			weight[ir][it]=1.;
	if(nr>50) residual_weighting(weight, nt, nr, wtn1, wtn2, woffn1, woffn2, gain);

	/* reflection fwi */
	rfwi=fwipar->rfwi;
	if(rfwi){
		Fd0=sf_input("Fd0"); /* initial density model */
		den00=sf_floatalloc(nzx);
		den0=sf_floatalloc2(padnz, padnx);
		sf_floatread(den00, nzx, Fd0);
		pad2d(den00, den0, nz, nx, nb);
		free(den00);
	}

	sx=s0_v+cpuid*ds_v;
	for(ir=0; ir<nr2[cpuid]; ir++){
		rx=r0_v[cpuid]+ir*dr_v;
		tdis=sqrtf(dx*dx*(rx-sx)*(rx-sx)+dz*dz*(rz-sz)*(rz-sz))/v0;
		for(it=0; it<nt; it++){
			if(it*dt<tdis+t0) weight[r02[cpuid]+ir][it]=0.;
		}
	}
		
	ww=array->ww;
	bcxp=acpar->bcxp;
	bczp=acpar->bczp;
	bcxv=acpar->bcxv;
	bczv=acpar->bczv;
	
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
	pad2d(array->vv, vv, nz, nx, nb);
	den = sf_floatalloc2(padnz, padnx);
	pad2d(array->dd, den, nz, nx, nb);

	return;
}

void gradient_v(float *x, float *fcost, float *grad)
/*< v-den: velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **v2d, **iden; 
	float **pp, **px, **pz, **p, **vx, **vz, **term, *rr, ***wp;
	float *sendbuf, *recvbuf;

	/* residual file */
	//sf_file Fres, Fwav;
	//Fres=sf_output("Fres");
	//Fwav=sf_output("Fwav");
	//sf_putint(Fres,"n1",nt);
	//sf_putint(Fres,"n2",nr);
	//sf_putint(Fwav, "n1", padnz);
	//sf_putint(Fwav, "n2", padnx);
	//sf_putint(Fwav, "n3", (nt-1)/50+1);

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	v2d  = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);
	pp=sf_floatalloc2(nt, nr);
	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	wp=sf_floatalloc3(nz, nx, wnt);

	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz]*vv[ix][iz];
			iden[ix][iz]=1./den[ix][iz];
		}
	}

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] + v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

//			if(is==ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);

			/* output data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wp[wit][ix][iz]=p[ix+nb][iz+nb];
				wit++;
			}
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=pp[ir][it]-dd[iturn][ir][it];
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

		/* check the data residual */
		//if(is==ns/2) sf_floatwrite(pp[0], nr*nt, Fres);
		//sf_fileclose(Fres);

		/* initialization */
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);
	
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=((idt-bcxp[ix])*px[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] - iden[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] - term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) 
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=-2./v2d[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp[wit+1][ix][iz]-wp[wit-1][ix][iz])/wdt*p[ix+nb][iz+nb]*temp;
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
	//sf_warning("cpuid=%d fcost=%15.13f", cpuid, *fcost);
	
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
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*pp);  free(pp);  free(*px);   free(px);
	free(*pz);  free(pz);  free(*p);    free(p);
	free(*vx);  free(vx);  free(*vz);   free(vz);
	free(*term);  free(term);  free(rr);   
	free(**wp); free(*wp); free(wp);
}

void gradient_d(float *x, float *fcost, float *grad)
/*< v-den: density gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **v2d, **iden; 
	float **pp, **px, **pz, **p, **vx, **vz, **term, *rr, ***wp, ***wvx, ***wvz;
	float *sendbuf, *recvbuf;

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, den, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	v2d  = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);

	pp=sf_floatalloc2(nt, nr);
	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	wp=sf_floatalloc3(nz, nx, wnt);
	wvx=sf_floatalloc3(nz, nx, wnt);
	wvz=sf_floatalloc3(nz, nx, wnt);

	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz]*vv[ix][iz];
			iden[ix][iz]=1./den[ix][iz];
		}
	}

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] + v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

//			if(is==ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);

			/* output data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						wp[wit][ix][iz]=p[ix+nb][iz+nb];
						wvx[wit][ix][iz]=vx[ix+nb][iz+nb];
						wvz[wit][ix][iz]=vz[ix+nb][iz+nb];
					}
				}
				wit++;
			}
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=pp[ir][it]-dd[iturn][ir][it];
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

		/* initialization */
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);
	
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=((idt-bcxp[ix])*px[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] - iden[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] - term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) 
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=-1./v2d[ix+nb][iz+nb]/den[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p[ix+nb][iz+nb]*temp;
							grad[ix*nz+iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz[ix+nb][iz+nb];
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
	//sf_warning("cpuid=%d fcost=%15.13f", cpuid, *fcost);
	
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
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*pp);  free(pp);  free(*px);   free(px);
	free(*pz);  free(pz);  free(*p);    free(p);
	free(*vx);  free(vx);  free(*vz);   free(vz);
	free(*term); free(term); free(rr);   
	free(**wp);  free(*wp);  free(wp);
	free(**wvx); free(*wvx); free(wvx);
	free(**wvz); free(*wvz); free(wvz);
}

void gradient_vhat(float *x, float *fcost, float *grad)
/*< v-i: velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **v2d, **iden; 
	float **pp, **px, **pz, **p, **vx, **vz, **term, *rr, ***wp, ***wvx, ***wvz;
	float *sendbuf, *recvbuf;

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	v2d  = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);

	pp=sf_floatalloc2(nt, nr);
	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	wp=sf_floatalloc3(nz, nx, wnt);
	wvx=sf_floatalloc3(nz, nx, wnt);
	wvz=sf_floatalloc3(nz, nx, wnt);

	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
			iden[ix][iz]=vv[ix][iz]/den[ix][iz];
		}
	}

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] + v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

//			if(is==ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);

			/* output data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						wp[wit][ix][iz]=p[ix+nb][iz+nb];
						wvx[wit][ix][iz]=vx[ix+nb][iz+nb];
						wvz[wit][ix][iz]=vz[ix+nb][iz+nb];
					}
				}
				wit++;
			}
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=pp[ir][it]-dd[iturn][ir][it];
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

		/* initialization */
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);
	
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=((idt-bcxp[ix])*px[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] - iden[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] - term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) 
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=-1./v2d[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p[ix+nb][iz+nb]*temp;
							temp=-1./iden[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx[ix+nb][iz+nb]*temp;
							grad[ix*nz+iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz[ix+nb][iz+nb]*temp;
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
	//sf_warning("cpuid=%d fcost=%15.13f", cpuid, *fcost);
	
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
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*pp);  free(pp);  free(*px);   free(px);
	free(*pz);  free(pz);  free(*p);    free(p);
	free(*vx);  free(vx);  free(*vz);   free(vz);
	free(*term); free(term); free(rr);   
	free(**wp);  free(*wp);  free(wp);
	free(**wvx); free(*wvx); free(wvx);
	free(**wvz); free(*wvz); free(wvz);
}

void gradient_vhat0(float *x, float *fcost, float *grad)
/*< v-i: rfwi velocity gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **v2d, **iden, **v2d0, **iden0, **pp, **term, *rr; 
	float **px, **pz, **p, **vx, **vz, ***wp, ***wvx, ***wvz;
	float **px0, **pz0, **p0, **vx0, **vz0, ***wp0, ***wvx0, ***wvz0;
	float *sendbuf, *recvbuf;

	sf_file Fres;
//	sf_file Fres, Fwav0, Fwav;
	Fres=sf_output("Fres");
//	Fwav0=sf_output("Fwav0");
//	Fwav=sf_output("Fwav");
	sf_putint(Fres,"n1",nt);
	sf_putint(Fres,"n2",nr);
//	sf_putint(Fwav, "n1", padnz);
//	sf_putint(Fwav, "n2", padnx);
//	sf_putint(Fwav, "n3", (nt-1)/50+1);
//	sf_putint(Fwav0, "n1", padnz);
//	sf_putint(Fwav0, "n2", padnx);
//	sf_putint(Fwav0, "n3", (nt-1)/50+1);

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	v2d   = sf_floatalloc2(padnz, padnx);
	iden  = sf_floatalloc2(padnz, padnx);
	v2d0  = sf_floatalloc2(padnz, padnx);
	iden0 = sf_floatalloc2(padnz, padnx);
	pp=sf_floatalloc2(nt, nr);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	wp=sf_floatalloc3(nz, nx, wnt);
	wvx=sf_floatalloc3(nz, nx, wnt);
	wvz=sf_floatalloc3(nz, nx, wnt);

	px0=sf_floatalloc2(padnz, padnx);
	pz0=sf_floatalloc2(padnz, padnx);
	p0 =sf_floatalloc2(padnz, padnx);
	vx0=sf_floatalloc2(padnz, padnx);
	vz0=sf_floatalloc2(padnz, padnx);
	wp0=sf_floatalloc3(nz, nx, wnt);
	wvx0=sf_floatalloc3(nz, nx, wnt);
	wvz0=sf_floatalloc3(nz, nx, wnt);

	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
			iden[ix][iz]=vv[ix][iz]/den[ix][iz];
			v2d0[ix][iz]=den0[ix][iz]*vv[ix][iz];
			iden0[ix][iz]=vv[ix][iz]/den0[ix][iz];
		}
	}

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));

		memset(px0[0], 0., padnzx*sizeof(float));
		memset(pz0[0], 0., padnzx*sizeof(float));
		memset(vx0[0], 0., padnzx*sizeof(float));
		memset(vz0[0], 0., padnzx*sizeof(float));
		memset( p0[0], 0., padnzx*sizeof(float));

		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			// calculate u
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] + v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);

			/* output data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p[rx][rz];
			}

			// calculate u0
			/* update px0 and load source */
			derivpx(vx0, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px0[ix][iz]=( (idt-bcxp[ix])*px0[ix][iz] + v2d0[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz0 */
			derivpz(vz0, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz0[ix][iz]=((idt-bczp[iz])*pz0[ix][iz] + v2d0[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p0[ix][iz]=px0[ix][iz]+pz0[ix][iz];
				}

			/* update vx0 */
			derivvx(p0, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx0[ix][iz]=((idt-bcxv[ix])*vx0[ix][iz] + term[ix][iz]*iden0[ix][iz]) / (idt+bcxv[ix]);

			/* update vz0 */
			derivvz(p0, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz0[ix][iz]=((idt-bczv[iz])*vz0[ix][iz] + term[ix][iz]*iden0[ix][iz]) / (idt+bczv[iz]);

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						wp[wit][ix][iz]=p[ix+nb][iz+nb]-p0[ix+nb][iz+nb];
						wvx[wit][ix][iz]=vx[ix+nb][iz+nb]-vx0[ix+nb][iz+nb];
						wvz[wit][ix][iz]=vz[ix+nb][iz+nb]-vz0[ix+nb][iz+nb];
						wp0[wit][ix][iz]=p0[ix+nb][iz+nb];
						wvx0[wit][ix][iz]=vx0[ix+nb][iz+nb];
						wvz0[wit][ix][iz]=vz0[ix+nb][iz+nb];
					}
				}
				wit++;
			}
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=pp[ir][it]-dd[iturn][ir][it];
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

//		/* check the data residual */
		if(is==ns/2) sf_floatwrite(pp[0], nr*nt, Fres);
		sf_fileclose(Fres);
		
		/* initialization */
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		
		memset(px0[0], 0., padnzx*sizeof(float));
		memset(pz0[0], 0., padnzx*sizeof(float));
		memset(vx0[0], 0., padnzx*sizeof(float));
		memset(vz0[0], 0., padnzx*sizeof(float));

		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);
	
			// calculate lambda
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=((idt-bcxp[ix])*px[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

		//	if(is==ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] - iden[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] - term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);
		
			// calculate lambda0
			/* update px0 and load source */
			derivpx(vx0, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px0[ix][iz]=((idt-bcxp[ix])*px0[ix][iz] - v2d0[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz0 */
			derivpz(vz0, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz0[ix][iz]=((idt-bczp[iz])*pz0[ix][iz] - v2d0[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p0[ix][iz]=px0[ix][iz]+pz0[ix][iz];
				}

	//		if(is==ns/2 && it%50==0) sf_floatwrite(p0[0], padnzx, Fwav0);

			/* update vx0 */
			derivvx(p0, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx0[ix][iz]=((idt-bcxv[ix])*vx0[ix][iz] - iden0[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz0 */
			derivvz(p0, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz0[ix][iz]=((idt-bczv[iz])*vz0[ix][iz] - iden0[ix][iz]*term[ix][iz]) / (idt+bczv[iz]);
		
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) 
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=-1./v2d[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p0[ix+nb][iz+nb]*temp;
							temp=-1./iden[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx0[ix+nb][iz+nb]*temp;
							grad[ix*nz+iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz0[ix+nb][iz+nb]*temp;

							temp=-1./v2d0[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp0 [wit+1][ix][iz]-wp0 [wit-1][ix][iz])/wdt*(p[ix+nb][iz+nb]-p0[ix+nb][iz+nb])*temp;
							temp=-1./iden0[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wvx0[wit+1][ix][iz]-wvx0[wit-1][ix][iz])/wdt*(vx[ix+nb][iz+nb]-vx0[ix+nb][iz+nb])*temp;
							grad[ix*nz+iz] += (wvz0[wit+1][ix][iz]-wvz0[wit-1][ix][iz])/wdt*(vz[ix+nb][iz+nb]-vz0[ix+nb][iz+nb])*temp;
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
	//sf_warning("cpuid=%d fcost=%15.13f", cpuid, *fcost);
	
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
		scaling=0.05/dmax;
		first=false;
	}

	/* smooth gradient */
	gradient_smooth2(grectx, grectz, nx, nz, waterz, scaling, grad);

	/* free allocated memory */
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*v2d0); free(v2d0); free(*iden0); free(iden0);
	free(*pp);  free(pp); free(*term); free(term); free(rr);   

	free(*px);   free(px); free(*pz);  free(pz);  
	free(*p);    free(p); 
	free(*vx);  free(vx); free(*vz);   free(vz);
	free(**wp);  free(*wp);  free(wp);
	free(**wvx); free(*wvx); free(wvx);
	free(**wvz); free(*wvz); free(wvz);

	free(*px0);   free(px0); free(*pz0); free(pz0);  
	free(*p0);    free(p0); 
	free(*vx0);  free(vx0); free(*vz0); free(vz0);
	free(**wp0);  free(*wp0);  free(wp0);
	free(**wvx0); free(*wvx0); free(wvx0);
	free(**wvz0); free(*wvz0); free(wvz0);
}

void gradient_i(float *x, float *fcost, float *grad)
/*< v-i: impedance gradient >*/
{
	int ix, iz, is, ir, it, wit, iturn;
	int sx, rx;

	float temp, dmax;
	float **v2d, **iden; 
	float **pp, **px, **pz, **p, **vx, **vz, **term, *rr, ***wp, ***wvx, ***wvz;
	float *sendbuf, *recvbuf;

	sf_file Fres;
	Fres=sf_output("Fres");
	sf_putint(Fres,"n1",nt);
	sf_putint(Fres,"n2",nr);
	
	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, den, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));

	v2d  = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);

	pp=sf_floatalloc2(nt, nr);
	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	wp=sf_floatalloc3(nz, nx, wnt);
	wvx=sf_floatalloc3(nz, nx, wnt);
	wvz=sf_floatalloc3(nz, nx, wnt);

	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
			iden[ix][iz]=vv[ix][iz]/den[ix][iz];
		}
	}

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(pp[0], 0., nr*nt*sizeof(float));
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*ww[it]) ) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] + v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

//			if(is==ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] + term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);

			/* output data */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				pp[r02[is]+ir][it]=p[rx][rz];
			}

			/* save wavefield */
			if(it%interval==0){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						wp[wit][ix][iz]=p[ix+nb][iz+nb];
						wvx[wit][ix][iz]=vx[ix+nb][iz+nb];
						wvz[wit][ix][iz]=vz[ix+nb][iz+nb];
					}
				}
				wit++;
			}
		} // end of time loop

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		wit--;
		
		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=pp[ir][it]-dd[iturn][ir][it];
				pp[ir][it] *= weight[ir][it];
				*fcost += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		iturn++;

		/* check the data residual */
		if(is==ns/2) sf_floatwrite(pp[0], nr*nt, Fres);
		sf_fileclose(Fres);

		/* initialization */
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			if(verb) sf_warning("Backward propagation is=%d; it=%d;", is+1, it);
	
			/* update px and load source */
			derivpx(vx, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr2[is]; ir++){
				rx=r0_v[is]+ir*dr_v;
				term[rx][rz] += pp[r02[is]+ir][it];
			}

#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					px[ix][iz]=((idt-bcxp[ix])*px[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bcxp[ix]);

			/* update pz */
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					pz[ix][iz]=((idt-bczp[iz])*pz[ix][iz] - v2d[ix][iz]*term[ix][iz]) / (idt+bczp[iz]);
					p[ix][iz]=px[ix][iz]+pz[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=((idt-bcxv[ix])*vx[ix][iz] - iden[ix][iz]*term[ix][iz]) / (idt+bcxv[ix]);

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
			private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=((idt-bczv[iz])*vz[ix][iz] - term[ix][iz]*iden[ix][iz]) / (idt+bczv[iz]);
			
			/* calculate gradient  */
			if(it%interval==0){
				if(wit != wnt-1 && wit != 0){ // avoid the first and last time step
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) 
#endif
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							temp=-1./v2d[ix+nb][iz+nb]/den[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p[ix+nb][iz+nb]*temp;
							temp=1./vv[ix+nb][iz+nb];
							grad[ix*nz+iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx[ix+nb][iz+nb]*temp;
							grad[ix*nz+iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz[ix+nb][iz+nb]*temp;
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* misfit reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, comm);
	//sf_warning("cpuid=%d fcost=%15.13f", cpuid, *fcost);
	
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
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*pp);  free(pp);  free(*px);   free(px);
	free(*pz);  free(pz);  free(*p);    free(p);
	free(*vx);  free(vx);  free(*vz);   free(vz);
	free(*term); free(term); free(rr);   
	free(**wp);  free(*wp);  free(wp);
	free(**wvx); free(*wvx); free(wvx);
	free(**wvz); free(*wvz); free(wvz);
}

void fwi(sf_file Fdat, sf_file Finv, sf_file Ferr, sf_file Fmod, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_d array, sf_fwi_d fwipar, sf_optim optpar, int para_type, bool verb1)
/*< fwi >*/
{
	int iter=0, flag;
	float fcost, threshold[2];
	float *x, *direction, *grad;
	sf_gradient gradient;
	FILE *fp;

	/* initialize */
	gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb1);

	/* gradient type */
	if(para_type==1){
		if(fwipar->grad_type==1){
			gradient=gradient_v;
			x=array->vv;
		}else{
			gradient=gradient_d;
			x=array->dd;
;
		}
	}else{
		if(fwipar->grad_type==1){
			if(rfwi) gradient=gradient_vhat0;
			else gradient=gradient_vhat;
			x=array->vv;
		}else{
			gradient=gradient_i;
			x=array->dd;
		}
	}

	/* calculate first gradient */
	grad=sf_floatalloc(nzx);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nzx, Fgrad);

	/* if onlygrad=y, program terminates */
	if(fwipar->onlygrad) return; 

	/* hard thresholding */
	if(fwipar->grad_type==1){
		threshold[0]=fwipar->v1;
		threshold[1]=fwipar->v2;
	}else{
		threshold[0]=fwipar->den1;
		threshold[1]=fwipar->den2;
	}

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");

	direction=sf_floatalloc(nzx);
	optpar->sk=sf_floatalloc2(nzx, optpar->npair);
	optpar->yk=sf_floatalloc2(nzx, optpar->npair);

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
	optpar->err[0]=optpar->fk;
	if(mpipar->cpuid==0) sf_floatwrite(x, nzx, Fmod);

	iter=0;
	if(mpipar->cpuid==0){
		l2norm(nzx, grad, &optpar->gk_norm);
		print_iteration(fp, iter, optpar);
	}

	/* optimization loop */
	for(iter=0; iter<optpar->niter; iter++){
		if(mpipar->cpuid==0) sf_warning("-------------------iter=%d----------------------", iter+1);
		
		optpar->ils=0;
		if(iter%optpar->repeat==0) optpar->alpha=1.;

		/* search direction */
		if(iter==0){
			reverse(nzx, grad, direction);
		}else{
			lbfgs_update(nzx, x, grad, optpar->sk, optpar->yk, optpar);
			lbfgs_direction(nzx, grad, direction, optpar->sk, optpar->yk, optpar);
		}

		/* line search */
		lbfgs_save(nzx, x, grad, optpar->sk, optpar->yk, optpar);
		line_search(nzx, x, grad, direction, gradient, optpar, threshold, &flag, mpipar->cpuid, 1);
		optpar->err[iter+1]=optpar->fk;
		if(mpipar->cpuid==0) sf_floatwrite(x, nzx, Fmod);
		
		if(mpipar->cpuid==0){
			l2norm(nzx, grad, &optpar->gk_norm);
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
	} // end of iteration

	if(mpipar->cpuid==0 && iter==optpar->niter){
		fprintf(fp, "Maximum Iteration Number Reached\n");
	}

	/* output vel & misfit */
	if(mpipar->cpuid==0) sf_floatwrite(x, nzx, Finv);
	if(mpipar->cpuid==0) sf_floatwrite(optpar->err, optpar->niter+1, Ferr);
	if(mpipar->cpuid==0) fclose(fp);

	return;
}
