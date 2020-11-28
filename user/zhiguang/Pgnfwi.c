/* Acoustic FWI with Gauss-Newton optimization */
/*
 Copyright (C) 2017 The University of Texas at Austin
 
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
int frectx, frectz; // wavefield
int waterz, wtn1, wtn2, woffn1, woffn2, grectx, grectz; // gradient
int drectx, drectz, nrepeat, ider=0; // smoothing kernel
int lniter, ITER=0; // GN CG iteration number

float dt, dt2, dx2, dz2; // wavefield
float wt1, wt2, woff1, woff2, gain, scaling, swap; // gradient; swap recording misfit without smoothing
float ***dd, **vv, *ww, *bc, **weight, threshold[2]; // arrays
float **p0, **p1, **p2, **term, **tmparray, *rr, **pp, ***wave; // multiple use

MPI_Comm comm=MPI_COMM_WORLD;

void gradient_init(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_s array, sf_fwi_s fwipar, bool verb1)
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
	
	dt=acpar->dt;
	dt2=dt*dt;
	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;

	/* smoothing kernel parameters */
	drectx=fwipar->drectx;
	drectz=fwipar->drectz;
	nrepeat=fwipar->nrepeat;

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
	pad2d(array->vv, vv, nz, nx, nb);

	/* hard thresholding */
	threshold[0]=fwipar->v1;
	threshold[1]=fwipar->v2;

	/* CG iteration number */
	lniter=fwipar->lniter;

	/* storage for forward wavefield */
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	pp=sf_floatalloc2(nt, nr);
	wave=sf_floatalloc3(nz, nx, nt);

	return;
}

void gradient_standard(float *x, float *fcost, float *grad)
/*< standard velocity gradient >*/
{
	int ix, iz, is, ir, it, iturn;
	int sx, rx;

	float temp, dmax;
	float *sendbuf, *recvbuf;

	/* residual file */
	sf_file Fres;
	Fres=sf_output("Fres");
	sf_putint(Fres,"n1",nt);
	sf_putint(Fres,"n2",nr);

	/* initialize fcost */
	*fcost=0.;
	/* update velocity */
	pad2d(x, vv, nz, nx, nb);
	/* initialize gradient */
	memset(grad, 0., nzx*sizeof(float));
	/* initialize data misfit */
	swap=0.;

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(pp[0], 0., nr*nt*sizeof(float));
		
		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

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
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}

			/* save weighted wavefield  */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz,temp) \
			shared(nx,nz,vv,wave,p1,it,dt2)
#endif
			for(ix=nb; ix<nx+nb; ix++){
				for(iz=nb; iz<nz+nb; iz++){
					temp=vv[ix][iz];
					temp=temp*temp*temp;
					temp=-2./temp;
					wave[it][ix-nb][iz-nb]=(p2[ix][iz]-2.*p1[ix][iz]+p0[ix][iz])/dt2*temp;
				}
			}

			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
		} // end of time loop

		/* calculate data residual and data misfit */
		for(ir=0; ir<nr; ir++){
			for(it=0; it<nt; it++){
				pp[ir][it]=dd[iturn][ir][it]-pp[ir][it];
				pp[ir][it] *= weight[ir][it];
				swap += 0.5*pp[ir][it]*pp[ir][it];
			}
		}
		smooth_misfit(pp, fcost, nr, nt, drectx, drectz, nrepeat, ider);

		iturn++;

		/* check the data residual */
		if(is==ns/2) sf_floatwrite(pp[0], nr*nt, Fres);
		sf_fileclose(Fres);

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
			
			/* load data residual */
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
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* calculate gradient  */
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) \
				shared(nx,nz,wave,p1,grad)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						grad[ix*nz+iz] += wave[it][ix][iz]*p1[ix+nb][iz+nb];
					}
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
	
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=&swap;
	}else{
		sendbuf=&swap;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(&swap, 1, MPI_FLOAT, 0, comm);

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
	//free(*p0); free(p0); free(*p1); free(p1);
	//free(*p2); free(p2); free(*pp); free(pp);
	//free(**wave); free(*wave); free(wave);
	//free(*rr); free(rr); free(*term); free(term);
}

void Hessian_vector_mul(float *p, float *Hp)
/*< Hessian-vector multiplication >*/
{
	int ix, iz, is, ir, it, iturn, rx;
	float *sendbuf, *recvbuf, temp;

	/* initialize Hp */
	memset(Hp, 0., nzx*sizeof(float));

	iturn=0;
	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(pp[0], 0., nr*nt*sizeof(float));
		
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

			/* load source */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(term,nx,nz,it)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					term[ix+nb][iz+nb] += wave[it][ix][iz]*p[ix*nz+iz];
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

		/* smooth data */
		temp=0.;
		smooth_misfit(pp, &temp, nr, nt, drectx, drectz, nrepeat, 0);

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
			
			/* load data residual */
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
			for(ix=4; ix<padnx-4; ix++){
				for(iz=4; iz<padnz-4; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* calculate Hessian-vector product  */
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz) \
				shared(nx,nz,wave,p1,Hp)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						Hp[ix*nz+iz] += wave[it][ix][iz]*p1[ix+nb][iz+nb];
					}
				}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, bc, padnx, padnz, nb);
			apply_sponge(p1, bc, padnx, padnz, nb);
		} // end of time loop
	}// end of shot loop
	MPI_Barrier(comm);
	
	/* Hp reduction */
	//MPI_ALLreduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, comm);
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=Hp;
	}else{
		sendbuf=Hp;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(Hp, nzx, MPI_FLOAT, 0, comm);

	return;
}

void gauss_newton(float *grad, float *direction)
/*< Gauss Newton equation with CG solver >*/
{
	int iter;
	float r2, r3, pHp, alpha, beta;
	float *r, *p, *Hp;
	FILE *fp=NULL;

	/* initialize solution */
	memset(direction, 0., nzx*sizeof(float));

	/* storage allocation */
	r=sf_floatalloc(nzx);
	p=sf_floatalloc(nzx);
	Hp=sf_floatalloc(nzx);

	/* initialize r and p */
	copy(nzx, grad, r);
	reverse(nzx, r, p);
	dot_product(nzx, r, r, &r2);

	/* open the printing file */
	if(cpuid==0){
		fp=fopen("cgiterate.txt","a");
		fprintf(fp, "****************************\n");
		fprintf(fp, "Iteration: %d\n", ITER);
		fclose(fp); /* get written to disk right away */
		fp=fopen("cgiterate.txt","a");
	}
	ITER++;

	for (iter=0; iter<lniter; iter++){

		if(cpuid==0){
			fprintf(fp, "    liter=%3d   misfit=%3.2e\n", iter, sqrtf(r2));
			fclose(fp); /* get written to disk right away */
			fp=fopen("cgiterate.txt","a");
		}

		Hessian_vector_mul(p, Hp);
		dot_product(nzx, p, Hp, &pHp);
		alpha=r2/pHp;

		combine_vectors(nzx, direction, p, 1., alpha, direction);
		combine_vectors(nzx, r, Hp, 1., alpha, r);

		r3=r2;
		dot_product(nzx, r, r, &r2);
		beta=r2/r3;

		combine_vectors(nzx, r, p, -1, beta, p);
	}

	if(cpuid==0){
		fprintf(fp, "    liter=%3d   misfit=%3.2e\n", iter, sqrtf(r2));
		fclose(fp); /* get written to disk right away */
	}

	free(r); free(p); free(Hp);
}

void gnfwi(sf_file Fdat, sf_file Finv, sf_file Ferr, sf_file Fmod, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_s array, sf_fwi_s fwipar, sf_optim optpar, bool verb1)
/*< fwi >*/
{
	float fcost, *x, *grad, *direction;
	sf_gradient gradient;

	/* initialize */
	gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb1);

	/* gradient type */
	gradient=gradient_standard;
	x=array->vv;
	
	/* calculate first gradient */
	grad=sf_floatalloc(nzx);
	gradient(x, &fcost, grad);

	/* solve for search direction */
	direction=sf_floatalloc(nzx);
	gauss_newton(grad, direction);

	/* output first gradient and GN update */
	if(cpuid==0) sf_floatwrite(grad, nzx, Fgrad);
	if(cpuid==0) sf_floatwrite(direction, nzx, Fgrad);

	/* if onlygrad=y, program terminates */
	if(fwipar->onlygrad) return; 

//	if(cpuid==0) fp=fopen("iterate.txt","a");
//
//	optpar->igrad=1;
//	optpar->ils=0;
//	optpar->f0=fcost;
//	optpar->fk=fcost;
//	optpar->alpha=1.;
//
//	/* initialize data error vector */
//	for(iter=0; iter<optpar->nerr; iter++){
//		optpar->err[iter]=0.;
//	}
//	optpar->err[0]=optpar->fk;
//	if (optpar->err_type==1) optpar->err[optpar->nerr/2]=swap;
//	if(cpuid==0) sf_floatwrite(x, nzx, Fmod);
//
//	iter=0;
//	if(cpuid==0){
//		l2norm(nzx, grad, &optpar->gk_norm);
//		print_iteration2(fp, iter, optpar);
//	}
//
//	/* optimization loop */
//	for(iter=0; iter<optpar->niter; iter++){
//		if(cpuid==0) sf_warning("-------------------iter=%d----------------------", iter+1);
//		
//		optpar->ils=0;
//		if(iter%optpar->repeat==0) optpar->alpha=1.;
//
//		/* search direction */
//		if(iter==0){
//			reverse(nzx, grad, direction);
//		}else{
//			lbfgs_update(nzx, x, grad, optpar->sk, optpar->yk, optpar);
//			lbfgs_direction(nzx, grad, direction, optpar->sk, optpar->yk, optpar);
//		}
//
//		/* line search */
//		lbfgs_save(nzx, x, grad, optpar->sk, optpar->yk, optpar);
//		line_search(nzx, x, grad, direction, gradient, optpar, threshold, &flag, cpuid, 1);
//		optpar->err[iter+1]=optpar->fk;
//		if(optpar->err_type==1) optpar->err[optpar->nerr/2+iter+1]=swap;
//		if(cpuid==0) sf_floatwrite(x, nzx, Fmod);
//		
//		if(cpuid==0){
//			l2norm(nzx, grad, &optpar->gk_norm);
//			print_iteration(fp, iter+1, optpar);
//			fclose(fp); /* get written to disk right away */
//			fp=fopen("iterate.txt","a");
//		}
//
//		if(cpuid==0 && flag==2){
//			fprintf(fp, "Line Search Failed\n");
//			break;
//		}
//
//		if(cpuid==0 && optpar->fk/optpar->f0 < optpar->conv_error){
//			fprintf(fp, "Convergence Criterion Reached\n");
//			break;
//		}
//	} // end of iteration
//
//	if(cpuid==0 && iter==optpar->niter){
//		fprintf(fp, "Maximum Iteration Number Reached\n");
//	}
//
//	/* output vel & misfit */
//	if(cpuid==0) sf_floatwrite(x, nzx, Finv);
//	if(cpuid==0) sf_floatwrite(optpar->err, optpar->nerr, Ferr);
//	if(cpuid==0) fclose(fp);
//
//	/* calculate prediction-corrector */
//	if(fwipar->tangent==1){
//		Ftangent=sf_output("Ftangent");
//		sf_putint(Ftangent, "n1", nz);
//		sf_putint(Ftangent, "n2", nx);
//		sf_putint(Ftangent, "n3", 1);
//
//		direction2=sf_floatalloc(nzx);
//		lbfgs_update(nzx, x, grad, optpar->sk, optpar->yk, optpar);
//
//		/* update sigma 1 */
//		ider=1;
//		gradient(x, &fcost, grad);
//		lbfgs_direction(nzx, grad, direction, optpar->sk, optpar->yk, optpar);
//		for(i=0; i<nzx; i++) direction2[i] = direction[i]*fwipar->sigma1;
//
//		/* update sigma 2 */
//		ider=2;
//		gradient(x, &fcost, grad);
//		lbfgs_direction(nzx, grad, direction, optpar->sk, optpar->yk, optpar);
//		for(i=0; i<nzx; i++) direction2[i] = direction2[i]+direction[i]*fwipar->sigma2;
//
//		if(cpuid==0) sf_floatwrite(direction2, nzx, Ftangent);
//	}

	return;
}
