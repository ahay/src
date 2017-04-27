/* Forward modeling */
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
/*^*/

void forward_modeling(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_s array, sf_encode encodepar, bool verb)
/*< acoustic forward modeling >*/
{
	int ix, iz, is, ir, it, isou;
	int sx, rx, sz, rz, frectx, frectz;
	int nz, nx, padnz, padnx, padnzx, nt, nr, nb;

	float dx2, dz2, dt2, dt;
	float **vv, **dd;
	float **p0, **p1, **p2, **term, **tmparray, **rr;

	FILE *swap;

	MPI_Comm comm=MPI_COMM_WORLD;

	swap=fopen("temswap.bin", "wb+");

	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nz=acpar->nz;
	nx=acpar->nx;
	nt=acpar->nt;
	nr=acpar->nr;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	frectx=soupar->frectx;
	frectz=soupar->frectz;

	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;
	dt2=acpar->dt*acpar->dt;
	dt=acpar->dt;

	vv = sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);

	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc2(padnzx, soupar->nsource);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(dd[0], 0., nr*nt*sizeof(float));
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map2(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		for(it=0; it<nt; it++){
			if(verb) sf_warning("Modeling is=%d; it=%d;", is+1, it);

			/* output data */
			for(ir=0; ir<acpar->nr2[is]; ir++){
				rx=acpar->r0_v[is]+ir*acpar->dr_v;
				dd[acpar->r02[is]+ir][it]=p1[rx][rz];
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* load source */
			for (isou=0; isou<soupar->nsource; isou++){
				if(encodepar==NULL){
					for(ix=0; ix<padnx; ix++){
						for(iz=0; iz<padnz; iz++){
							term[ix][iz] += rr[isou][ix*padnz+iz]*array->ww[it];
						}
					}
				}else{
					if(it>encodepar->shift[is][isou]){
						for(ix=0; ix<padnx; ix++){
							for(iz=0; iz<padnz; iz++){
								term[ix][iz] += rr[isou][ix*padnz+iz]*array->ww[it-encodepar->shift[is][isou]]*encodepar->sign[is][isou];
							} // iz
						} // ix
					} // it >shift
				} // encodepar==NULL
			} // end of isou

			/* update */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}
			
			/* swap wavefield pointer of different time steps */
			tmparray=p0; p0=p1; p1=p2; p2=tmparray;

			/* boundary condition */
			apply_sponge(p0, acpar->bc, padnx, padnz, nb);
			apply_sponge(p1, acpar->bc, padnx, padnz, nb);
		} // end of time loop

		fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
		fwrite(dd[0], sizeof(float), nr*nt, swap);
	}// end of shot loop
	fclose(swap);
	MPI_Barrier(comm);

	/* transfer data to Fdat */
	if(mpipar->cpuid==0){
		swap=fopen("temswap.bin", "rb");
		for(is=0; is<acpar->ns; is++){
			fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
			fread(dd[0], sizeof(float), nr*nt, swap);
			sf_floatwrite(dd[0], nr*nt, Fdat);
		}
		fclose(swap);
		remove("temswap.bin");
	}
	MPI_Barrier(comm);
	
	/* release allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*vv); free(vv);
	free(*dd); free(dd);
	free(*rr); free(rr); free(*term); free(term);
}

void forward_modeling_q(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_q array, bool verb)
/*< visco-acoustic forward modeling >*/
{
	int ix, iz, is, ir, it;
	int sx, rx, sz, rz, frectx, frectz;
	int nz, nx, padnz, padnx, padnzx, nt, nr, nb;

	float dx2, dz2, dt2, idt;
	float **vv, **tau, **taus, **dd;
	float **p0, **p1, **p2, **r1, **r2, **term, **tmparray, *rr;

	FILE *swap;

	MPI_Comm comm=MPI_COMM_WORLD;

	swap=fopen("temswap.bin", "wb+");

	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nz=acpar->nz;
	nx=acpar->nx;
	nt=acpar->nt;
	nr=acpar->nr;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	frectx=soupar->frectx;
	frectz=soupar->frectz;

	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;
	dt2=acpar->dt*acpar->dt;
	idt=1./acpar->dt;

	vv = sf_floatalloc2(padnz, padnx);
	tau= sf_floatalloc2(padnz, padnx);
	taus=sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);

	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	r1=sf_floatalloc2(padnz, padnx);
	r2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(dd[0], 0., nr*nt*sizeof(float));
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		for(it=0; it<nt; it++){
			if(verb) sf_warning("Modeling is=%d; it=%d;", is+1, it);

			/* output data */
			for(ir=0; ir<acpar->nr2[is]; ir++){
				rx=acpar->r0_v[is]+ir*acpar->dr_v;
				dd[acpar->r02[is]+ir][it]=p1[rx][rz];
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);

			/* calculate r, load source and update wavefield */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					r2[ix][iz]=
						(tau[ix][iz]/taus[ix][iz]*term[ix][iz]
						 + (idt-0.5/taus[ix][iz])*r1[ix][iz])
						/(idt+0.5/taus[ix][iz]);
					term[ix][iz]=term[ix][iz]*(1.+tau[ix][iz]) - (r2[ix][iz]+r1[ix][iz])*0.5 + rr[ix*padnz+iz]*array->ww[it];
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

		fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
		fwrite(dd[0], sizeof(float), nr*nt, swap);
	}// end of shot loop
	fclose(swap);
	MPI_Barrier(comm);

	/* transfer data to Fdat */
	if(mpipar->cpuid==0){
		swap=fopen("temswap.bin", "rb");
		for(is=0; is<acpar->ns; is++){
			fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
			fread(dd[0], sizeof(float), nr*nt, swap);
			sf_floatwrite(dd[0], nr*nt, Fdat);
		}
		fclose(swap);
		remove("temswap.bin");
	}
	MPI_Barrier(comm);
	
	/* release allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*r1); free(r1);
	free(*r2); free(r2); free(*vv); free(vv);
	free(*tau); free(tau); free(*taus); free(taus);
	free(*dd); free(dd); free(rr); 
	free(*term); free(term);
}

void forward_modeling_d(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_d array, int para_type, bool verb)
/*< Variable-density acoustic forward modeling >*/
{
	int ix, iz, is, ir, it;
	int sx, rx, sz, rz, frectx, frectz;
	int nz, nx, padnz, padnx, padnzx, nt, nr, ns, nb;
	int numprocs, cpuid, iturn, nspad;

	float dx, dz, dt, idt;
	float **vv, **den, **v2d, **iden, **dd, *bcxp, *bczp, *bcxv, *bczv;
	float **px, **pz, **p, **vx, **vz, **term, *rr, *sendbuf, *recvbuf, ***ddall;

	//sf_file Fwav;
	MPI_Comm comm=MPI_COMM_WORLD;

	//Fwav=sf_output("Fwav");

	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nz=acpar->nz;
	nx=acpar->nx;
	nt=acpar->nt;
	nr=acpar->nr;
	ns=acpar->ns;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	frectx=soupar->frectx;
	frectz=soupar->frectz;

	dx=acpar->dx;
	dz=acpar->dz;
	dt=acpar->dt;
	idt=1./acpar->dt;

	numprocs=mpipar->numprocs;
	cpuid=mpipar->cpuid;
	if(ns%numprocs==0) nspad=ns;
	else nspad=(ns/numprocs+1)*numprocs;
	if(cpuid==0) ddall=sf_floatalloc3(nt, nr, nspad);

	//sf_putint(Fwav, "n1", padnz);
	//sf_putint(Fwav, "n2", padnx);
	//sf_putint(Fwav, "n3", (nt-1)/50+1);

	vv = sf_floatalloc2(padnz, padnx);
	den = sf_floatalloc2(padnz, padnx);
	v2d = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);

	bcxp=acpar->bcxp;
	bczp=acpar->bczp;
	bcxv=acpar->bcxv;
	bczv=acpar->bczv;

	px=sf_floatalloc2(padnz, padnx);
	pz=sf_floatalloc2(padnz, padnx);
	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->dd, den, nz, nx, nb);

	if(para_type==1){
		for (ix=0; ix<padnx; ix++)
			for (iz=0; iz<padnz; iz++){
				v2d[ix][iz]=den[ix][iz]*vv[ix][iz]*vv[ix][iz];
				iden[ix][iz]=1./den[ix][iz];
			}
	}else{
		for (ix=0; ix<padnx; ix++)
			for (iz=0; iz<padnz; iz++){
				v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
				iden[ix][iz]=vv[ix][iz]/den[ix][iz];
			}
	}

	if(cpuid==0) sf_warning("nspad=%d, numprocs=%d, ns=%d", nspad, numprocs, ns);

	for(iturn=0; iturn*numprocs<nspad; iturn++){
		is=iturn*numprocs+cpuid;
		sf_warning("###### is=%d ######", is+1);

		memset(dd[0], 0., nr*nt*sizeof(float));
		memset(px[0], 0., padnzx*sizeof(float));
		memset(pz[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));

		if(is<ns){
			sx=acpar->s0_v+is*acpar->ds_v;
			source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

			for(it=0; it<nt; it++){
				if(verb) sf_warning("Modeling is=%d; it=%d;", is+1, it);

				/* update px and load source */
				derivpx(vx, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
				private(ix,iz) 
#endif
				for (ix=4; ix<padnx-4; ix++)
					for (iz=4; iz<padnz-4; iz++)
						px[ix][iz]=( (idt-bcxp[ix])*px[ix][iz] + v2d[ix][iz]*(term[ix][iz]+rr[ix*padnz+iz]*array->ww[it]) ) / (idt+bcxp[ix]);

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
				for(ir=0; ir<acpar->nr2[is]; ir++){
					rx=acpar->r0_v[is]+ir*acpar->dr_v;
					dd[acpar->r02[is]+ir][it]=p[rx][rz];
				}

				//	if(is==acpar->ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			} // end of time loop
		} // if is<ns

		if(cpuid==0){
			sendbuf=dd[0];
			recvbuf=ddall[iturn*numprocs][0];
		}else{
			sendbuf=dd[0];
			recvbuf=NULL;
		}
		MPI_Gather(sendbuf, nt*nr, MPI_FLOAT, recvbuf, nt*nr, MPI_FLOAT, 0, comm);
	}// end of shot loop
	MPI_Barrier(comm);
	//sf_fileclose(Fwav);

	/* transfer data to Fdat */
	if(cpuid==0) sf_floatwrite(ddall[0][0], nr*nt*ns, Fdat);
	MPI_Barrier(comm);
	
	/* release allocated memory */
	free(*px); free(px); free(*pz); free(pz);
	free(*p); free(p); free(*vx); free(vx); 
	free(*vz); free(vz); free(rr);
	free(*vv); free(vv); free(*den); free(den);
	free(*dd); free(dd); free(*term); free(term);
	free(*v2d); free(v2d); free(*iden); free(iden);
	if(cpuid==0) free(**ddall); free(*ddall); free(ddall);
}

void forward_modeling_dq(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec_dq array, int para_type, bool verb)
/*< Variable-density visco-acoustic forward modeling >*/
{
	int ix, iz, is, ir, it;
	int sx, rx, sz, rz, frectx, frectz;
	int nz, nx, padnz, padnx, padnzx, nt, nr, nb;

	float dx, dz, dt, idt;
	float **vv, **den, **v2d, **iden, **tau, **taus, **dd;
	float **p, **vx, **vz, **r1, **r2, **term, **term1, *rr, *ww, *bc;

	FILE *swap;
	sf_file Fwav;

	MPI_Comm comm=MPI_COMM_WORLD;

	swap=fopen("temswap.bin", "wb+");
	Fwav=sf_output("Fwav");

	padnz=acpar->padnz;
	padnx=acpar->padnx;
	padnzx=padnz*padnx;
	nz=acpar->nz;
	nx=acpar->nx;
	nt=acpar->nt;
	nr=acpar->nr;
	nb=acpar->nb;
	sz=acpar->sz;
	rz=acpar->rz;
	frectx=soupar->frectx;
	frectz=soupar->frectz;

	dx=acpar->dx;
	dz=acpar->dz;
	dt=acpar->dt;
	idt=1./dt;

	sf_putint(Fwav, "n1", padnz);
	sf_putint(Fwav, "n2", padnx);
	sf_putint(Fwav, "n3", (nt-1)/50+1);

	vv = sf_floatalloc2(padnz, padnx);
	den = sf_floatalloc2(padnz, padnx);
	v2d = sf_floatalloc2(padnz, padnx);
	iden = sf_floatalloc2(padnz, padnx);
	tau = sf_floatalloc2(padnz, padnx);
	taus= sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);

	p =sf_floatalloc2(padnz, padnx);
	vx=sf_floatalloc2(padnz, padnx);
	vz=sf_floatalloc2(padnz, padnx);
	r1=sf_floatalloc2(padnz, padnx);
	r2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	term1=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);
	bc=acpar->bc;
	ww=array->ww;

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);
	pad2d(array->dd, den, nz, nx, nb);
	pad2d(array->tau, tau, nz, nx, nb);
	pad2d(array->taus, taus, nz, nx, nb);

	if(para_type==1){
		for (ix=0; ix<padnx; ix++)
			for (iz=0; iz<padnz; iz++){
				v2d[ix][iz]=den[ix][iz]*vv[ix][iz]*vv[ix][iz];
				iden[ix][iz]=1./den[ix][iz];
			}
	}else{
		for (ix=0; ix<padnx; ix++)
			for (iz=0; iz<padnz; iz++){
				v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
				iden[ix][iz]=vv[ix][iz]/den[ix][iz];
			}
	}

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(dd[0], 0., nr*nt*sizeof(float));
		memset( p[0], 0., padnzx*sizeof(float));
		memset(vx[0], 0., padnzx*sizeof(float));
		memset(vz[0], 0., padnzx*sizeof(float));
		memset(r1[0], 0., padnzx*sizeof(float));
		memset(r2[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		memset(term1[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		for(it=0; it<nt; it++){
			if(verb) sf_warning("Modeling is=%d; it=%d;", is+1, it);

			/* get term and update r */
			derivpx(vx, term1, padnx, padnz, dx);
			derivpz(vz, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					term[ix][iz] += term1[ix][iz];
					r2[ix][iz] = ( (taus[ix][iz]*idt-0.5)*r1[ix][iz] - v2d[ix][iz]*tau[ix][iz]*term[ix][iz] )/ (taus[ix][iz]*idt+0.5);
				}

			/* update p and pass r */
#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++){
					p[ix][iz]= p[ix][iz] + dt*( v2d[ix][iz]*(1.+tau[ix][iz]) * (rr[ix*padnz+iz]*ww[it]-term[ix][iz]) - (r1[ix][iz]+r2[ix][iz])*0.5 );
					r1[ix][iz]=r2[ix][iz];
				}

			/* update vx */
			derivvx(p, term, padnx, padnz, dx);
#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vx[ix][iz]=vx[ix][iz] - dt*iden[ix][iz]*term[ix][iz];

			/* update vz */
			derivvz(p, term, padnx, padnz, dz);
#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) 
#endif
			for (ix=4; ix<padnx-4; ix++)
				for (iz=4; iz<padnz-4; iz++)
					vz[ix][iz]=vz[ix][iz] - dt*iden[ix][iz]*term[ix][iz];

			/* output data */
			for(ir=0; ir<acpar->nr2[is]; ir++){
				rx=acpar->r0_v[is]+ir*acpar->dr_v;
				dd[acpar->r02[is]+ir][it]=p[rx][rz];
			}

			if(is==acpar->ns/2 && it%50==0) sf_floatwrite(p[0], padnzx, Fwav);

			/* boundary condition */
			apply_sponge(p,  acpar->bc, padnx, padnz, nb);
			apply_sponge(vx, acpar->bc, padnx, padnz, nb);
			apply_sponge(vz, acpar->bc, padnx, padnz, nb);
			apply_sponge(r1, acpar->bc, padnx, padnz, nb);
		} // end of time loop

		fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
		fwrite(dd[0], sizeof(float), nr*nt, swap);
	}// end of shot loop
	fclose(swap);
	MPI_Barrier(comm);

	sf_fileclose(Fwav);

	/* transfer data to Fdat */
	if(mpipar->cpuid==0){
		swap=fopen("temswap.bin", "rb");
		for(is=0; is<acpar->ns; is++){
			fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
			fread(dd[0], sizeof(float), nr*nt, swap);
			sf_floatwrite(dd[0], nr*nt, Fdat);
		}
		fclose(swap);
		remove("temswap.bin");
	}
	MPI_Barrier(comm);
	
	/* release allocated memory */
	free(*vv); free(vv); free(*den); free(den);
	free(*v2d); free(v2d); free(*iden); free(iden);
	free(*tau); free(tau); free(*taus); free(taus);

	free(*p); free(p); free(*vx); free(vx); 
	free(*vz); free(vz); free(rr);
	free(*r1); free(r1); free(*r2); free(r2);
	free(*term); free(term); free(*term1); free(term1);
}
