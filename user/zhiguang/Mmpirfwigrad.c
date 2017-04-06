/* Calculate acoustic Reflection FWI gradient with the prepared adjoint source (velocity-impedance scale separation) */
/*
 Copyright (C) 2017 University of Texas at Austin
 
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

int main(int argc, char* argv[])
{
	bool verb;
	int i, ix, iz, is, ir, it, wit;
	int nz, nx, padnx, padnz, padnzx, nzx, nb; 
	int ns, nr, sz, rz, sx, rx, ds_v, s0_v, dr_v, r0_v; 
	int nt, frectx, frectz, interval, wnt;
	float dz, z0, dx, x0, dt, t0, ds, s0, dr, r0, idt, wdt, coef, temp;

	float **vv, **den, **den0, *ww, *bcxp, *bczp, *bcxv, *bczv, **mm;
	float **v2d, **iden, **v2d0, **iden0, **pp, **term, *rr; 
	float **px, **pz, **p, **vx, **vz, ***wp, ***wvx, ***wvz;
	float **px0, **pz0, **p0, **vx0, **vz0, ***wp0, ***wvx0, ***wvz0;

	int cpuid, numprocs;
	float *sendbuf, *recvbuf;
	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Fd, Fd0, Fw, Fadj, Fgrad;

	sf_init(argc, argv);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &cpuid);
	MPI_Comm_size(comm, &numprocs);

	Fv=sf_input("Fvel"); /* velocity model */
	Fd=sf_input("Fd"); /* density/impedance model */
	Fd0=sf_input("Fd0"); /* smooth density/impedance model */
	Fw=sf_input("Fwavelet"); /* wavelet */
	Fadj=sf_input("Fadj"); /* adjoint source */
	Fgrad=sf_output("Fgrad"); /* gradient */

	if(!sf_histint(Fv, "n1", &nz)) sf_error("No n1= in Fv");
	if(!sf_histint(Fv, "n2", &nx)) sf_error("No n2= in Fv");
	if(!sf_histfloat(Fv, "d1", &dz)) sf_error("No d1= in Fv");
	if(!sf_histfloat(Fv, "d2", &dx)) sf_error("No d2= in Fv");
	if(!sf_histfloat(Fv, "o1", &z0)) sf_error("No o1= in Fv");
	if(!sf_histfloat(Fv, "o2", &x0)) sf_error("No o2= in Fv");
	if(!sf_histint(Fw, "n1", &nt)) sf_error("No n1= in Fw");
	if(!sf_histfloat(Fw, "d1", &dt)) sf_error("No d1= in Fw");
	if(!sf_histfloat(Fw, "o1", &t0)) sf_error("No o1= in Fw");

	if(!sf_getint("ns", &ns)) sf_error("shot number required"); /* shot number */
	if(!sf_getfloat("ds", &ds)) sf_error("shot interval required"); /* shot interval */
	if(!sf_getfloat("s0", &s0)) sf_error("shot origin required"); /* shot origin */
	if(!sf_getint("sz", &sz)) sz=3; /* source depth */
	if(!sf_getint("nr", &nr)) nr=nx; /* number of receiver */
	if(!sf_getfloat("dr", &dr)) dr=dx; /* receiver interval */
	if(!sf_getfloat("r0", &r0)) r0=x0; /* receiver origin */
	if(!sf_getint("rz", &rz)) rz=3; /* receiver depth */

	if(!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
	if(!sf_getint("frectx", &frectx)) frectx=2; /* source smoothing in x */
	if(!sf_getint("frectz", &frectz)) frectz=2; /* source smoothing in z */
	if(!sf_getint("nb", &nb)) nb=20; /* PML boundary width */
	if(!sf_getfloat("coef", &coef)) coef=5.; /* maximum velocity of the medium */
	if(!sf_getint("interval", &interval)) interval=1; /* wavefield storing interval */

	/* padding variables */
	padnx=nx+2*nb;
	padnz=nz+2*nb;
	padnzx=padnx*padnz;
	nzx=nz*nx;
	wnt=(nt-1)/interval+1;
	idt=1./dt;
	wdt=dt*interval*2.;

	/* absorbing boundary coefficients */
	bcxp=sf_floatalloc(padnx);
	bczp=sf_floatalloc(padnz);
	bcxv=sf_floatalloc(padnx);
	bczv=sf_floatalloc(padnz);
	memset(bcxp, 0., padnx*sizeof(float));
	memset(bczp, 0., padnz*sizeof(float));
	memset(bcxv, 0., padnx*sizeof(float));
	memset(bczv, 0., padnz*sizeof(float));

	for(i=0; i<nb; i++){
		temp=1.5*coef/nb*logf(10000.)/nb/nb;

		bcxp[nb-i-1]=temp/dx*(i+1)*(i+1)/2;
		bcxp[nb+nx+i]=temp/dx*(i+0.5)*(i+0.5)/2;

		bczp[nb-i-1]=temp/dz*(i+1)*(i+1)/2;
		bczp[nb+nz+i]=temp/dz*(i+0.5)*(i+0.5)/2;

		bcxv[nb-i-1]=temp/dx*(i+0.5)*(i+0.5)/2;
		bcxv[nb+nx+i]=temp/dx*(i+1)*(i+1)/2;

		bczv[nb-i-1]=temp/dz*(i+0.5)*(i+0.5)/2;
		bczv[nb+nz+i]=temp/dz*(i+1)*(i+1)/2;
	}

	/* acquisition parameters */
	ds_v=ds/dx+0.5;
	s0_v=(s0-x0)/dx+0.5+nb;
	sz += nb;
	dr_v=dr/dx+0.5;
	r0_v=r0/dx+0.5+nb;
	rz += nb;

	/* memory allocation */
	//float **vv, **den, **den0, *ww, **mm, *bcxp, *bczp, *bcxv, *bczv;
	//float **v2d, **iden, **v2d0, **iden0, **pp, **term, *rr; 
	//float **px, **pz, **p, **vx, **vz, ***wp, ***wvx, ***wvz;
	//float **px0, **pz0, **p0, **vx0, **vz0, ***wp0, ***wvx0, ***wvz0;
	vv=sf_floatalloc2(padnz,padnx);
	den=sf_floatalloc2(padnz,padnx);
	den0=sf_floatalloc2(padnz,padnx);
	ww=sf_floatalloc(nt);
	mm=sf_floatalloc2(nz, nx);

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

	/* read files */
	sf_floatread(mm[0], nzx, Fv);
	pad2d(mm[0], vv, nz, nx, nb);
	sf_floatread(mm[0], nzx, Fd);
	pad2d(mm[0], den, nz, nx, nb);
	sf_floatread(mm[0], nzx, Fd0);
	pad2d(mm[0], den0, nz, nx, nb);
	memset(mm[0], 0., nzx*sizeof(float));
	sf_floatread(ww, nt, Fw);

	/* dimension set up */
	sf_putint(Fgrad, "n1", nz);
	sf_putfloat(Fgrad, "d1", dz);
	sf_putfloat(Fgrad, "o1", z0);
	sf_putstring(Fgrad, "label1", "Depth");
	sf_putstring(Fgrad, "unit1", "km");
	sf_putint(Fgrad, "n2", nx);
	sf_putfloat(Fgrad, "d2", dx);
	sf_putfloat(Fgrad, "o2", x0);
	sf_putstring(Fgrad, "label2", "Distance");
	sf_putstring(Fgrad, "unit2", "km");

	/* prepare v2d and iden */
	for (ix=0; ix<padnx; ix++){
		for (iz=0; iz<padnz; iz++){
			v2d[ix][iz]=den[ix][iz]*vv[ix][iz];
			iden[ix][iz]=vv[ix][iz]/den[ix][iz];
			v2d0[ix][iz]=den0[ix][iz]*vv[ix][iz];
			iden0[ix][iz]=vv[ix][iz]/den0[ix][iz];
		}
	}

	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		/* read adjoint source */
		sf_seek(Fadj, is*nr*nt*sizeof(float), SEEK_SET);
		sf_floatread(pp[0], nr*nt, Fadj);

		/* source location */
		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		/*********************u and u0*********************************/
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

		wit=0;
		/* forward propagation */
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation is=%d; it=%d;", is+1, it);

			// calculate u *********************************
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
		
		/*********************u and u0*********************************/
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
			for(ir=0; ir<nr; ir++){
				rx=r0_v+ir*dr_v;
				term[rx][rz] += pp[ir][it];
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
		
			// calculate lambda0
			/* update px0 and load source */
			derivpx(vx0, term, padnx, padnz, dx);

			/* load adjoint source */
			for(ir=0; ir<nr; ir++){
				rx=r0_v+ir*dr_v;
				term[rx][rz] += pp[ir][it];
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
							//temp=-1./v2d[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							//mm[ix][iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p0[ix+nb][iz+nb]*temp;
							//temp=-1./iden[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							//mm[ix][iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx0[ix+nb][iz+nb]*temp;
							//mm[ix][iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz0[ix+nb][iz+nb]*temp;

							temp=-1./v2d0[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							mm[ix][iz] += (wp0 [wit+1][ix][iz]-wp0 [wit-1][ix][iz])/wdt*(p[ix+nb][iz+nb]-p0[ix+nb][iz+nb])*temp;
							mm[ix][iz] += (wp [wit+1][ix][iz]-wp [wit-1][ix][iz])/wdt*p0[ix+nb][iz+nb]*temp;

							temp=-1./iden0[ix+nb][iz+nb]/vv[ix+nb][iz+nb];
							mm[ix][iz] += (wvx0[wit+1][ix][iz]-wvx0[wit-1][ix][iz])/wdt*(vx[ix+nb][iz+nb]-vx0[ix+nb][iz+nb])*temp;
							mm[ix][iz] += (wvx[wit+1][ix][iz]-wvx[wit-1][ix][iz])/wdt*vx0[ix+nb][iz+nb]*temp;

							mm[ix][iz] += (wvz0[wit+1][ix][iz]-wvz0[wit-1][ix][iz])/wdt*(vz[ix+nb][iz+nb]-vz0[ix+nb][iz+nb])*temp;
							mm[ix][iz] += (wvz[wit+1][ix][iz]-wvz[wit-1][ix][iz])/wdt*vz0[ix+nb][iz+nb]*temp;
						} // end of iz
					} // end of ix
				} // avoid first and last time step
				wit--;
			} // end of 'if it%interval==0'

		} // end of time loop

	}// end of shot loop
	MPI_Barrier(comm);

	/* gradient reduction */
	if(cpuid==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=mm[0];
	}else{
		sendbuf=mm[0];
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nzx, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(mm[0], nzx, MPI_FLOAT, 0, comm);

	if(cpuid==0) sf_floatwrite(mm[0], nzx, Fgrad);

	MPI_Finalize();
	exit(0);
}
