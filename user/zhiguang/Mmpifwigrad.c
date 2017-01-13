/* Calculate FWI gradient */
/*
 Copyright (C) 2016 University of Texas at Austin
 
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

const float c0=-205./72, c1=8./5, c2=-1./5, c3=8./315, c4=-1./560;
void pad2d(float *vec, float **array, int nz, int nx, int nb);
void source_map(int sx, int sz, int rectx, int rectz, int padnx, int padnz, int padnzx, float *rr);
void laplace(float **p1, float **term, int padnx, int padnz, float dx2, float dz2);
void apply_sponge(float **p, float *bc, int padnx, int padnz, int nb);

int main(int argc, char* argv[])
{
	int ix, iz, is, ir, it;
	int nz, nx, nt, ns, nr, sz, rz, sx, rx, ds_v, s0_v, dr_v, r0_v, nb, padnx, padnz, padnzx, nzx, frectx, frectz;
	float dz, z0, dx, x0, dt, t0, ds, s0, dr, r0, temp, coef, dx2, dz2, dt2;
	float **vv, *ww, **p0, **p1, **p2, **term, **tmparray, *rr, ***wave, **pp, **mm, *bc;

	int cpuid, numprocs;
	float *sendbuf, *recvbuf;
	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Fw, Fadj, Fgrad;

	sf_file Fwfl1, Fwfl2;

	sf_init(argc, argv);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &cpuid);
	MPI_Comm_size(comm, &numprocs);

	Fv=sf_input("Fvel"); /* velocity model */
	Fw=sf_input("Fwavelet"); /* wavelet */
	Fadj=sf_input("Fadj"); /* adjoint source */
	Fgrad=sf_output("Fgrad"); /* gradient */

	Fwfl1=sf_output("Fwfl1");
	Fwfl2=sf_output("Fwfl2");

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
	if(!sf_getint("sz", &sz)) sz=5; /* source depth */
	if(!sf_getint("nr", &nr)) nr=nx; /* number of receiver */
	if(!sf_getfloat("dr", &dr)) dr=dx; /* receiver interval */
	if(!sf_getfloat("r0", &r0)) r0=x0; /* receiver origin */
	if(!sf_getint("rz", &rz)) rz=5; /* receiver depth */

	if(!sf_getint("frectx", &frectx)) frectx=2; /* source smoothing in x */
	if(!sf_getint("frectz", &frectz)) frectz=2; /* source smoothing in z */
	if(!sf_getint("nb", &nb)) nb=100; /* boundary width */
	if(!sf_getfloat("coef", &coef)) coef=0.002; /* absorbing boundary coefficient */

	/* absorbing boundary coefficients */
	bc=sf_floatalloc(nb);
	for(it=0; it<nb; it++){
		temp=coef*(nb-it);
		bc[it]=expf(-temp*temp);
	}

	/* acquisition parameters */
	padnx=nx+2*nb;
	padnz=nz+2*nb;
	padnzx=padnx*padnz;
	nzx=nz*nx;
	ds_v=ds/dx+0.5;
	s0_v=(s0-x0)/dx+0.5+nb;
	sz += nb;
	dr_v=dr/dx+0.5;
	r0_v=r0/dx+0.5+nb;
	rz += nb;
	dx2=dx*dx;
	dz2=dz*dz;
	dt2=dt*dt;

	/* memory allocation */
	vv=sf_floatalloc2(padnz,padnx);
	ww=sf_floatalloc(nt);
	pp=sf_floatalloc2(nt, nr);
	mm=sf_floatalloc2(nz, nx);
	rr=sf_floatalloc(padnzx);
	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	wave=sf_floatalloc3(nz, nx, nt);

	sf_floatread(mm[0], nzx, Fv);
	pad2d(mm[0], vv, nz, nx, nb);
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

	sf_putint(Fwfl1, "n1", padnz);
	sf_putint(Fwfl1, "n2",padnx);
	sf_putint(Fwfl1, "n3",(nt-1)/50+1);
	sf_putint(Fwfl2, "n1", padnz);
	sf_putint(Fwfl2, "n2",padnx);
	sf_putint(Fwfl2, "n3",(nt-1)/50+1);

	for(is=cpuid; is<ns; is+=numprocs){
		if(cpuid==0) sf_warning("###### is=%d ######", is+1);

		/* read adjoint source */
		sf_seek(Fadj, is*nr*nt*sizeof(float), SEEK_SET);
		sf_floatread(pp[0], nr*nt, Fadj);

		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		
		sx=s0_v+is*ds_v;
		source_map(sx, sz, frectx, frectz, padnx, padnz, padnzx, rr);

		/* forward propagation */
		for(it=0; it<nt; it++){

			/* save wavefield */
#ifdef _OPENMP 
#pragma omp parallel for \
			private(ix,iz) \
			shared(wave,p1,it,nb,nx,nz)
#endif
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[it][ix][iz]=p1[ix+nb][iz+nb];

			if(is==ns/2 && it%50==0) sf_floatwrite(p1[0],padnzx, Fwfl1);
			
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
				} // iz
			} // ix

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

		/* initialization */
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		memset(term[0], 0., padnzx*sizeof(float));
		
		/* backward propagation */
		for(it=nt-1; it>=0; it--){
			
			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* load data residual */
			for(ir=0; ir<nr; ir++){
				rx=r0_v+ir*dr_v;
				term[rx][rz] -= pp[ir][it];
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
			
			if(is==ns/2 && it%50==0) sf_floatwrite(p1[0],padnzx, Fwfl2);

			/* calculate gradient  */
			if(it>0 && it<nt-1){
#ifdef _OPENMP 
#pragma omp parallel for \
				private(ix,iz,temp) \
				shared(nx,nz,vv,wave,p1,it,dt2,mm)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						temp=vv[ix+nb][iz+nb];
						temp=temp*temp*temp;
						temp=-2./temp;
						mm[ix][iz] += (wave[it+1][ix][iz]-2.*wave[it][ix][iz]+wave[it-1][ix][iz])/dt2*p1[ix+nb][iz+nb]*temp;
					}
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
	
	sf_fileclose(Fwfl1);
	sf_fileclose(Fwfl2);

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

void pad2d(float *vec, float **array, int nz, int nx, int nb)
/*< convert a vector to an array >*/
{
	int ix, iz;
	
	for(ix=0; ix<nx; ix++){
		for(iz=0; iz<nz; iz++){
			array[ix+nb][iz+nb]=vec[ix*nz+iz];
		}
	}

    for (ix=nb; ix<nx+nb; ix++){
		for (iz=0; iz<nb; iz++){
			array[ix][iz]=array[ix][nb];
			array[ix][iz+nz+nb]=array[ix][nz+nb-1];
		}
	}

	for (ix=0; ix<nb; ix++){
		for (iz=0; iz<nz+2*nb; iz++){
			array[ix][iz]=array[nb][iz];
			array[ix+nx+nb][iz]=array[nx+nb-1][iz];
		}
	}
}

void source_map(int sx, int sz, int rectx, int rectz, int padnx, int padnz, int padnzx, float *rr)
/*< generate source map >*/
{
	int i, j, i0;
	int n[2], s[2], rect[2];
	bool diff[2], box[2];
	sf_triangle tr;

	n[0]=padnz; n[1]=padnx;
	s[0]=1; s[1]=padnz;
	rect[0]=rectz; rect[1]=rectx;
	diff[0]=false; diff[1]=false;
	box[0]=false; box[1]=false;

	for(i=0; i<padnzx; i++)
		rr[i]=0.;
	j=sx*padnz+sz;
	rr[j]=1.;

	for (i=0; i<2; i++){
		if(rect[i] <=1) continue;
		tr=sf_triangle_init(rect[i], n[i], box[i]);
		for(j=0; j<padnzx/n[i]; j++){
			i0=sf_first_index(i,j,2,n,s);
			sf_smooth2(tr,i0,s[i],diff[i],rr);
		}
		sf_triangle_close(tr);
	}
}

void laplace(float **p1, float **term, int padnx, int padnz, float dx2, float dz2)
/*< laplace operator >*/
{
	int ix, iz;

#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) \
	shared(padnx,padnz,p1,term)
#endif

	for (ix=4; ix<padnx-4; ix++){
		for (iz=4; iz<padnz-4; iz++){
			term[ix][iz] = 
				(c0*p1[ix][iz]
				+c1*(p1[ix+1][iz]+p1[ix-1][iz])
				+c2*(p1[ix+2][iz]+p1[ix-2][iz])
				+c3*(p1[ix+3][iz]+p1[ix-3][iz])
				+c4*(p1[ix+4][iz]+p1[ix-4][iz]))/dx2 
				+(c0*p1[ix][iz]
				+c1*(p1[ix][iz+1]+p1[ix][iz-1])
				+c2*(p1[ix][iz+2]+p1[ix][iz-2])
				+c3*(p1[ix][iz+3]+p1[ix][iz-3])
				+c4*(p1[ix][iz+4]+p1[ix][iz-4]))/dz2;
		}
	}
}

void apply_sponge(float **p, float *bc, int padnx, int padnz, int nb)
/*< apply absorbing boundary condition >*/
{
	int ix, iz;

	for (ix=0; ix<padnx; ix++){
		for(iz=0; iz<nb; iz++){	// top ABC
			p[ix][iz]=bc[iz]*p[ix][iz];
		}
		for(iz=padnz-nb; iz<padnz; iz++){ // bottom ABC			
			p[ix][iz]=bc[padnz-iz-1]*p[ix][iz];
		} 
	}

	for (iz=0; iz<padnz; iz++){
		for(ix=0; ix<nb; ix++){ // left ABC			
			p[ix][iz]=bc[ix]*p[ix][iz];
		}
		for(ix=padnx-nb; ix<padnx; ix++){ // right ABC			
			p[ix][iz]=bc[padnx-ix-1]*p[ix][iz];
		}
	}
}
