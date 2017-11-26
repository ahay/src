/* Acoustic wave equation forward modeling with MPI and OpenMP */
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

static int nx, nz, nt, nb, padnx, padnz, padnzx;
static float dx, dz, dt, dx2, dz2, dt2;
float **rho, *bc, **px, **pz, **rhox, **rhoz;

const float c1=4./5, c2=-1./5, c3=4./105, c4=-1./280;
const float d0=-205./72, d1=8./5, d2=-1./5, d3=8./315, d4=-1./560;

void pad2d(float **in, float **out);
void apply_sponge(float **p);
void source_map(int sx, int sz, int rectx, int rectz, float *rr);
void cal_grad(float **irho);
void cal_div(float **p1, float **term);

int main(int argc, char* argv[])
{
	bool verb, wfl;
	int ix, iz, is, ir, it, ib;
	int ns, nr, sx, rx, sz, rz, rectx, rectz;
	int ds_v, dr_v, s0_v, r0_v;
	int cpuid, numprocs;

	float ds, dr, s0, r0, x0;
	float tmp, coef;
	float **vv, *ww, **dd, **irho;
	float **read, **p0, **p1, **p2, *rr, **term, **tmparray;

	FILE *swap;

	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Frho, Fw, Fdat;

	sf_init(argc,argv);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &cpuid);
	MPI_Comm_size(comm, &numprocs);

#ifdef _OPENMP
	omp_init();
#endif

	Fv = sf_input("--input"); /* veloctiy model */
	Frho = sf_input("rho"); /* density */
	Fw = sf_input("wavelet"); /* source wavelet */
	Fdat = sf_output("--output"); /* shot data */

	swap=fopen("temswap.bin","wb+");

	if (!sf_histint(Fv, "n1", &nz)) sf_error("No n1= in input"); 
	if (!sf_histint(Fv, "n2", &nx)) sf_error("No n2= in input");
	if (!sf_histfloat(Fv, "d1", &dz)) sf_error("No d1= in input");
	if (!sf_histfloat(Fv, "d2", &dx)) sf_error("No d2= in input");
	if (!sf_histfloat(Fv, "o2", &x0)) sf_error("No o2= in input");
	if (!sf_histint(Fw, "n1", &nt)) sf_error("No n1= in Fw"); /* number of time steps */
	if (!sf_histfloat(Fw, "d1", &dt)) sf_error("No d1= in Fw"); /* time sampling interval */

	/* command-line parameters */
	if (!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
	if (!sf_getint("ns", &ns)) sf_error("shot number required"); /* shot number */
	if (!sf_getfloat("ds", &ds)) sf_error("shot interval required"); /* shot interval */
	if (!sf_getfloat("s0", &s0)) sf_error("shot origin required"); /* shot origin */
	if (!sf_getint("sz", &sz)) sz=5; /* source depth */

	if (!sf_getint("nr", &nr)) sf_error("receiver number required"); /* number of receiver */
	if (!sf_getfloat("dr", &dr)) sf_error("receiver interval required"); /* receiver interval */
	if (!sf_getfloat("r0", &r0)) sf_error("receiver origin required"); /* receiver origin */
	if (!sf_getint("rz", &rz)) rz=sz; /* receiver depth */

	if (!sf_getint("nb", &nb)) nb=80; /* boundary width */
	if (!sf_getfloat("coef", &coef)) coef=0.003; /* absorbing boundary coefficient */

	if (!sf_getint("rectx", &rectx)) rectx=2; /* source smooothing parameter */
	if (!sf_getint("rectz", &rectz)) rectz=2; /* source smooothing parameter */

	/* correct the coordinate */
	padnx=nx+2*nb;
	padnz=nz+2*nb;
	padnzx=padnz*padnx;
        x0=x0-nb*dx;
        sz=sz+nb;
        rz=rz+nb;

        ds_v=ds/dx+0.5;
        s0_v=(s0-x0)/dx+0.5;
        dr_v=dr/dx+0.5;

        /* check parameters */
        if(cpuid==0){
            sf_warning("numprocs=%d", numprocs);
            sf_warning("nx=%d padnx=%d nz=%d padnz=%d nt=%d nr=%d ns=%d", nx, padnx, nz, padnz, nt, nr, ns);
            sf_warning("dx=%.3f dz=%.3f  dt=%.3f dr=%.3f ds=%.3f", dx, dz, dt, dr, ds);
            sf_warning("x0=%.3f r0=%.3f s0=%.3f sz=%d rz=%d", x0, r0, s0, sz, rz);
            sf_warning("dr_v=%d ds_v=%d s0_v=%d", dr_v, ds_v, s0_v);
        }

        /* set up the dimension of output data file */
        sf_putint(Fdat, "n1", nt);
        sf_putfloat(Fdat, "d1", dt);
        sf_putfloat(Fdat, "o1", 0.);
        sf_putstring(Fdat, "label1", "Time");
        sf_putstring(Fdat, "unit1", "s");
        sf_putint(Fdat, "n2", nr);
        sf_putfloat(Fdat, "d2", dr);
        sf_putfloat(Fdat, "o2", r0);
        sf_putstring(Fdat, "label2", "Offset");
        sf_putstring(Fdat, "unit2", "km");
        sf_putint(Fdat, "n3", ns);
        sf_putfloat(Fdat, "d3", ds);
        sf_putfloat(Fdat, "o3", s0);
        sf_putstring(Fdat, "label3", "Shot");
        sf_putstring(Fdat, "unit3", "km");

        /* allocate variables */
        vv  = sf_floatalloc2(padnz, padnx);
        rho = sf_floatalloc2(padnz, padnx);
        ww  = sf_floatalloc(nt);
        bc  = sf_floatalloc(nb);
        rr  = sf_floatalloc(padnzx);
	dd=sf_floatalloc2(nt, nr);

	read= sf_floatalloc2(nz, nx);
	p0  = sf_floatalloc2(padnz, padnx);
	p1  = sf_floatalloc2(padnz, padnx);
	p2  = sf_floatalloc2(padnz, padnx);
	term= sf_floatalloc2(padnz, padnx);
	irho= sf_floatalloc2(padnz, padnx);

	px  = sf_floatalloc2(padnz, padnx);
	pz  = sf_floatalloc2(padnz, padnx);
	rhox  = sf_floatalloc2(padnz, padnx);
	rhoz  = sf_floatalloc2(padnz, padnx);

	/* read parameters */
	sf_floatread(read[0], nz*nx, Fv);
	pad2d(read, vv);
	sf_floatread(read[0], nz*nx, Frho);
        pad2d(read, rho);
        sf_floatread(ww, nt, Fw);

        /* boundary absorbing coefficients */
        for(ib=0; ib<nb; ib++){
            tmp=coef*(nb-ib);
            bc[ib]=expf(-tmp*tmp);
        }

        /* gradient of density */
        for (ix=0; ix<padnx; ix++){
            for (iz=0; iz<padnz; iz++){
                irho[ix][iz]=1./rho[ix][iz];
            }
        }
        cal_grad(irho);

        dt2=dt*dt;
        dx2=dx*dx;
        dz2=dz*dz;

        /* shot loop */
        for (is=cpuid; is<ns; is+=numprocs){

            sf_warning("### is=%d ###", is+1);

            sx=s0_v+is*ds_v;
            r0_v=(s0+r0-x0)/dx+0.5;

            /* get source weighting map for stability */
            source_map(sx, sz, rectx, rectz, rr);

            /* initialize wavefield variables */
            memset(dd[0], 0, nr*nt*sizeof(float));
            memset(p0[0], 0, padnz*padnx*sizeof(float));
            memset(p1[0], 0, padnz*padnx*sizeof(float));
            memset(p2[0], 0, padnz*padnx*sizeof(float));

            /* time loop */
            for(it=0; it<nt; it++){

                if(verb) sf_warning(" Modeling is=%d; it=%d; ", is+1, it);

                /* output data */
#ifdef _OPENMP
#pragma omp parallel for \
                private(ir,rx) \
                shared(nr,r0_v,dr_v,dd,p1)
#endif
                for(ir=0; ir<nr; ir++){
                    rx=r0_v+ir*dr_v;
                    dd[ir][it]=p1[rx][rz];
                }

                /* calculate divergence */
                cal_div(p1, term);

                /* load source */
#ifdef _OPENMP
#pragma omp parallel for \
                private(ix,iz) \
                shared(padnx, padnz, rr, ww, term, it)
#endif
                for(ix=0; ix<padnx; ix++){
                    for(iz=0; iz<padnz; iz++){ 
                        term[ix][iz] += rr[ix*padnz+iz]*ww[it];
                    }
                }

                /* update pressure */
#ifdef _OPENMP
#pragma omp parallel for \
                private(ix,iz) \
                shared(padnx, padnz, p2, p1, p0, vv, dt2, term)
#endif
                for (ix=0; ix<padnx; ix++){
                    for (iz=0; iz<padnz; iz++){
                        p2[ix][iz]=
                            2*p1[ix][iz] - p0[ix][iz]
                            + vv[ix][iz]*dt2*term[ix][iz];
                    }
                }

                /* swap wavefield pointer of different time steps */
                tmparray=p0; p0=p1; p1=p2; p2=tmparray;

                apply_sponge(p0);
                apply_sponge(p1);
            } // end of time loop

            fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
            fwrite(dd[0], sizeof(float), nr*nt, swap);

        }// end of shot loop
        fclose(swap);
        MPI_Barrier(comm);

        /* transfer data to Fdat */
        if(cpuid==0){
            swap=fopen("temswap.bin", "rb");
            for(is=0; is<ns; is++){
                fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
                fread(dd[0], sizeof(float), nr*nt, swap);
                sf_floatwrite(dd[0], nr*nt, Fdat);
            }
            fclose(swap);
            remove("temswap.bin");
        }
        MPI_Barrier(comm);

        free(*vv); free(vv);
        free(*rho); free(rho);
        free(bc);
        free(*px); free(px);
        free(*pz); free(pz);
        free(*rhox); free(rhox);
        free(*rhoz); free(rhoz);

        MPI_Finalize();
        exit(0);
}

void pad2d(float **in, float **out)
    /*< expand domain >*/
{
    int iz, ix;

    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            out[nb+ix][nb+iz]=in[ix][iz];
        }
    }

	for (ix=nb; ix<nx+nb; ix++){
		for (iz=0; iz<nb; iz++){
			out[ix][iz]=out[ix][nb];
			out[ix][iz+nz+nb]=out[ix][nz+nb-1];
		}
	}

	for (ix=0; ix<nb; ix++){
		for (iz=0; iz<padnz; iz++){
			out[ix][iz]=out[nb][iz];
			out[ix+nx+nb][iz]=out[nx+nb-1][iz];
		}
	}
}

void apply_sponge(float **p)
	/*< apply absorbing boundary condition >*/
{
	int iz, ix;

#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) \
	shared(padnx, padnz, nb, p, bc)
#endif
	for (ix=0; ix<padnx; ix++){
		for(iz=0; iz<nb; iz++){	// top ABC
			p[ix][iz]=bc[iz]*p[ix][iz];
		}
		for(iz=nz+nb; iz<padnz; iz++){ // bottom ABC			
			p[ix][iz]=bc[padnz-iz-1]*p[ix][iz];
		} 
	}

#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) \
	shared(padnx, padnz, nb, p, bc)
#endif
	for (iz=0; iz<padnz; iz++){
		for(ix=0; ix<nb; ix++){ // left ABC			
			p[ix][iz]=bc[ix]*p[ix][iz];
		}	
		for(ix=nx+nb; ix<padnx; ix++){ // right ABC			
			p[ix][iz]=bc[padnx-ix-1]*p[ix][iz];
		}	
	}
}

void source_map(int sx, int sz, int rectx, int rectz, float *rr)
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

#ifdef _OPENMP
#pragma omp parallel for \
	private(i) \
	shared(padnzx, rr)
#endif
	for (i=0; i<padnzx; i++)
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

void cal_div(float **p1, float **term)
	/*< calculate divergence term >*/
{
	int ix, iz;

#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) \
	shared(padnz, padnx, px, pz, p1, dx, dz, rho, rhox, rhoz)
#endif
	for (ix=4; ix<padnx-4; ix++){
		for (iz=4; iz<padnz-4; iz++){
			px[ix][iz] = (c1*(p1[ix+1][iz]-p1[ix-1][iz])
					+c2*(p1[ix+2][iz]-p1[ix-2][iz])
					+c3*(p1[ix+3][iz]-p1[ix-3][iz])
					+c4*(p1[ix+4][iz]-p1[ix-4][iz]))
				/dx*rho[ix][iz]*rhox[ix][iz];

			pz[ix][iz] = (c1*(p1[ix][iz+1]-p1[ix][iz-1])
					+c2*(p1[ix][iz+2]-p1[ix][iz-2])
					+c3*(p1[ix][iz+3]-p1[ix][iz-3])
					+c4*(p1[ix][iz+4]-p1[ix][iz-4]))
				/dz*rho[ix][iz]*rhoz[ix][iz];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for \
	private(ix,iz) \
	shared(padnz, padnx, term, p1, dx2, dz2, px, pz)
#endif
	for (ix=4; ix<padnx-4; ix++){
		for (iz=4; iz<padnz-4; iz++){
			term[ix][iz] = 
				(d0*p1[ix][iz]
				 +d1*(p1[ix+1][iz]+p1[ix-1][iz])
				 +d2*(p1[ix+2][iz]+p1[ix-2][iz])
				 +d3*(p1[ix+3][iz]+p1[ix-3][iz])
				 +d4*(p1[ix+4][iz]+p1[ix-4][iz]))/dx2 
				+(d0*p1[ix][iz]
						+d1*(p1[ix][iz+1]+p1[ix][iz-1])
						+d2*(p1[ix][iz+2]+p1[ix][iz-2])
						+d3*(p1[ix][iz+3]+p1[ix][iz-3])
						+d4*(p1[ix][iz+4]+p1[ix][iz-4]))/dz2
				+px[ix][iz]+pz[ix][iz];
		}
	}
}

void cal_grad(float **irho)
	/*< calculate the gradient of the inverse of density >*/
{
	int ix, iz;

	for (ix=4; ix<padnx-4; ix++){
		for (iz=4; iz<padnz-4; iz++){
			rhox[ix][iz] = (c1*(irho[ix+1][iz]-irho[ix-1][iz])
					+c2*(irho[ix+2][iz]-irho[ix-2][iz])
					+c3*(irho[ix+3][iz]-irho[ix-3][iz])
					+c4*(irho[ix+4][iz]-irho[ix-4][iz]))/dx;

			rhoz[ix][iz] = (c1*(irho[ix][iz+1]-irho[ix][iz-1])
					+c2*(irho[ix][iz+2]-irho[ix][iz-2])
					+c3*(irho[ix][iz+3]-irho[ix][iz-3])
					+c4*(irho[ix][iz+4]-irho[ix][iz-4]))/dz;
		}
	}
}
