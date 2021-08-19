/* 2D Visco-acoustic Forward Modeling, FWI, and RTM based on SLS model */
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

int omp_init()
/*< init OMP parameters >*/
{
    int ompnth;
    int ompchunk;
    
#ifdef _OPENMP
    int ompath;
#endif

    /* OMP data chunk size */
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    
#ifdef _OPENMP
    /* OMP available threads */
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#else
    ompnth=0;
#endif
    
    return ompnth;
}

#endif

#include <time.h>
#include <stdlib.h>

typedef struct sf_mpipar{
	int cpuid;
	int numprocs;
} sf_mpi;
/*^*/

typedef struct sf_source{
	float fhi;
	float flo;
	int rectx;
	int rectz;
} *sf_sou;
/*^*/

typedef struct sf_acquisition{
	// model dimension
	int nx;
	int nz;
	float dx;
	float dz;
	float x0;
	float z0;
	// wavelet dimension
	int nt;
	float dt;
	float t0;
	// absorbing boundary condition
	int nb;
	float coef;
	float *bc;
	// padding
	int padnx;
	int padnz;
	float padx0;
	float padz0;
	// acquisition type
	int acqui_type;
	// shot
	int ns;
	float ds;
	float s0;
	int sz;
	int ds_v;
	int s0_v;
	// receiver
	int nr;
	float dr;
	float r0;
	int rz;
	int dr_v;
	int *r0_v;
	int *r02;
	int *nr2;
	// reference frequency
	float f0;
	// wavefield storing interval
	int interval;
} *sf_acqui;
/*^*/

typedef struct sf_1darray{
	float *vv;
	float *qq;
	float *tau;
	float *taus;
	float *ww;
} *sf_vec;
/*^*/

typedef struct sf_fwipar{
	bool onlygrad;
	int grad_type;
	int misfit_type;
	int opt_type;
	// data residual weighting
	float wt1;
	float wt2;
	float woff1;
	float woff2;
        bool oreo;
	// water layer depth
	int waterz;
	int waterzb;
	// gradient smoothing parameters
	int rectx;
	int rectz;
} *sf_fwi;
/*^*/

typedef struct sf_optimization {
	int niter;
	float conv_error;
	int npair;
	int nls;
	int igrad;
	int ipair;
	int ils;
        int repeat;
	float c1;
	float c2;
	float factor;
	float alpha;
	float f0;
	float fk;
	float gk_norm;
	float **sk, **yk;
        /* bound constraints */
        float v1;
        float v2;
} *sf_optim;
/*^*/

typedef struct sf_passive{
    bool inv;
    bool onlysrc;
    bool onlyvel;
    bool sw;
    bool ctr;
    bool prec;
    bool hidesrc;
    int niter;
    int ngrp;
    int size;
    int rectz;
    int rectx;
    int rectt;
    int repeat;
    float perc;
    float hard;
} *sf_pas;
/*^*/


typedef void (*sf_gradient)(float*,float*,float*);
/*^*/

const float c0=-205./72, c1=8./5, c2=-1./5, c3=8./315, c4=-1./560;

void preparation(sf_file Fv, sf_file Fq, sf_file Fw, sf_acqui acpar, sf_sou soupar, sf_vec array)
/*< read data, initialize variables and prepare acquisition geometry >*/
{
	int i, nb, nzx;
	float sx, xend, rbegin, rend;
	float *taue, tmp;

	int nplo=3, nphi=3, nt;
	float eps=0.0001;
	sf_butter blo=NULL, bhi=NULL;

	/* absorbing boundary coefficients */
	nb=acpar->nb;
	acpar->bc=sf_floatalloc(nb);
	for(i=0; i<nb; i++){
		tmp=acpar->coef*(nb-i);
		acpar->bc[i]=expf(-tmp*tmp);
	}

	/* padding variables */
	acpar->padnx=acpar->nx+2*nb;
	acpar->padnz=acpar->nz+2*nb;
	acpar->padx0=acpar->x0-nb*acpar->dx;
	acpar->padz0=acpar->z0-nb*acpar->dz;

	/* acquisition parameters */
	acpar->ds_v=acpar->ds/acpar->dx+0.5;
	acpar->s0_v=acpar->s0/acpar->dx+0.5+nb;
	acpar->sz += nb;

	acpar->dr_v=acpar->dr/acpar->dx+0.5;
	acpar->r0_v=sf_intalloc(acpar->ns);
	acpar->r02=sf_intalloc(acpar->ns);
	acpar->nr2=sf_intalloc(acpar->ns);
	acpar->rz += nb;
	xend=acpar->x0+(acpar->nx-1)*acpar->dx;
	if(acpar->acqui_type==1){
		for(i=0; i<acpar->ns; i++){
			acpar->r0_v[i]=(acpar->r0-acpar->x0)/acpar->dx+0.5+nb;
			acpar->r02[i]=0;
			acpar->nr2[i]=acpar->nr;
		}
	}else{
		for(i=0; i<acpar->ns; i++){
			sx=acpar->s0+acpar->ds*i;
			rbegin=(sx+acpar->r0 <acpar->x0)? acpar->x0 : sx+acpar->r0;
			rend=sx+acpar->r0 +(acpar->nr-1)*acpar->dr;
			rend=(rend < xend)? rend : xend;
			acpar->r0_v[i]=rbegin/acpar->dx+0.5+nb;
			acpar->r02[i]=(rbegin-sx-acpar->r0)/acpar->dx+0.5;
			acpar->nr2[i]=(rend-rbegin)/acpar->dx+1.5;
		}
	}

	/* read model parameters */
	nzx=acpar->nz*acpar->nx;
	nt=acpar->nt;
	array->vv=sf_floatalloc(nzx);
	array->qq=sf_floatalloc(nzx);
	array->tau=sf_floatalloc(nzx);
	array->taus=sf_floatalloc(nzx);
	array->ww=sf_floatalloc(nt);
	taue=sf_floatalloc(nzx);

	sf_floatread(array->vv, nzx, Fv);
	sf_floatread(array->qq, nzx, Fq);
	sf_floatread(array->ww, nt, Fw);

	/* calculate tau */
	for(i=0; i<nzx; i++){
		tmp=sqrtf(array->qq[i]*array->qq[i]+1);
		taue[i]=(tmp+1)/(2.*SF_PI*acpar->f0*array->qq[i]);
		array->taus[i]=(tmp-1)/(2.*SF_PI*acpar->f0*array->qq[i]);
		array->tau[i]=taue[i]/array->taus[i]-1;
	}

	/* bandpass the wavelet */
	soupar->flo *= acpar->dt;
	soupar->fhi *= acpar->dt;
	if(soupar->flo > eps) blo=sf_butter_init(false, soupar->flo, nplo);
	if(soupar->fhi < 0.5-eps) bhi=sf_butter_init(true, soupar->fhi, nphi);

	if(NULL != blo){
		sf_butter_apply(blo, nt, array->ww);
		sf_reverse(nt, array->ww);
		sf_butter_apply(blo, nt, array->ww);
		sf_reverse(nt, array->ww);
		sf_butter_close(blo);
	}
	if(NULL != bhi){
		sf_butter_apply(bhi, nt, array->ww);
		sf_reverse(nt, array->ww);
		sf_butter_apply(bhi, nt, array->ww);
		sf_reverse(nt, array->ww);
		sf_butter_close(bhi);
	}
	
	free(taue);
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

void laplace(float **p1, float **term, int padnx, int padnz, float dx2, float dz2)
/*< laplace operator >*/
{
	int ix, iz;

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

void residual_weighting(float **ww, int nt, int nr, int wtn1, int wtn2, int woffn1, int woffn2, bool oreo)
/*< data residual weighting >*/
{
	int it, ir;
	float w[10];

	for(it=0; it<10; it++){
		w[it]=sinf(0.5*SF_PI*(it+1)/11.);
	}

        if (oreo) {

            for(ir=0; ir<nr; ir++){
                for(it=0; it<nt; it++){
                    ww[ir][it]=1.;
                }
            }

            for(ir=woffn1; ir<=woffn2; ir++){
                for(it=wtn1; it<=wtn2; it++){
                    ww[ir][it]=0.;
                }
            }

            for(ir=0; ir<10; ir++){
                for(it=wtn1; it<=wtn2; it++){
                    if(woffn1-ir>0 ) ww[woffn1-ir][it] *= w[ir];
                    if(woffn2+ir<nr) ww[woffn2+ir][it] *= w[ir];
                }
            }

            for(it=0; it<10; it++){
                for(ir=0; ir<woffn1; ir++){
                    ww[ir][wtn1+it] *= w[it];
                    ww[ir][wtn2-it] *= w[it];
                }
                for(ir=woffn2+1; ir<nr; ir++){
                    ww[ir][wtn1+it] *= w[it];
                    ww[ir][wtn2-it] *= w[it];
                }
            }

        } else {

            for(ir=0; ir<nr; ir++){
                for(it=0; it<nt; it++){
                    ww[ir][it]=0.;
                }
            }

            for(ir=woffn1; ir<=woffn2; ir++){
                for(it=wtn1; it<=wtn2; it++){
                    ww[ir][it]=1.;
                }
            }

            for(ir=0; ir<10; ir++){
                for(it=wtn1; it<=wtn2; it++){
                    ww[woffn1+ir][it] *= w[ir];
                    ww[woffn2-ir][it] *= w[ir];
                }
            }

            for(it=0; it<10; it++){
                for(ir=woffn1; ir<=woffn2; ir++){
                    ww[ir][wtn1+it] *= w[it];
                    ww[ir][wtn2-it] *= w[it];
                }
            }

        }
}

void gradient_smooth2(int rectx, int rectz, int nx, int nz, int waterz, float scaling, float *grad)
/*< smooth gradient, zero bathymetry layer and normalization >*/
{
	int i, j, i0, nzx;
	int n[2], s[2], rect[2];
	bool diff[2], box[2];
	sf_triangle tr;

	nzx=nz*nx;
	n[0]=nz; n[1]=nx;
	s[0]=1; s[1]=nz;
	rect[0]=rectz; rect[1]=rectx;
	diff[0]=false; diff[1]=false;
	box[0]=false; box[1]=false;

	for (i=0; i<2; i++){
		if(rect[i] <=1) continue;
		tr=sf_triangle_init(rect[i], n[i], box[i]);
		for(j=0; j<nzx/n[i]; j++){
			i0=sf_first_index(i,j,2,n,s);
			sf_smooth2(tr,i0,s[i],diff[i],grad);
		}
		sf_triangle_close(tr);
	}

        if (waterz>=0) {
            for(i=0; i<waterz; i++)
                for(j=0; j<nx; j++)
                    grad[i+j*nz]=0.;
        } else {
            for(i=nz+waterz; i<nz; i++)
                for(j=0; j<nx; j++)
                    grad[i+j*nz]=0.;
        }

	for(i=0; i<nzx; i++)
		grad[i] *= scaling;
}

void gradient_smooth2b(int rectx, int rectz, int nx, int nz, int waterz, int waterzb, float scaling, float *grad)
/*< smooth gradient, zero bathymetry layer and normalization >*/
{
	int i, j, i0, nzx;
	int n[2], s[2], rect[2];
	bool diff[2], box[2];
	sf_triangle tr;

	nzx=nz*nx;
	n[0]=nz; n[1]=nx;
	s[0]=1; s[1]=nz;
	rect[0]=rectz; rect[1]=rectx;
	diff[0]=false; diff[1]=false;
	box[0]=false; box[1]=false;

	for (i=0; i<2; i++){
		if(rect[i] <=1) continue;
		tr=sf_triangle_init(rect[i], n[i], box[i]);
		for(j=0; j<nzx/n[i]; j++){
			i0=sf_first_index(i,j,2,n,s);
			sf_smooth2(tr,i0,s[i],diff[i],grad);
		}
		sf_triangle_close(tr);
	}

        for(i=0; i<waterz; i++)
            for(j=0; j<nx; j++)
                grad[i+j*nz]=0.;

        for(i=nz-waterzb; i<nz; i++)
            for(j=0; j<nx; j++)
                grad[i+j*nz]=0.;

	for(i=0; i<nzx; i++)
		grad[i] *= scaling;
}

void l2norm(int n, float *a, float *norm)
/*< L2 norm of a vector >*/
{
	int i;
	*norm=0.;
	for(i=0; i<n; i++){
		*norm += a[i]*a[i];
	}
	*norm=sqrtf(*norm);
}

void reverse(int n, float *a, float *b)
/*< reverse the sign of the vector >*/
{
	int i;
	for (i=0; i<n; i++)
		b[i]=-a[i];
}

void copy(int n, float *a, float *b)
/*< copy vector >*/
{
	int i;
	for(i=0; i<n; i++)
		b[i]=a[i];
}

void dot_product(int n, float *a, float *b, float *product)
/*< dot product of two vectors >*/
{
	int i;
	*product=0.;
	for (i=0; i<n; i++)
		*product += a[i]*b[i];
}

void print_iteration(FILE *fp, int iter, sf_optim opt)
/*< print out iteration information >*/
{
	if(iter%10==0){
		fprintf(fp,"*********************************************\n");
		fprintf(fp,"Maximum Iteration: %d\n", opt->niter);
		fprintf(fp,"Convergence Error: %3.2e\n", opt->conv_error);
		fprintf(fp,"*********************************************\n");
		fprintf(fp,"Niter  Misfit  Rel_Misfit  Grad_Norm  Alpha   Num_Pair  Num_LS  Total_Grad\n");
	}
	fprintf(fp,"%3d   %3.2e  %3.2e   %3.2e  %3.2e  %3d       %3d      %4d\n", iter, opt->fk, opt->fk/opt->f0, opt->gk_norm, opt->alpha, opt->ipair, opt->ils, opt->igrad);
}

sf_pas passive_init(sf_acqui acpar)
/*< read data, initialize variables and prepare acquisition geometry >*/
{
    sf_pas paspar;

    paspar = (sf_pas) sf_alloc(1,sizeof(*paspar));

    if(!sf_getbool("inv", &paspar->inv)) paspar->inv=false; /* inversion flag */
    if(!sf_getbool("onlysrc", &paspar->onlysrc)) paspar->onlysrc=false;  /* only invert for source (vel known), active when inv=y */
    if(!sf_getbool("onlyvel", &paspar->onlyvel)) paspar->onlyvel=false;  /* only invert for vel (source known), active when inv=y */
    if(!sf_getbool("sw", &paspar->sw)) paspar->sw=false;  /* sliding window normalization */
    if(!sf_getbool("ctr", &paspar->ctr)) paspar->ctr=false; /* cross-correlation time-reversal imaging */
    if(!sf_getbool("precsrc", &paspar->prec)) paspar->prec=false; /* source inversion preconditioning */
    if(!sf_getbool("hidesrc", &paspar->hidesrc)) paspar->hidesrc=false; /* hide source footprint in fwi */
    if(!sf_getint("nitersrc", &paspar->niter)) paspar->niter=1;   /* num of iter'ns for source inversion */
    if(!sf_getint("ngrp", &paspar->ngrp))   paspar->ngrp=1; /* number of sub-groups of receivers */
    if(!sf_getint("size", &paspar->size))   paspar->size=0; /* sliding window radius */
    if(!sf_getint("rectzsrc", &paspar->rectz)) paspar->rectz=1; /* source smoothing in z before masking */
    if(!sf_getint("rectxsrc", &paspar->rectx)) paspar->rectx=1; /* source smoothing in x before masking */
    if(!sf_getint("recttsrc", &paspar->rectt)) paspar->rectt=50; /* source smoothing in t before masking */
    if(!sf_getint("repeatsrc", &paspar->repeat)) paspar->repeat=1; /* source smoothing repeatation times */
    if(!sf_getfloat("perc", &paspar->perc)) paspar->perc=1.0f; /* padding percentatge for swnorm */
    if(!sf_getfloat("hard", &paspar->hard)) paspar->hard=0.1f; /* hard thresholding for masking */

    if(paspar->onlyvel && paspar->onlysrc) sf_error("Error: onlyvel and onlysrc cannot both be true!");
    if(paspar->inv) {
        if(paspar->onlyvel) sf_warning("inverting for velocity only ...");
        else if(paspar->onlysrc) sf_warning("inverting for source only ...");
        else sf_warning("inverting for both velocity and source ...");
    } else sf_warning("forward modeling ...");

    return paspar;
}


void forward_modeling_a(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, bool verb)
/*< acoustic forward modeling >*/
{
	int ix, iz, is, ir, it;
	int sx, rx, sz, rz, rectx, rectz;
	int nz, nx, padnz, padnx, padnzx, nt, nr, nb;

	float dx2, dz2, dt2;
	float **vv, **dd;
	float **p0, **p1, **p2, **term, **tmparray, *rr;

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
	rectx=soupar->rectx;
	rectz=soupar->rectz;

	dx2=acpar->dx*acpar->dx;
	dz2=acpar->dz*acpar->dz;
	dt2=acpar->dt*acpar->dt;

	vv = sf_floatalloc2(padnz, padnx);
	dd=sf_floatalloc2(nt, nr);

	p0=sf_floatalloc2(padnz, padnx);
	p1=sf_floatalloc2(padnz, padnx);
	p2=sf_floatalloc2(padnz, padnx);
	term=sf_floatalloc2(padnz, padnx);
	rr=sf_floatalloc(padnzx);

	/* padding and convert vector to 2-d array */
	pad2d(array->vv, vv, nz, nx, nb);

	for(is=mpipar->cpuid; is<acpar->ns; is+=mpipar->numprocs){
		sf_warning("###### is=%d ######", is+1);

		memset(dd[0], 0., nr*nt*sizeof(float));
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		
		sx=acpar->s0_v+is*acpar->ds_v;
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

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
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					term[ix][iz] += rr[ix*padnz+iz]*array->ww[it];
				}
			}

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
			if (!fread(dd[0], sizeof(float), nr*nt, swap))
				abort();
			sf_floatwrite(dd[0], nr * nt, Fdat);
		}
		fclose(swap);
		remove("temswap.bin");
	}
	MPI_Barrier(comm);
	
	/* release allocated memory */
	free(*p0); free(p0); free(*p1); free(p1);
	free(*p2); free(p2); free(*vv); free(vv);
	free(*dd); free(dd);
	free(rr); free(*term); free(term);
}

void forward_modeling(sf_file Fdat, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, bool verb)
/*< visco-acoustic forward modeling >*/
{
	int ix, iz, is, ir, it;
	int sx, rx, sz, rz, rectx, rectz;
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
	rectx=soupar->rectx;
	rectz=soupar->rectz;

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
		source_map(sx, sz, rectx, rectz, padnx, padnz, padnzx, rr);

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
			if (!fread(dd[0], sizeof(float), nr*nt, swap))
				abort();
			sf_floatwrite(dd[0], nr * nt, Fdat);
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


void lbfgs_save(int n, float *x, float *grad, float **sk, float **yk, sf_optim opt)
/*< save current model and gradient >*/
{
	int i;
	if(opt->ipair < opt->npair){
		copy(n, x, sk[opt->ipair]);
		copy(n, grad, yk[opt->ipair]);
		opt->ipair += 1;
	}else{
		for(i=0; i<opt->npair-1; i++){
			copy(n, sk[i+1], sk[i]);
			copy(n, yk[i+1], yk[i]);
		}
		copy(n, x, sk[opt->npair-1]);
		copy(n, grad, yk[opt->npair-1]);
	}
}

void lbfgs_update(int n, float *x, float *grad, float **sk, float **yk, sf_optim opt)
/*< update current sk and yk >*/
{
	int i, j;
	j=opt->ipair-1;
	for(i=0; i<n; i++){
		sk[j][i]=x[i]-sk[j][i];
		yk[j][i]=grad[i]-yk[j][i];
	}
}

void lbfgs_direction(int n, float *grad, float *r, float **sk, float **yk, sf_optim opt)
/*< calculate search direction >*/
{
	int i, j;
	float *rho, *q, *alpha, tmp, tmp1, gamma, beta;

	// safeguard
	l2norm(n, sk[opt->ipair-1], &tmp);
	l2norm(n, yk[opt->ipair-1], &tmp1);
	if(tmp==0. || tmp1==0.){
		reverse(n, grad, r);
		return;
	}
	
	q=sf_floatalloc(n);
	rho=sf_floatalloc(opt->ipair);
	alpha=sf_floatalloc(opt->ipair);

	copy(n, grad, q);
	
	// first loop
	for(i=opt->ipair-1; i>=0; i--){
		
		// calculate rho
		dot_product(n, yk[i], sk[i], &tmp);  
		rho[i]=1./tmp;

		dot_product(n, sk[i], q, &tmp);
		alpha[i]=rho[i]*tmp;
		for(j=0; j<n; j++)
			q[j] -= alpha[i]*yk[i][j];
	}

	// initial Hessian
	dot_product(n, yk[opt->ipair-1], yk[opt->ipair-1], &tmp);
	gamma=1./tmp/rho[opt->ipair-1];
	for (j=0; j<n; j++){
		r[j]=gamma*q[j];
	}

	// second loop
	for(i=0; i<opt->ipair; i++){
		dot_product(n, yk[i], r, &tmp);
		beta=tmp*rho[i];
		tmp=alpha[i]-beta;
		for(j=0; j<n; j++)
			r[j] += tmp*sk[i][j];
	}

	// opposite direction of H^(-1)*grad(f)
	for(j=0; j<n; j++)
		r[j]=-r[j];

	// deallocate variables
	free(q);
	free(alpha);
	free(rho);
}

void clip(float *x, int n, float min, float max)
/*< clip data >*/
{
    int i;

    for(i=0; i<n; i++){
        if(x[i]<min) x[i]=min;
        if(x[i]>max) x[i]=max;
    }
}

void line_search(int n, float *x, float *grad, float *direction, sf_gradient gradient, sf_optim opt, int *flag, int cpuid)
/*< line search (Wolfe condition) >*/
{
	int i, j;
	float m1, m2, m3, fcost, alpha1=0., alpha2=0.;
	float *xk;

	xk=sf_floatalloc(n);
	copy(n, x, xk);
	dot_product(n, grad, direction, &m1);
	m2=m1*opt->c2;
	m1=m1*opt->c1;
	
	for(i=0; i<opt->nls; i++){
		
		opt->ils += 1;
		for(j=0; j<n; j++)
			x[j] =xk[j] + opt->alpha*direction[j];

                clip(x, n, opt->v1, opt->v2);

		gradient(x, &fcost, grad);
		opt->igrad += 1;
		dot_product(n, grad, direction, &m3);
		
		if(cpuid==0){
			sf_warning("line search i=%d",i+1);
			sf_warning("alpha1=%g alpha2=%g alpha=%g",alpha1, alpha2, opt->alpha);
			sf_warning("fcost=%g fk=%g fk+c1*alpha*m1=%g m3=%g c2*m1=%g ",fcost, opt->fk, opt->fk+opt->alpha*m1, m3, m2);
		}

		if(fcost <= opt->fk + opt->alpha*m1 && m3 >= m2){
			opt->fk=fcost;
			*flag=0;
			break;
		}else if (fcost > opt->fk + opt->alpha*m1){
			alpha2=opt->alpha;
			opt->alpha=0.5*(alpha1+alpha2);
		}else{
			alpha1=opt->alpha;
			if(alpha2 == 0.)
				opt->alpha *= opt->factor;
			else
				opt->alpha = 0.5*(alpha1+alpha2);
		}

	}
	
	if(i==opt->nls){
		if(fcost <= opt->fk)
			*flag=1;
		else
			*flag=2;
	}

	free(xk);
}

#ifndef _triutil_h

#define NOP 4 /* derivative operator half-size */
#define C0 -205.0f/72.0f
#define C1 +8.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 +8.0f/315.0f
#define C4 -1.0f/560.0f
#define Lap(a,ix,iz,sx,sz,v)  ( ( C4*(a[ix+4][iz  ] + a[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ] + a[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ] + a[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ] + a[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]) )*sx            +      \
                                ( C4*(a[ix  ][iz+4] + a[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3] + a[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2] + a[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1] + a[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]) )*sz )*v[ix][iz]
#define LapT(a,ix,iz,sx,sz,v) ( ( C4*(a[ix+4][iz  ]*v[ix+4][iz  ] + a[ix-4][iz  ]*v[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ]*v[ix+3][iz  ] + a[ix-3][iz  ]*v[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ]*v[ix+2][iz  ] + a[ix-2][iz  ]*v[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ]*v[ix+1][iz  ] + a[ix-1][iz  ]*v[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sx                          +      \
                                ( C4*(a[ix  ][iz+4]*v[ix  ][iz+4] + a[ix  ][iz-4]*v[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3]*v[ix  ][iz+3] + a[ix  ][iz-3]*v[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2]*v[ix  ][iz+2] + a[ix  ][iz-2]*v[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1]*v[ix  ][iz+1] + a[ix  ][iz-1]*v[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sz )
/*^*/

typedef struct tri2 *tri2d;
/*^*/

struct tri2{
    bool verb,abc;
    int  nt, nx, nz, nb, depth, nxpad, nzpad, nzxpad;
    float dt2, idz2, idx2, cb;
};
/*^*/

#endif

/***************************************************************/
tri2d tri2d_make(bool verb, bool abc,
                int nt, int nx, int nz, int nb, int depth,
                float dt, float dx, float dz, float cb)
/*< initialize tri2d utilities >*/
{   
    tri2d tri;

    tri = (tri2d) sf_alloc(1,sizeof(*tri));

    tri->verb  = verb;
    tri->abc   = abc;
    tri->nt    = nt;
    tri->nx    = nx;
    tri->nz    = nz;
    tri->nb    = nb;
    tri->depth = depth;
    tri->dt2   = dt*dt;
    tri->idx2  = 1.0f/(dx*dx);
    tri->idz2  = 1.0f/(dz*dz);
    tri->cb    = cb; 

    tri->nxpad = nx+2*nb;
    tri->nzpad = nz+2*nb;
    tri->nzxpad= tri->nzpad*tri->nxpad;
    tri->depth = depth+nb;

    return tri;
}

/***************************************************************/
static tri2d tr;
static float **vvpad, **u0, **u1, **u2, **tmp;

/***************************************************************/
/* absorbing boundary */
static float *decay=NULL;

void abc_cal(int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< calculate absorbing coefficients >*/
{
    int ib;
    if(!nb) return;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
    for(ib=0; ib<nb; ib++){
        w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
    }
}

void abc_init(tri2d tri)
/*< initialization >*/
{
    if(tri->nb) decay =  sf_floatalloc(tri->nb);
    abc_cal(tri->nb,tri->cb,decay);
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(NULL!=decay) free(decay);
    decay = NULL;
}

void abc_apply(float *a /*2-D matrix*/,
               tri2d tri) 
/*< boundary decay>*/
{
    int iz, ix;

#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix)
{
#endif
    /* top & bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < tri->nb; iz++) {  
        for (ix=0; ix < tri->nxpad; ix++) {
	  a[tri->nzpad*ix +              iz] *= decay[iz];
	  a[tri->nzpad*ix + tri->nzpad-1-iz] *= decay[iz];
        }
    }
    /* left & right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < tri->nzpad; iz++) {  
        for (ix=0; ix < tri->nb; ix++) {
	  a[tri->nzpad*              ix  + iz] *= decay[ix];
	  a[tri->nzpad*(tri->nxpad-1-ix) + iz] *= decay[ix];
        }
    }
#ifdef _OPENMP
}
#endif
}


void timerev_init(bool verb, bool abc,
                  int nt, int nx, int nz, int nb, int depth,
                  float dt, float dx, float dz, float cb,
                  float **vv)
/*< initialize >*/
{
    int ix, iz;

    tr = tri2d_make(verb, abc, nt, nx, nz, nb, depth, dt, dx, dz, cb);

    /* set Laplacian coefficients */
    vvpad = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u0    = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u1    = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u2    = sf_floatalloc2(tr->nzpad, tr->nxpad);

    /* pad boundary */
    for     (ix=0; ix<tr->nx; ix++)
        for (iz=0; iz<tr->nz; iz++)
            vvpad[ix+tr->nb][iz+tr->nb] = vv[ix][iz]*vv[ix][iz]*tr->dt2;
    for     (ix=0; ix<tr->nxpad; ix++){
        for (iz=0; iz<tr->nb;    iz++){
            vvpad[ix][          iz  ] = vvpad[ix][          tr->nb  ];
            vvpad[ix][tr->nzpad-iz-1] = vvpad[ix][tr->nzpad-tr->nb-1];
        }
    }
    for     (ix=0; ix<tr->nb;    ix++){
        for (iz=0; iz<tr->nzpad; iz++){
            vvpad[          ix  ][iz]=vvpad[          tr->nb  ][iz];
            vvpad[tr->nxpad-ix-1][iz]=vvpad[tr->nxpad-tr->nb-1][iz];
        }
    }

    /* absorbing boundary condition */
    abc_init(tr);

}

void timerev_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< time reversal imaging linear operator >*/
{
    int ix, iz, it;
    float **dd, ***ww;

    if (nm!=tr->nz*tr->nx*tr->nt || nd!=tr->nt*tr->nx) sf_error("%s: wrong dimensions",__FILE__);
    sf_adjnull(adj, add, nm, nd, mod, dat);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nxpad; ix++)
                for (iz=0; iz<tr->nzpad; iz++)
                {
                    u0[ix][iz] = 0.0f;
                    u1[ix][iz] = 0.0f;
                    u2[ix][iz] = 0.0f;
                }
 
    /* map 1d to 2d */
    dd = (float**) sf_alloc (tr->nx,sizeof(float*)); 
    dd[0] = dat;
    for (ix=1; ix<tr->nx; ix++) dd[ix] = dd[0]+ix*tr->nt;

    /* map 1d to 3d */
    ww = (float***) sf_alloc (tr->nt,sizeof(float**));
    ww[0] = (float**) sf_alloc (tr->nx*tr->nt,sizeof(float*));
    ww[0][0] = mod;
    for (ix=1; ix<tr->nx*tr->nt; ix++) ww[0][ix] = ww[0][0]+ix*tr->nz; 
    for (it=1; it<tr->nt; it++) ww[it] = ww[0]+it*tr->nx;

    if (adj) { /* migration */
        
        for (it=tr->nt-1; it>-1; it--){
           
            /* 4 - apply abc */
            if (tr->abc) abc_apply(u1[0],tr);
            if (tr->abc) abc_apply(u0[0],tr);

            /* 3 - image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    ww[it][ix][iz] += u1[ix+tr->nb][iz+tr->nb]*vvpad[ix+tr->nb][iz+tr->nb];

            /* 2 - time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = LapT(u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;

            /* 1 - inject data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix=tr->nb; ix<tr->nb+tr->nx; ix++)
                u1[ix][tr->depth] += dd[ix-tr->nb][it];
 
        } /* it loop */

    } else { /* modeling */
    	
    	for (it=0; it<tr->nt; it++){

             /* 1 - record data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix=tr->nb; ix<tr->nb+tr->nx; ix++)
                dd[ix-tr->nb][it] += u1[ix][tr->depth];
           
            /* 2 - time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = Lap(u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;

            /* 3 - inject source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    u1[ix+tr->nb][iz+tr->nb] += ww[it][ix][iz]*vvpad[ix+tr->nb][iz+tr->nb];

            /* 4 - apply abc */
            if (tr->abc) abc_apply(u0[0],tr);
            if (tr->abc) abc_apply(u1[0],tr);

        } /* it loop */
        
    }

    free(dd);
    free(*ww); free(ww);
}

void timerev_close()
/*< finalize >*/
{
    abc_close();
    free(tr);
    free(*vvpad); free(vvpad);
    free(*u0); free(u0);
    free(*u1); free(u1);
    free(*u2); free(u2);
}

/***************************************************************/
void ctimerev(int ngrp, float ***ww, float **dd)
/*< correlative time reversal imaging condition >*/
{
    int ix, iz, it, ig, counter, *beg, *end;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
    for         (it=0; it<tr->nt; it++)
        for     (ix=0; ix<tr->nx; ix++)
            for (iz=0; iz<tr->nz; iz++)
                ww[it][ix][iz] = 1.0f;

    /* set start and end index */
    beg=sf_intalloc(ngrp);
    end=sf_intalloc(ngrp);
    counter = 0;
    for (ig=0; ig<ngrp; ig++) {
        beg[ig] = counter;
        counter += tr->nx/ngrp;
        end[ig] = counter;
    }
    end[ngrp-1] = tr->nx;
    if (tr->verb) {
        for (ig=0; ig<ngrp; ig++) {
            sf_warning("beg[%d]=%d",ig,beg[ig]);
            sf_warning("end[%d]=%d",ig,end[ig]);
        }
    }

    for (ig=0; ig<ngrp; ig++) { /* loop over subgroups of receivers */

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nxpad; ix++)
                for (iz=0; iz<tr->nzpad; iz++)
                {
                    u0[ix][iz] = 0.0f;
                    u1[ix][iz] = 0.0f;
                    u2[ix][iz] = 0.0f;
                }

        for (it=tr->nt-1; it>-1; it--){
            if (tr->verb) sf_warning("Time reversal: %d/%d;", it, 0);

            /* time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = Lap (u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;
            if (tr->abc) abc_apply(u1[0],tr);
            if (tr->abc) abc_apply(u0[0],tr);

            /* inject data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix=tr->nb+beg[ig]; ix<tr->nb+end[ig]; ix++)
                u1[ix][tr->depth] += dd[ix-tr->nb][it];

            /* image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    ww[it][ix][iz] *= u1[ix+tr->nb][iz+tr->nb];

        } /* it loop */
        if (tr->verb) sf_warning(".");

    } /* ig loop */


}

/***************************************************************/
/* sliding window normalization */

void threshold(bool step, int n, float hard, float *dat)
/*< in-place hard thresholding >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++) {
        if (dat[i]<hard) {
            if (step) dat[i] = 0.0f;
            else { 
                /*dat[i] /= hard;*/
                dat[i] = 0.5*(1.+cosf(SF_PI*(dat[i]/hard-1.))); /* Tukey window */
            }
        } else dat[i] = 1.0f;
    }
}

void absval(int n, float *dat)
/*< in-place absolute value >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++)
        dat[i] = fabs(dat[i]);
}

void autopow(int n, float p, float *dat)
/*< in-place auto-correlation with abs >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++)
        dat[i] = powf(fabs(dat[i]),p);
}

float maxval(int n, float *dat)
/*< maximum absolute value and variance (optional) >*/
{
    /* no parallelism */
    float dd, max=0;
    int i;

    for (i=0; i<n; i++) { 
        dd = dat[i];
        if (max<dd) max = dd;
    }

    return max;
}

void scale(float a, int n, float *dat)
/*< scale an array >*/
{
    /* no parallelism */
    int i;

    for (i=0; i<n; i++)
        dat[i] *= a;
}

void swnorm(bool verb, bool sw, int nz, int nx, int nt, int size, float perc, float *dat)
/*< local (sliding-window) normalization >*/
{
    int i, nzx, nzxt;
    float *dat0,den,factor;
    float max_all,pad;

    nzx = nz*nx;
    nzxt = nzx*nt;

    max_all = maxval(nzxt,dat);
    pad = max_all*perc/100.0f;
    if(verb) sf_warning("max_all=%g",max_all);

    dat0 = sf_floatalloc(nzxt);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<nzxt; i++) dat0[i] = dat[i];

    if (!sw) {

        scale(1./max_all,nzxt,dat);

    } else {

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=0; i<size; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(i+1+size),dat0);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=size; i<nt-size; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(2*size+1),dat0+(i-size)*nzx);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=nt-size; i<nt; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(size+1+(nt-1-i)),dat0+(i-size)*nzx);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }
        if (verb) sf_warning(".");

    }

    free(dat0);
}

/*********************************************************/
/* smoothing */
void smooth(int n1, int n2, int n3,
            int rect1, int rect2, int rect3,
            int nrep,
	    float *dat)
/*< Generate reflectivity map with smoothing >*/
{   
    int i, j, i0, irep, n123; 
    sf_triangle trig;
    int n[3],s[3],rect[3];

    n[0]=n1; n[1]=n2; n[2]=n3;
    s[0]=1;  s[1]=n1; s[2]=n1*n2;
    rect[0]=rect1; rect[1]=rect2; rect[2]=rect3;
    n123=n1*n2*n3;
    
    /* 2-d triangle smoothing */
    for (i=0;i<3;i++) {
        if (rect[i] <= 1) continue;
        trig = sf_triangle_init (rect[i],n[i],false);
        for (j=0; j < n123/n[i]; j++) {
            i0 = sf_first_index (i,j,3,n,s);
            for (irep=0; irep < nrep; irep++) {
                sf_smooth2 (trig,i0,s[i],false,dat);
            }
        }
        sf_triangle_close(trig);
    }
}


/*Qfwi_gradient*/
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
#if MPI_VERSION >= 2
	    sendbuf=MPI_IN_PLACE;
#else /* will fail */
	    sendbuf=NULL;
#endif
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
#if MPI_VERSION >= 2	    
		sendbuf=MPI_IN_PLACE;
#else /* will fail */
		sendbuf=NULL;
#endif
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
#if MPI_VERSION >= 2	    
		sendbuf=MPI_IN_PLACE;
#else /* will fail */
		sendbuf=NULL;
#endif
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
#if MPI_VERSION >= 2		    
		sendbuf=MPI_IN_PLACE;
#else /* will fail */
		sendbuf=NULL;
#endif		
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
        char filename[20]="tempbin",srdn[100000];
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
                if (nz*nx*nt != fread(wwt[0][0], sizeof(float), nz*nx*nt, temp))
		    sf_error ("%s: trouble reading:",__FILE__);
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
#if MPI_VERSION >= 2
		sendbuf=MPI_IN_PLACE;
#else /* will fail */
		sendbuf=NULL;
#endif		
		recvbuf=fcost;
	}else{
		sendbuf=fcost;
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);
	MPI_Bcast(fcost, 1, MPI_FLOAT, 0, comm);

	/* gradient reduction */
	if(cpuid==0){
#if MPI_VERSION >= 2	    
		sendbuf=MPI_IN_PLACE;
#else /* will fail */
		sendbuf=NULL;
#endif
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

void fwi(sf_file Fdat, sf_file Finv, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_optim optpar, bool verb, int media)
/*< fwi >*/
{
	int iter=0, flag;
	int nz, nx, nzx, nm=0;
	float fcost;
	float *x=NULL, *direction, *grad;
	sf_gradient gradient=NULL;
	FILE *fp=NULL;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
	if(fwipar->grad_type==1) {
		if(media==1) gradient=gradient_av;
		else gradient=gradient_v;
		nm=nzx;
		x=array->vv;
	}

	/* initialize */
	gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb);

	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	if(fwipar->onlygrad) return; // program terminates 

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");
	direction=sf_floatalloc(nm);
	optpar->sk=sf_floatalloc2(nm, optpar->npair);
	optpar->yk=sf_floatalloc2(nm, optpar->npair);

	optpar->igrad=0;
	optpar->ipair=0;
	optpar->alpha=1.;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
	if(mpipar->cpuid==0){
		l2norm(nm, grad, &optpar->gk_norm);
		print_iteration(fp, iter, optpar);
	}

	/* optimization loop */
	for(iter=0; iter<optpar->niter; iter++){
		if(mpipar->cpuid==0) sf_warning("--------iter=%d---------", iter);

		optpar->ils=0;

		//reverse(nm, grad, direction);
		if(iter==0){
			reverse(nm, grad, direction);
		}else{
			lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
			lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
		} 

		lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
		line_search(nm, x, grad, direction, gradient, optpar, &flag, mpipar->cpuid);
		
		if(mpipar->cpuid==0){
			l2norm(nm, grad, &optpar->gk_norm);
			print_iteration(fp, iter+1, optpar);
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

	if(mpipar->cpuid==0){
		sf_floatwrite(x, nm, Finv);
	}
	if(mpipar->cpuid==0) fclose(fp);
}

void lstri(sf_file Fdat, sf_file Fmwt, sf_file Fsrc, sf_mpi *mpipar, sf_acqui acpar, sf_vec array, sf_pas paspar, bool verb)
/*< passive source inversion >*/
{
    float **dd, ***ww, ***mwt;
    int nturn, iturn, is, rdn;
    char filename[20]="tempbin",srdn[20];
    FILE *temp;
    MPI_Comm comm=MPI_COMM_WORLD;

    if (mpipar->cpuid==0) {
        srand(time(NULL));
        rdn = rand()%1000000000;
        sprintf(srdn,"%d",rdn);
        strcat(filename,srdn);
    }
    MPI_Bcast(filename, 20, MPI_CHAR, 0, comm);
    if(verb && mpipar->cpuid==0) sf_warning("filename=%s",filename);

    temp=fopen(filename, "wb+");

    /*
#ifdef _OPENMP
#pragma omp parallel
    {
        sf_warning("id=%d, nthreads=%d",omp_get_thread_num(),omp_get_num_threads());
    }
#endif */
    dd = sf_floatalloc2(acpar->nt, acpar->nx);
    ww = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
    if (paspar->inv) mwt = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
    else mwt = NULL;

    if(acpar->ns%mpipar->numprocs==0) nturn=acpar->ns/mpipar->numprocs;
    else nturn=acpar->ns/mpipar->numprocs+1;

    for(iturn=0; iturn<nturn; iturn++){
        is=iturn*mpipar->numprocs+mpipar->cpuid;
        if(is<acpar->ns){
            if (paspar->inv) {
                /* read data */
                sf_seek(Fdat, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                sf_floatread(dd[0], acpar->nt*acpar->nx, Fdat);
            } else {
                /* read source */
                sf_seek(Fsrc, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                sf_floatread(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
            }

            /* do the computation */
            lstri_op(dd, NULL, ww, mwt, acpar, array, paspar, verb);

            if (paspar->inv) {
                /* write source */
                fseeko(temp, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                fwrite(ww[0][0], sizeof(float), acpar->nz*acpar->nx*acpar->nt, temp);
                if (NULL!=Fmwt && is==0) sf_floatwrite(mwt[0][0], acpar->nz*acpar->nx*acpar->nt, Fmwt);
            } else {
                /* write data */
                fseeko(temp, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                fwrite(dd[0], sizeof(float), acpar->nt*acpar->nx, temp);
            }

        } /* if is<ns */
    }
    fclose(temp);
    MPI_Barrier(comm);
    
    if(mpipar->cpuid==0) {
        temp=fopen(filename, "rb");
        if (paspar->inv) {
            for(is=0; is<acpar->ns; is++){
                fseeko(temp, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                if (!fread(ww[0][0], sizeof(float), acpar->nz*acpar->nx*acpar->nt, temp))
                    abort();
                sf_floatwrite(ww[0][0], acpar->nz * acpar->nx * acpar->nt, Fsrc);
            }
        } else {
            for(is=0; is<acpar->ns; is++){
                fseeko(temp, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                if (!fread(dd[0], sizeof(float), acpar->nt*acpar->nx, temp))
                    abort();
                sf_floatwrite(dd[0], acpar->nt * acpar->nx, Fdat);
            }
        }
        fclose(temp);
        remove(filename);
    }
    MPI_Barrier(comm);

    /* close */
    free(*dd); free(dd); 
    free(**ww); free(*ww); free(ww);
    if (paspar->inv) { free(**mwt); free(*mwt); free(mwt); }

}

void pfwi(sf_file Fdat, sf_file Finv, sf_file Fgrad, sf_file Fmwt, sf_file Fsrc, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_optim optpar, sf_pas paspar, bool verb)
/*< passive fwi >*/
{
	int iter=0, flag;
	int nz, nx, nzx, nm=0;
	float fcost;
	float *x=NULL, *direction, *grad;
	sf_gradient gradient=NULL;
	FILE *fp=NULL;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
        gradient=gradient_pas_av;
        nm=nzx;
        x=array->vv;

	/* initialize */
	gradient_pas_init(Fdat, Fsrc, Fmwt, mpipar, soupar, acpar, array, fwipar, paspar, verb);

	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	if(fwipar->onlygrad) return; /* program terminates */

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");
	direction=sf_floatalloc(nm);
	optpar->sk=sf_floatalloc2(nm, optpar->npair);
	optpar->yk=sf_floatalloc2(nm, optpar->npair);

	optpar->igrad=0;
	optpar->ipair=0;
	optpar->alpha=1.;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
        if(mpipar->cpuid==0) {
            l2norm(nm, grad, &optpar->gk_norm);
            print_iteration(fp, iter, optpar);
        }

	/* optimization loop */
        for(iter=0; iter<optpar->niter; iter++){
            if(mpipar->cpuid==0) sf_warning("--------iter=%d---------", iter);

            if (iter%optpar->repeat==0) optpar->alpha=1.;

            optpar->ils=0;

            if(iter==0){
                reverse(nm, grad, direction);
            }else{
                lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
                lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
            } 

            lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
            line_search(nm, x, grad, direction, gradient, optpar, &flag, mpipar->cpuid);
            if(mpipar->cpuid==0) {
                l2norm(nm, grad, &optpar->gk_norm);
                print_iteration(fp, iter+1, optpar);
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

	if(mpipar->cpuid==0){
            sf_floatwrite(x, nm, Finv);
        }

	if(mpipar->cpuid==0) fclose(fp);

}



void rtm_a(sf_file Fdat, sf_file Fimg, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, bool verb)
/*< acoustic rtm >*/
{
	int ix, iz, is, ir, it, wit;
	int sx, rx, sz, rz, rectx, rectz;
	int nz, nx, nzx, padnz, padnx, padnzx, nt, nr, nb, wnt;

	float dx2, dz2, dt2;
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
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						wave[wit][ix][iz]=p1[ix+nb][iz+nb];
				wit++;
			}

			/* laplacian operator */
			laplace(p1, term, padnx, padnz, dx2, dz2);
			
			/* load source */
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					term[ix][iz] += rr[ix*padnz+iz]*array->ww[it];
				}
			}

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

		/* check */
		if(wit != wnt) sf_error("Incorrect number of wavefield snapshots");
		/* read data */
		sf_seek(Fdat, is*nr*nt*sizeof(float), SEEK_SET);
		sf_floatread(dd[0], nr*nt, Fdat);
		/* initialization */
		memset(p0[0], 0., padnzx*sizeof(float));
		memset(p1[0], 0., padnzx*sizeof(float));
		memset(p2[0], 0., padnzx*sizeof(float));
		
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
			for(ix=0; ix<padnx; ix++){
				for(iz=0; iz<padnz; iz++){
					p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*vv[ix][iz]*dt2*term[ix][iz];
				}
			}

			/* calculate image */
			if(it%acpar->interval==0){
				for(ix=0; ix<nx; ix++)
					for(iz=0; iz<nz; iz++)
						mm[ix][iz] += wave[wit-1][ix][iz]*p1[ix+nb][iz+nb];
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



int main(int argc, char* argv[])
{
	int function, media, ntmp;
	bool verb;

	sf_mpi mpipar;
	sf_sou soupar;
	sf_acqui acpar;
	sf_vec array;
	sf_fwi fwipar = NULL;
	sf_optim optpar=NULL;
	sf_pas paspar=NULL;

	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Fq, Fw, Fdat, Fimg, Finv=NULL, Fgrad=NULL, Fsrc, Fmwt=NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &mpipar.cpuid);
	MPI_Comm_size(comm, &mpipar.numprocs);

	sf_init(argc, argv);

#ifdef _OPENMP
        omp_init();
#endif

	Fv=sf_input("Fvel"); /* velocity model */
	Fq=sf_input("Fq"); /* quality factor */
	Fw=sf_input("Fwavelet"); /* wavelet */

	soupar=(sf_sou)sf_alloc(1, sizeof(*soupar));
	acpar=(sf_acqui)sf_alloc(1, sizeof(*acpar));
	array=(sf_vec)sf_alloc(1, sizeof(*array));

	/* parameters I/O */
	if(!sf_getint("media", &media)) media=1;
	/* if 1, acoustic media; if 2, visco-acoustic media */
	if(!sf_getint("function", &function)) function=2;
	/* if 1, forward modeling; if 2, FWI; if 3, RTM */

	if(!sf_histint(Fv, "n1", &acpar->nz)) sf_error("No n1= in Fv");
	if(!sf_histint(Fv, "n2", &acpar->nx)) sf_error("No n2= in Fv");
	if(!sf_histfloat(Fv, "d1", &acpar->dz)) sf_error("No d1= in Fv");
	if(!sf_histfloat(Fv, "d2", &acpar->dx)) sf_error("No d2= in Fv");
	if(!sf_histfloat(Fv, "o1", &acpar->z0)) sf_error("No o1= in Fv");
	if(!sf_histfloat(Fv, "o2", &acpar->x0)) sf_error("No o2= in Fv");
	if(!sf_histint(Fw, "n1", &acpar->nt)) sf_error("No n1= in Fw");
	if(!sf_histfloat(Fw, "d1", &acpar->dt)) sf_error("No d1= in Fw");
	if(!sf_histfloat(Fw, "o1", &acpar->t0)) sf_error("No o1= in Fw");

	if(!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
	if(!sf_getint("nb", &acpar->nb)) acpar->nb=100; /* boundary width */
	if(!sf_getfloat("coef", &acpar->coef)) acpar->coef=0.003; /* absorbing boundary coefficient */

	if(!sf_getint("acqui_type", &acpar->acqui_type)) acpar->acqui_type=1;
	/* if 1, fixed acquisition; if 2, marine acquisition; if 3, symmetric acquisition */
	if(!sf_getint("ns", &acpar->ns)) sf_error("shot number required"); /* shot number */
	if(!sf_getfloat("ds", &acpar->ds)) sf_error("shot interval required"); /* shot interval */
	if(!sf_getfloat("s0", &acpar->s0)) sf_error("shot origin required"); /* shot origin */
	if(!sf_getint("sz", &acpar->sz)) acpar->sz=5; /* source depth */
	if(!sf_getint("nr", &acpar->nr)) acpar->nr=acpar->nx; /* number of receiver */
	if(!sf_getfloat("dr", &acpar->dr)) acpar->dr=acpar->dx; /* receiver interval */
	if(!sf_getfloat("r0", &acpar->r0)) acpar->r0=acpar->x0; /* receiver origin */
	if(!sf_getint("rz", &acpar->rz)) acpar->rz=1; /* receiver depth */

	if(!sf_getfloat("f0", &acpar->f0)) sf_error("reference frequency required"); /* reference frequency */
	if(!sf_getint("interval", &acpar->interval)) acpar->interval=1; /* wavefield storing interval */

	if(!sf_getfloat("fhi", &soupar->fhi)) soupar->fhi=0.5/acpar->dt; 
	if(!sf_getfloat("flo", &soupar->flo)) soupar->flo=0.; 
	soupar->rectx=2; 
	soupar->rectz=2; 

	/* get prepared */
	preparation(Fv, Fq, Fw, acpar, soupar, array);

        switch (function) {

            case 1: /* Modeling */

		Fdat=sf_output("output"); /* shot data */

		/* dimension set up */
		sf_putint(Fdat, "n1", acpar->nt);
		sf_putfloat(Fdat, "d1", acpar->dt);
		sf_putfloat(Fdat, "o1", acpar->t0);
		sf_putstring(Fdat, "label1", "Time");
		sf_putstring(Fdat, "unit1", "s");
		sf_putint(Fdat, "n2", acpar->nr);
		sf_putfloat(Fdat, "d2", acpar->dr);
		sf_putfloat(Fdat, "o2", acpar->r0);
		sf_putstring(Fdat, "label2", "Receiver");
		sf_putstring(Fdat, "unit2", "km");
		sf_putint(Fdat, "n3", acpar->ns);
		sf_putfloat(Fdat, "d3", acpar->ds);
		sf_putfloat(Fdat, "o3", acpar->s0);
		sf_putstring(Fdat, "label3", "Shot");
		sf_putstring(Fdat, "unit3", "km");

		if(media==1) forward_modeling_a(Fdat, &mpipar, soupar, acpar, array, verb);
		else forward_modeling(Fdat, &mpipar, soupar, acpar, array, verb);

		sf_fileclose(Fdat);

                break;

            case 2: /* FWI */

		fwipar=(sf_fwi)sf_alloc(1, sizeof(*fwipar));
		if(!sf_getbool("onlygrad", &fwipar->onlygrad)) fwipar->onlygrad=false; /* only want gradident */
		fwipar->grad_type=1;
		fwipar->misfit_type=1;
		fwipar->opt_type=1;
                if(!sf_getfloat("wt1", &fwipar->wt1)) fwipar->wt1=acpar->t0;
                if(!sf_getfloat("wt2", &fwipar->wt2)) fwipar->wt2=acpar->t0+(acpar->nt-1)*acpar->dt;
                if(!sf_getfloat("woff1", &fwipar->woff1)) fwipar->woff1=acpar->r0;
                if(!sf_getfloat("woff2", &fwipar->woff2)) fwipar->woff2=acpar->r0+(acpar->nr-1)*acpar->dr;
                if(!sf_getbool("oreo", &fwipar->oreo)) fwipar->oreo=false; /* keep oreo or keep cream */
		if(!sf_getint("waterz", &fwipar->waterz)) fwipar->waterz=51; /* water layer depth */
		if(!sf_getint("grectx", &fwipar->rectx)) fwipar->rectx=3; /* gradient smoothing radius in x */
		if(!sf_getint("grectz", &fwipar->rectz)) fwipar->rectz=3; /* gradient smoothing radius in z */

		Fdat=sf_input("Fdat"); /* input data */
		if(!fwipar->onlygrad) Finv=sf_output("output"); /* FWI result */
		Fgrad=sf_output("Fgrad"); /* FWI gradient at first iteration */

		/* dimension set up */
		if(Finv != NULL){
			sf_putint(Finv, "n1", acpar->nz);
			sf_putfloat(Finv, "d1", acpar->dz);
			sf_putfloat(Finv, "o1", acpar->z0);
			sf_putstring(Finv, "label1", "Depth");
			sf_putstring(Finv, "unit1", "km");
			sf_putint(Finv, "n2", acpar->nx);
			sf_putfloat(Finv, "d2", acpar->dx);
			sf_putfloat(Finv, "o2", acpar->x0);
			sf_putstring(Finv, "label2", "Distance");
			sf_putstring(Finv, "unit2", "km");
			if(fwipar->grad_type==3) sf_putint(Finv, "n3", 2);
		}
		sf_putint(Fgrad, "n1", acpar->nz);
		sf_putfloat(Fgrad, "d1", acpar->dz);
		sf_putfloat(Fgrad, "o1", acpar->z0);
		sf_putstring(Fgrad, "label1", "Depth");
		sf_putstring(Fgrad, "unit1", "km");
		sf_putint(Fgrad, "n2", acpar->nx);
		sf_putfloat(Fgrad, "d2", acpar->dx);
		sf_putfloat(Fgrad, "o2", acpar->x0);
		sf_putstring(Fgrad, "label2", "Distance");
		sf_putstring(Fgrad, "unit2", "km");
		if(fwipar->grad_type==3) sf_putint(Fgrad, "n3", 2);

		if(!fwipar->onlygrad){
			optpar=(sf_optim)sf_alloc(1, sizeof(*optpar));
			if(!sf_getint("niter", &optpar->niter)) sf_error("iteration number required"); /* iteration number */
			if(!sf_getfloat("conv_error", &optpar->conv_error)) sf_error("convergence error required"); /* final convergence error */
			optpar->npair=20; /* number of l-BFGS pairs */
			optpar->nls=20; /* line search number */
			if(!sf_getfloat("c1", &optpar->c1)) optpar->c1=1e-4;
			if(!sf_getfloat("c2", &optpar->c2)) optpar->c2=0.9;
			optpar->factor=10;
                        if(!sf_getfloat("v1", &optpar->v1)) optpar->v1=0.;
                        if(!sf_getfloat("v2", &optpar->v2)) optpar->v2=10.;
		}

		fwi(Fdat, Finv, Fgrad, &mpipar, soupar, acpar, array, fwipar, optpar, verb, media);

		if(!fwipar->onlygrad) sf_fileclose(Finv);
		sf_fileclose(Fgrad);

                break;

            case 3: /* RTM */

		Fdat=sf_input("Fdat"); /* input data */
		Fimg=sf_output("output"); /* rtm image */

		/* dimension set up */
		sf_putint(Fimg, "n1", acpar->nz);
		sf_putfloat(Fimg, "d1", acpar->dz);
		sf_putfloat(Fimg, "o1", acpar->z0);
		sf_putstring(Fimg, "label1", "Depth");
		sf_putstring(Fimg, "unit1", "km");
		sf_putint(Fimg, "n2", acpar->nx);
		sf_putfloat(Fimg, "d2", acpar->dx);
		sf_putfloat(Fimg, "o2", acpar->x0);
		sf_putstring(Fimg, "label2", "Distance");
		sf_putstring(Fimg, "unit2", "km");

		if(media==1) rtm_a(Fdat, Fimg, &mpipar, soupar, acpar, array, verb);
		else rtm(Fdat, Fimg, &mpipar, soupar, acpar, array, verb);

		sf_fileclose(Fimg);
                
                break;

            case 4: /* Passive FWI */

                paspar = passive_init(acpar);

                if (paspar->inv) {

                    Fdat=sf_input("Fdat");
                    if(!sf_histint(Fdat, "n3", &acpar->ns)) acpar->ns=1;

                    if (!paspar->onlysrc) {
                        fwipar=(sf_fwi)sf_alloc(1, sizeof(*fwipar));
                        if(!sf_getbool("onlygrad", &fwipar->onlygrad)) fwipar->onlygrad=false; /* only want gradident */
                        fwipar->grad_type=1;
                        fwipar->misfit_type=1;
                        fwipar->opt_type=1;
                        if(!sf_getfloat("wt1", &fwipar->wt1)) fwipar->wt1=acpar->t0;
                        if(!sf_getfloat("wt2", &fwipar->wt2)) fwipar->wt2=acpar->t0+(acpar->nt-1)*acpar->dt;
                        if(!sf_getfloat("woff1", &fwipar->woff1)) fwipar->woff1=acpar->r0;
                        if(!sf_getfloat("woff2", &fwipar->woff2)) fwipar->woff2=acpar->r0+(acpar->nr-1)*acpar->dr;
                        if(!sf_getbool("oreo", &fwipar->oreo)) fwipar->oreo=false; /* keep oreo or keep cream */
                        if(!sf_getint("waterz", &fwipar->waterz)) fwipar->waterz=0; /* water layer depth */
                        if(!sf_getint("waterzb", &fwipar->waterzb)) fwipar->waterzb=0; /* water layer depth from bottom up */
                        if(!sf_getint("grectx", &fwipar->rectx)) fwipar->rectx=3; /* gradient smoothing radius in x */
                        if(!sf_getint("grectz", &fwipar->rectz)) fwipar->rectz=3; /* gradient smoothing radius in z */

                        if(!fwipar->onlygrad) Finv=sf_output("output"); /* FWI result */
                        Fgrad=sf_output("Fgrad"); /* FWI gradient at first iteration */

                        /* dimension set up */
                        if(Finv != NULL){
                            sf_putint(Finv, "n1", acpar->nz);
                            sf_putfloat(Finv, "d1", acpar->dz);
                            sf_putfloat(Finv, "o1", acpar->z0);
                            sf_putstring(Finv, "label1", "Depth");
                            sf_putstring(Finv, "unit1", "km");
                            sf_putint(Finv, "n2", acpar->nx);
                            sf_putfloat(Finv, "d2", acpar->dx);
                            sf_putfloat(Finv, "o2", acpar->x0);
                            sf_putstring(Finv, "label2", "Distance");
                            sf_putstring(Finv, "unit2", "km");
                            /*if(fwipar->grad_type==3) sf_putint(Finv, "n3", 2);*/
                        }
                        sf_putint(Fgrad, "n1", acpar->nz);
                        sf_putfloat(Fgrad, "d1", acpar->dz);
                        sf_putfloat(Fgrad, "o1", acpar->z0);
                        sf_putstring(Fgrad, "label1", "Depth");
                        sf_putstring(Fgrad, "unit1", "km");
                        sf_putint(Fgrad, "n2", acpar->nx);
                        sf_putfloat(Fgrad, "d2", acpar->dx);
                        sf_putfloat(Fgrad, "o2", acpar->x0);
                        sf_putstring(Fgrad, "label2", "Distance");
                        sf_putstring(Fgrad, "unit2", "km");
                        /*if(fwipar->grad_type==3) sf_putint(Fgrad, "n3", 2);*/

                        if(!fwipar->onlygrad){
                            optpar=(sf_optim)sf_alloc(1, sizeof(*optpar));
                            if(!sf_getint("niter", &optpar->niter)) sf_error("iteration number required"); /* iteration number */
                            if(!sf_getint("repeat", &optpar->repeat)) optpar->repeat=1; /* repeat resetting alpha */
                            if(!sf_getfloat("conv_error", &optpar->conv_error)) sf_error("convergence error required"); /* final convergence error */
                            optpar->npair=20; /* number of l-BFGS pairs */
                            optpar->nls=20; /* line search number */
                            if(!sf_getfloat("c1", &optpar->c1)) optpar->c1=1e-4;
                            if(!sf_getfloat("c2", &optpar->c2)) optpar->c2=0.9;
                            optpar->factor=10;
                            if(!sf_getfloat("v1", &optpar->v1)) optpar->v1=0.;
                            if(!sf_getfloat("v2", &optpar->v2)) optpar->v2=10.;
                        }
                    } /* if !onlysrc */

                    if (!paspar->onlyvel) {
                        Fsrc=sf_output("Fsrc");
                        sf_putint   (Fsrc, "n1", acpar->nz);
                        sf_putfloat (Fsrc, "o1", acpar->z0);
                        sf_putfloat (Fsrc, "d1", acpar->dz);
                        sf_putstring(Fsrc, "label1", "Depth");
                        sf_putstring(Fsrc, "unit1" , "km");
                        sf_putint   (Fsrc, "n2", acpar->nx);
                        sf_putfloat (Fsrc, "o2", acpar->x0);
                        sf_putfloat (Fsrc, "d2", acpar->dx);
                        sf_putstring(Fsrc, "label2", "Distance");
                        sf_putstring(Fsrc, "unit2" , "km");
                        sf_putint   (Fsrc, "n3", acpar->nt);
                        sf_putfloat (Fsrc, "o3", acpar->t0);
                        sf_putfloat (Fsrc, "d3", acpar->dt);
                        sf_putstring(Fsrc, "label3", "Time");
                        sf_putstring(Fsrc, "unit3" , "s");
                        sf_putint   (Fsrc, "n4", acpar->ns);
                        sf_putfloat (Fsrc, "d4", 1.0f);
                        sf_putfloat (Fsrc, "o4", 0.0f);
                        sf_putstring(Fsrc, "label4", "Stage");

                        Fmwt=sf_output("Fmwt"); /* output data */
                        sf_putint   (Fmwt, "n1", acpar->nz);
                        sf_putfloat (Fmwt, "o1", acpar->z0);
                        sf_putfloat (Fmwt, "d1", acpar->dz);
                        sf_putstring(Fmwt, "label1", "Depth");
                        sf_putstring(Fmwt, "unit1" , "km");
                        sf_putint   (Fmwt, "n2", acpar->nx);
                        sf_putfloat (Fmwt, "o2", acpar->x0);
                        sf_putfloat (Fmwt, "d2", acpar->dx);
                        sf_putstring(Fmwt, "label2", "Distance");
                        sf_putstring(Fmwt, "unit2" , "km");
                        sf_putint   (Fmwt, "n3", acpar->nt);
                        sf_putfloat (Fmwt, "o3", acpar->t0);
                        sf_putfloat (Fmwt, "d3", acpar->dt);
                        sf_putstring(Fmwt, "label3", "Time");
                        sf_putstring(Fmwt, "unit3" , "s");
                    } else {
                        Fsrc=sf_input("Fsrc");
                        if(!sf_histint(Fsrc, "n4", &ntmp)) ntmp=1;
                        if (ntmp!=acpar->ns) sf_error("Shot dimension mismatch!");
                    }

                } else { /* modeling */
                    Fsrc=sf_input("Fsrc");
                    if(!sf_histint(Fsrc, "n4", &acpar->ns)) acpar->ns=1;

                    Fdat=sf_output("output"); /* output data */
                    sf_putint   (Fdat, "n1", acpar->nt);
                    sf_putfloat (Fdat, "o1", acpar->t0);
                    sf_putfloat (Fdat, "d1", acpar->dt);
                    sf_putstring(Fdat, "label1", "Time");
                    sf_putstring(Fdat, "unit1" , "s");
                    sf_putint   (Fdat, "n2", acpar->nx);
                    sf_putfloat (Fdat, "o2", acpar->x0);
                    sf_putfloat (Fdat, "d2", acpar->dx);
                    sf_putstring(Fdat, "label2", "Distance");
                    sf_putstring(Fdat, "unit2" , "km");
                    sf_putint   (Fdat, "n3", acpar->ns);
                    sf_putfloat (Fdat, "d3", 1.0f);
                    sf_putfloat (Fdat, "o3", 0.0f);
                    sf_putstring(Fdat, "label3", "Stage");
                }

                if (paspar->inv) {
                    if (paspar->onlysrc) { /* only inverting for source */
                        lstri(Fdat, Fmwt, Fsrc, &mpipar, acpar, array, paspar, verb);
                        sf_fileclose(Fsrc);
                        sf_fileclose(Fmwt);
                    } else { /* inverting for velocity ( and source ) */
                        pfwi(Fdat, Finv, Fgrad, Fmwt, Fsrc, &mpipar, soupar, acpar, array, fwipar, optpar, paspar, verb);
                        if(!fwipar->onlygrad) sf_fileclose(Finv);
                        sf_fileclose(Fgrad);
                        if (!paspar->onlyvel) {
                            sf_fileclose(Fsrc); 
                            sf_fileclose(Fmwt);
                        }
                    }
                } else {
                    lstri(Fdat, Fmwt, Fsrc, &mpipar, acpar, array, paspar, verb);
                    sf_fileclose(Fdat);
                }

                break;

            default:
                sf_warning("Please specify a valid function");

	} /* switch */

	MPI_Finalize();
	exit(0);
}
