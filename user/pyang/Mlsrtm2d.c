/* Least-squares reverse time migration (LSRTM) with exact adjoint
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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

  Reference: Ji, Jun. "An exact adjoint operation pair in time extrapolation 
	and its application in least-squares reverse-time migration." 
	Geophysics 74.5 (2009): H27-H33.
*/
#include <rsf.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int nb, nz, nx, nt, nzpad, nxpad;
static float dz, dx, dt, c0, c11, c12, c21, c22;
static float *bndr;
static float **vv, **p0, **p1, **p2, **ptr=NULL;

void expand2d(float** b, float** a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz,ix;

    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    b[nb+ix][nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<nxpad; ix++) {
	for (iz=0; iz<nb;    iz++) {
	    b[ix][      iz  ] = b[ix][nb  ];
	    b[ix][nzpad-iz-1] = b[ix][nzpad-nb-1];
	}
    }

    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[ix 	 ][iz] = b[nb  		][iz];
	    b[nxpad-ix-1 ][iz] = b[nxpad-nb-1	][iz];
	}
    }
}


void window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz,ix;
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}
void rtm2d_adjnull(bool adj, bool add, float **mod, float **dat)
/*< rtm2d adjnull: mod=[0] if adj, else dat=[0] >*/
{
    if(add) return;
    
    if(adj) {
	memset(mod[0],0,nz*nx*sizeof(float));
    } else {
	memset(dat[0],0,nt*nx*sizeof(float));
    }
}

void fd2d_init(float d1, float d2)
/*< initialize 4-th order fd coefficients >*/
{
	float t;

	t = 1.0/(d1*d1);
	c11 = 4.0*t/3.0;
	c12= -t/12.0;

	t = 1.0/(d2*d2);
	c21 = 4.0*t/3.0;
	c22= -t/12.0;
	c0=-2.0*(c11+c12+c21 + c22);
}

void bndr_init()
/*< initialize sponge absorbing boundary coefficients >*/
{
	for(int i=0;i<nb;i++){
		float t=0.015*(nb-i);
		t=expf(-t*t);
		bndr[i]=powf(t,10.0);
	}
}

void rtm2d_init(int nb_, int nz_, int nx_, int nt_, 
	float dz_, float dx_, float dt_, float **v0)
/* allocate and initialize variables */
{
#ifdef _OPENMP
    	omp_init();
#endif
	nb=nb_;
	nz=nz_;
	nx=nx_;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nt=nt_;
	dz=dz_;
	dx=dx_;
	dt=dt_;
	
	bndr=sf_floatalloc(nb);
	vv=sf_floatalloc2(nzpad, nxpad);
	p0=sf_floatalloc2(nzpad, nxpad);
	p1=sf_floatalloc2(nzpad, nxpad);
	p2=sf_floatalloc2(nzpad, nxpad);

	bndr_init();
	expand2d(vv, v0);	
	for(int ix=0;ix<nxpad;ix++){
	    for(int iz=0;iz<nzpad;iz++){
		float tmp=vv[ix][iz]*dt;
		vv[ix][iz]=tmp*tmp;// velocity transformation
	    }
	}
	memset(p0[0],0,nzpad*nxpad*sizeof(float));
	memset(p1[0],0,nzpad*nxpad*sizeof(float));
	memset(p2[0],0,nzpad*nxpad*sizeof(float));
	fd2d_init(dz, dx);
}


void rtm2d_lop(bool adj, bool add, float **mod, float **dat)
/*< 2d rtm linear operator: 
modeling (adj==false); migration (adj==true)
>*/
{
	int iz, ix, it;
	rtm2d_adjnull(adj, false, mod, dat);

	if(!adj){/* forward modeling */
	    	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++)
		{
			p1[ix+nb][iz+nb]=mod[ix][iz];
		}

	    	for(ix=0; ix<nx; ix++) dat[ix][0]=p1[ix+nb][0+nb];

	   	for(it=1;it<nt;it++)
		{
			for (ix=0; ix < nxpad; ix++) 
			for (iz=0; iz < nzpad; iz++) 
			{
				float u = c0*p1[ix][iz];
				if(iz >= 1) u += c11*p1[ix][iz-1];
				if(iz >= 2) u += c12*p1[ix][iz-2];
				if(iz < nzpad-1) u += c11*p1[ix][iz+1];
				if(iz < nzpad-2) u += c12*p1[ix][iz+2];
				if(ix >= 1) u += c21*p1[ix-1][iz];
				if(ix >= 2) u += c22*p1[ix-2][iz];
				if(ix < nxpad-1) u += c21*p1[ix+1][iz];
				if(ix < nxpad-2) u += c22*p1[ix+2][iz];
				p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*u;
			}
			ptr=p0; p0=p1; p1=p2; p2=ptr;

			/* apply absorbing boundary condition */
			for(ix=0; ix<nxpad; ix++)
			{
				for(iz=0;iz<nb;iz++){	// top ABC			
					p0[ix][iz]=bndr[iz]*p0[ix][iz];
					p1[ix][iz]=bndr[iz]*p1[ix][iz];
				}	
				for(iz=nz;iz<nzpad;iz++){// bottom ABC			
					p0[ix][iz]=bndr[nzpad-iz-1]*p0[ix][iz];
					p1[ix][iz]=bndr[nzpad-iz-1]*p1[ix][iz];
				}	
			}
			for(iz=0; iz<nzpad; iz++)
			{
				for(ix=0;ix<nb;ix++){	// left ABC			
					p0[ix][iz]=bndr[ix]*p0[ix][iz];
					p1[ix][iz]=bndr[ix]*p1[ix][iz];
				}	
				for(ix=nx;ix<nxpad;ix++){// right ABC			
					p0[ix][iz]=bndr[nxpad-iz-1]*p0[ix][iz];
					p1[ix][iz]=bndr[nxpad-iz-1]*p1[ix][iz];
				}	
			}

			/* record seismogram */
			for(ix=0; ix<nx; ix++)	dat[ix][it]=p1[ix+nb][0+nb];
		}
	}else{/* migration, adjoint of modeling */
		for(ix=0; ix<nx; ix++)
		{
			p0[ix+nb][0+nb]=dat[ix][nt-1];
			p1[ix+nb][0+nb]=dat[ix][nt-2];
			p2[ix+nb][0+nb]=dat[ix][nt-3];
		}

		for(it=nt-2;it>=0;it--)
		{
			/* apply absorbing boundary condition */
			for(ix=0; ix<nxpad; ix++)
			{
				for(iz=0;iz<nb;iz++){	// top ABC			
					p0[ix][iz]=bndr[iz]*p0[ix][iz];
					p1[ix][iz]=bndr[iz]*p1[ix][iz];
				}	
				for(iz=nz;iz<nzpad;iz++){// bottom ABC			
					p0[ix][iz]=bndr[nzpad-iz-1]*p0[ix][iz];
					p1[ix][iz]=bndr[nzpad-iz-1]*p1[ix][iz];
				}	
			}
			for(iz=0; iz<nzpad; iz++)
			{
				for(ix=0;ix<nb;ix++){	// left ABC			
					p0[ix][iz]=bndr[ix]*p0[ix][iz];
					p1[ix][iz]=bndr[ix]*p1[ix][iz];
				}	
				for(ix=nx;ix<nxpad;ix++){// right ABC			
					p0[ix][iz]=bndr[nxpad-iz-1]*p0[ix][iz];
					p1[ix][iz]=bndr[nxpad-iz-1]*p1[ix][iz];
				}	
			}


			for(ix=0; ix<nxpad; ix++)
			for(iz=0; iz<nzpad-2; iz++)
			{
				p2[ix][iz]-=p0[ix][iz];
			}


			for(ix=0; ix<nxpad; ix++)
			for(iz=0; iz<nzpad; iz++)
			{
				float u = c0*vv[ix][iz]*p1[ix][iz];
				if(iz >= 1) u += c11*vv[ix][iz-1]*p1[ix][iz-1];
				if(iz >= 2) u += c12*vv[ix][iz-2]*p1[ix][iz-2];
				if(iz < nzpad-1) u += c11*vv[ix][iz+1]*p1[ix][iz+1];
				if(iz < nzpad-2) u += c12*vv[ix][iz+2]*p1[ix][iz+2];
				if(ix >= 1) u += c21*vv[ix-1][iz]*p1[ix-1][iz];
				if(ix >= 2) u += c22*vv[ix-2][iz]*p1[ix-2][iz];
				if(ix < nxpad-1) u += c21*vv[ix+1][iz]*p1[ix+1][iz];
				if(ix < nxpad-2) u += c22*vv[ix+2][iz]*p1[ix+2][iz];
				p2[ix][iz]=2*p1[ix][iz]+u;
/*
				p1[ix][iz]+=c22*(vv[ix-2][iz]*p1[ix-2][iz]+vv[ix+2][iz]*p1[ix+2][iz])+
					c21*(vv[ix-1][iz]*p1[ix-1][iz]+vv[ix+1][iz]*p1[ix+1][iz])+
					c20*vv[ix][iz]*p1[ix][iz]+
					c12*(vv[ix][iz-2]*p1[ix][iz-2]+vv[ix][iz+2]*p1[ix][iz+2])+
					c11*(vv[ix][iz-1]*p1[ix][iz-1]+vv[ix][iz+1]*p1[ix][iz+1])+
					c10*vv[ix][iz]*p1[ix][iz]+ 2.0*p1[ix][iz];
*/
			}
			ptr=p0; p0=p1; p1=p2; p2=ptr;

			/* inject data */
			memset(p2[0], 0, nzpad*nxpad*sizeof(float));
			if (it>1){
				for(ix=0; ix<nx; ix++) p2[ix+nb][0+nb]=dat[ix][it-2];
			}
		}
	}
}

void rtm2d_close()
/*< free the allocated memory >*/
{
	free(bndr);
	free(*vv); free(vv);	
	free(*p0); free(p0);
	free(*p1); free(p1);
	free(*p2); free(p2);
}

float inner_product2d(float **a, float **b, int n1, int n2)
/* compute inner product of vectorized 2d variables */
{
	float sum=0;
	for(int i2=0;i2<n2;i2++)
	for(int i1=0;i1<n1;i1++)
		sum+=a[i2][i1]*b[i2][i1];
	return sum;
}

void dot_test() 
/*< The dot product test to see if the adjoint is coded correctly.
   In the output dot1[0] shpould be equal to dot1[1] 
   (within machine precision),and dot2[0] should be equal to dot2[1]. 
>*/
{
    int nb_=30;
    int nz_=200;
    int nx_=200;
    int nt_=1800;
    float dz_=0.004;//unit=km
    float dx_=0.004;//unit=km
    float dt_=0.001;//unit=s
    float **v0=sf_floatalloc2(nz_, nx_);

    float *dot1, *dot2;
    float **mod1, **mod2, **dat1, **dat2;
    dot1 = sf_floatalloc(2);
    dot2 = sf_floatalloc(2);
    mod1 = sf_floatalloc2(nz,nx);
    mod2 = sf_floatalloc2(nz,nx);
    dat1 = sf_floatalloc2(nt,nx);
    dat2 = sf_floatalloc2(nt,nx);

    srand((int)time(0));
    for(int ix=0;ix<nx;ix++)
    { 
    	for(int iz=0;iz<nz;iz++)
    	{
	    v0[ix][iz]=500.0*((double)rand()/(double)RAND_MAX) +1500;
	    mod1[ix][iz]=rand();
    	}
    	for(int it=0;it<nt;it++)
    	{
	    dat2[ix][it]=rand();
    	}
    }
    rtm2d_init(nb_, nz_, nx_, nt_, dz_, dx_, dt_, v0);

    /* < L m1, d2 > = < m1, L* d2 > (w/o add) */
    rtm2d_lop(false, false, mod1, dat1);
    dot1[0] = inner_product2d(dat1, dat2, nt, nx);

    rtm2d_lop(true, false, mod2, dat2);
    dot1[1] = inner_product2d(mod1, mod2, nz, nx);

    /* < L m1, d2 > = < m1, L* d2 > (w/  add) */
    rtm2d_lop(false, true, mod1, dat1);
    dot2[0] = inner_product2d(dat1, dat2, nt, nx); 

    rtm2d_lop(true, true, mod2, dat2);
    dot2[1] = inner_product2d(mod1, mod2, nz, nx);

    /* dot-test report */ 
    fprintf(stderr,"\n < L m1, d2 > =?= < m1, L' d2 > \n");
    /* < L m1, d2 > = < m1, L' d2 > (w/o add) */
    fprintf(stderr,"%13.7f =?= %13.7f\n",dot1[0],dot1[1]);	
    /* < L m1, d2 > = < m1, L' d2 > (w/  add) */
    fprintf(stderr,"%13.7f =?= %13.7f\n",dot2[0],dot2[1]);


    free(dot1);
    free(dot2);
    free(*mod1); free(mod1);
    free(*mod2); free(mod2);
    free(*dat1); free(dat1);
    free(*dat2); free(dat2);
    rtm2d_close();
    free(*v0); free(v0);
}



void rtm2d_inversion(int nb_, int nz_, int nx_, int nt_, 
	float dz_, float dx_, float dt_, 
	float **v0, float **mod, float **dat, float tol, int niter)
/*< LSRTM with conjugate gradient method >*/
{
	int ix, iz, it, iter;
	float beta, alpha, g0, gn, gnp;
	float **rr, **mm, **gm, **gr, **sm;

	rr=sf_floatalloc2(nt_, nx_);
	gr=sf_floatalloc2(nt_, nx_);
	mm=sf_floatalloc2(nz_, nx_);
	gm=sf_floatalloc2(nz_, nx_);
	sm=sf_floatalloc2(nz_, nx_);

	for(ix=0;ix<nx_;ix++)
	for(it=0;it<nt_;it++)
	{
		rr[ix][it]=-dat[ix][it];
	}
	memset(gr[0],0,nt_*nx_*sizeof(float));
	memset(mm[0],0,nz_*nx_*sizeof(float));
	memset(gm[0],0,nz_*nx_*sizeof(float));
	memset(sm[0],0,nz_*nx_*sizeof(float));
    	rtm2d_init(nb_, nz_, nx_, nt_, dz_, dx_, dt_, v0);


	for(iter=0;iter<niter;iter++)
	{
		rtm2d_lop(true, false, rr, gm);
		gn=inner_product2d(gm, gm, nz_, nx_);
		if (iter==0){
			beta=0.0;
			g0=gn;
		}else{
			beta=gn/gnp;
			if(beta<tol || gn/g0<tol) break;
		}
		gnp=gn;

		for(ix=0;ix<nx_;ix++)
		for(iz=0;iz<nz_;iz++)
			sm[ix][iz]=gm[ix][iz]+beta*sm[ix][iz];
		rtm2d_lop(false, false, sm, gr);

		alpha=-gn/inner_product2d(gr, gr, nt_, nx_);

		for(ix=0; ix<nx_; ix++)
		{
			for(iz=0; iz<nz_; iz++) mm[ix][iz]+=alpha*sm[ix][iz];
			for(it=0; it<nt_; it++) rr[ix][it]+=alpha*gr[ix][iz];
		}
	}

	for(ix=0; ix<nx_; ix++)
	for(iz=0; iz<nz_; iz++)
		mod[ix][iz]=mm[ix][iz];

	free(*rr); free(rr);
	free(*gr); free(gr);
	free(*mm); free(mm);
	free(*gm); free(gm);
	free(*sm); free(sm);	

	rtm2d_close();
}

int main(int argc, char* argv[])
{
	int niter=50;
	float tol=1.e-6;
	int nb_=30;
	int nz_=200;
	int nx_=200;
	int nt_=1800;
	float dz_=0.004;//unit=km
	float dx_=0.004;//unit=km
	float dt_=0.001;//unit=s
	float **v0,**mod,**dat;
	v0=sf_floatalloc2(nz_, nx_);
	for(int ix=0; ix<nx_; ix++)
	for(int iz=0; iz<nz_; iz++)
	{
		v0[ix][iz]=2.0;//unit=km/s
	}

	mod=sf_floatalloc2(nz_,nx_);	
	dat=sf_floatalloc2(nt_,nx_);

	/* input dataset for dat */

	rtm2d_inversion(nb_, nz_, nx_, nt_, dz_, dx_, dt_, v0, mod, dat, tol, niter);


	free(*v0); free(v0);
	sf_close();
    	exit(0);
}

