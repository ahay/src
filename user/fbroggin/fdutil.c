#include <rsf.h>
#include "fdutil.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef _fdutil_h

typedef struct fdm2 *fdm2d;
/*^*/

typedef struct lcoef2 *lint2d;
/*^*/

typedef struct abc2 *abcone2d;
/*^*/

typedef struct spon *sponge;
/*^*/

struct fdm2{
	int nb;
	int   nz,nzpad;
	int   nx,nxpad;
	float oz,ozpad;
	float ox,oxpad;
	float dz;
	float dx;
	bool verb;
	bool free;
	int ompchunk;
	int nop;
};
/*^*/

struct lcoef2{
    int n;
    float *w00;
    float *w01;
    float *w10;
    float *w11;
    int *jz;
    int *jx;
};
/*^*/

struct abc2{
    bool free;
    float *bzl;
    float *bzh;
    float *bxl;
    float *bxh;
};
/*^*/

struct spon{
    float *w;
};
/*^*/

#endif

static float ** bell;

static int    nbell;

/*------------------------------------------------------------*/
fdm2d fdutil_init(	bool verb_, 
					bool free_,
					sf_axis az_, 
					sf_axis ax_, 
					int     nb_,
					int ompchunk_,
					int nop_) 
/*< init fdm utilities >*/
{ 
    fdm2d fdm;
    fdm = (fdm2d)sf_alloc(1,sizeof(*fdm));

    fdm->free=free_;
    fdm->verb=verb_;

    fdm->nb=nb_;

	fdm->nop=nop_;

    fdm->nz=sf_n(az_);
    fdm->nx=sf_n(ax_);

    fdm->dz=sf_d(az_);
    fdm->dx=sf_d(ax_);

    fdm->oz=sf_o(az_);
    fdm->ox=sf_o(ax_);

    fdm->nzpad=sf_n(az_)+2*fdm->nb+2*fdm->nop-1;
    fdm->nxpad=sf_n(ax_)+2*fdm->nb+2*fdm->nop-1;
	
    fdm->ozpad=sf_o(az_)-(fdm->nb+fdm->nop-1)*fdm->dz;
    fdm->oxpad=sf_o(ax_)-(fdm->nb+fdm->nop-1)*fdm->dx;

    fdm->ompchunk=ompchunk_;

    return fdm;
}

/*------------------------------------------------------------*/
void expand(float *a, 
			float *b, 
			fdm2d fdm)

/*< expand domain >*/
{
    int iz,ix;

	/* Copy the initial matrix a in the central part of the bigger matrix b */
	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(a,b)
	#endif
    for (ix=0; ix<fdm->nx; ix++) {
		for (iz=0; iz<fdm->nz; iz++) {
	    	b[(fdm->nb+fdm->nop-1+ix)*fdm->nzpad+(fdm->nb+fdm->nop-1+iz)] = a[ix*fdm->nz+iz];
		}
    }

	/* Initialize the top values of b */
	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(b)    
	#endif
    for (ix=0; ix<fdm->nxpad; ix++) {
		for (iz=0; iz<fdm->nb+fdm->nop-1; iz++) {
	    	b[ix*fdm->nzpad+iz] = b[ix*fdm->nzpad+fdm->nb+fdm->nop-1];
		}
    }

	/* Initialize the bottom values of b */
	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(b)    
	#endif
    for (ix=0; ix<fdm->nxpad; ix++) {
		for (iz=0; iz<fdm->nb+fdm->nop; iz++) {
	   		b[ix*fdm->nzpad+fdm->nzpad-iz-1] = b[ix*fdm->nzpad+fdm->nzpad-fdm->nb-fdm->nop-1];
		}
    }

	/* Initialize the left values of b */
	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(b)
	#endif
    for (ix=0; ix<fdm->nb+fdm->nop-1; ix++) {
		for (iz=0; iz<fdm->nzpad; iz++) {
	    	b[ix*fdm->nzpad+iz] = b[(fdm->nb+fdm->nop-1)*fdm->nzpad+iz];
		}
    }

	/* Initialize the right values of b */
	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(b)
	#endif
    for (ix=0; ix<fdm->nb+fdm->nop; ix++) {
		for (iz=0; iz<fdm->nzpad; iz++) {
	    	b[(fdm->nxpad-ix-1)*fdm->nzpad+iz] = b[(fdm->nxpad-fdm->nb-fdm->nop-1)*fdm->nzpad+iz];
		}
    }

}

/*------------------------------------------------------------*/
void cut2d( float *a,
			float *b,
			fdm2d  fdm,
			sf_axis cz, 
			sf_axis cx)
			
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix;
    int fz,fx;
    
    fx = (int)(floor)((sf_o(cx)-fdm->oxpad)/fdm->dx);
    fz = (int)(floor)((sf_o(cz)-fdm->ozpad)/fdm->dz);
    
    /*fprintf(stderr,"%d %d %d %d %d %d\n",fdm->nxpad,fdm->nzpad,sf_n(cx),sf_n(cz),fx,fz);*/

	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,iz) \
    shared(a,b)
	#endif
    for (ix=0; ix<sf_n(cx); ix++) {
		for (iz=0; iz<sf_n(cz); iz++) {
	    	b[ix*sf_n(cz)+iz] = a[(fx+ix)*fdm->nzpad+(fz+iz)];
		}
    }
}

/*------------------------------------------------------------*/
void bfill(float** b, 
	   fdm2d fdm)
/*< fill boundaries >*/
{
    int iz,ix;
    
    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<fdm->nb;    iz++) {
	    b[ix][           iz  ] = b[ix][           fdm->nb  ];
	    b[ix][fdm->nzpad-iz-1] = b[ix][fdm->nzpad-fdm->nb-1];
	}
    }
    
    for     (ix=0; ix<fdm->nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
	    b[           ix  ][iz] = b[           fdm->nb  ][iz];
	    b[fdm->nxpad-ix-1][iz] = b[fdm->nxpad-fdm->nb-1][iz];
	}
    }
}

/*------------------------------------------------------------*/
lint2d lint2d_make(	int na, 
					pt2d *aa, 
					fdm2d fdm)
		   
/*< init 2D linear interpolation >*/
{
    lint2d ca;
    int    ia;
    float f1,f2;
    
    ca = (lint2d) sf_alloc(1,sizeof(*ca));

    ca->n = na;

    ca->w00 = sf_floatalloc(na);
    ca->w01 = sf_floatalloc(na);
    ca->w10 = sf_floatalloc(na);
    ca->w11 = sf_floatalloc(na);

    ca->jz  = sf_intalloc(na);
    ca->jx  = sf_intalloc(na);

    for (ia=0; ia<na; ia++) {
	
		if(aa[ia].z >= fdm->ozpad && 
			aa[ia].z <  fdm->ozpad + (fdm->nzpad-1)*fdm->dz &&
			aa[ia].x >= fdm->oxpad && 
			aa[ia].x <  fdm->oxpad + (fdm->nxpad-1)*fdm->dx   ) {
				
			ca->jz[ia] = NINT((aa[ia].z-fdm->ozpad)/fdm->dz);
			ca->jx[ia] = NINT((aa[ia].x-fdm->oxpad)/fdm->dx);
				
			f1 = (aa[ia].z-fdm->ozpad)/fdm->dz - ca->jz[ia];
			f2 = (aa[ia].x-fdm->oxpad)/fdm->dx - ca->jx[ia];
		}
		else {
			ca->jz[ia] = 0; 
			ca->jx[ia] = 0;
			
			f1 = 1; 
			f2 = 0;
		}
		
		ca->w00[ia] = (1-f1)*(1-f2);
		ca->w01[ia] = (  f1)*(1-f2);
		ca->w10[ia] = (1-f1)*(  f2);
		ca->w11[ia] = (  f1)*(  f2);
		
    }
	
    return ca;
}

/*------------------------------------------------------------*/
void lint2d_hold(float**uu,
		 float *ww,
		 lint2d ca)
/*< hold fixed value in field >*/
{
    int   ia;
    float wa;

#ifdef _OPENMP
#pragma omp for private(ia,wa)
#endif
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] = wa;
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] = wa;
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] = wa;
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] = wa;
    }
}

/*------------------------------------------------------------*/
void lint2d_inject( float *uu,
					float *ww,
					float scal,
					lint2d ca,
				fdm2d fdm)
		   			
/*< inject into wavefield >*/
{
    int   ia;
    float wa;

#ifdef _OPENMP
#pragma omp for \
    private(ia)
#endif
	for (ia=0; ia<ca->n; ia++) {
		
		wa = ww[ia];
		
		uu[ca->jx[ia]*fdm->nzpad+ca->jz[ia] ] 		+= wa * scal * ca->w00[ia]; /* IMPORTANT - THE SOURCE IS ADDED (NOT SUBTRACTED LIKE IN PAUL SAVA'S CODE */
		uu[ca->jx[ia]*fdm->nzpad+(ca->jz[ia]+1)] 	+= wa * scal * ca->w01[ia];
		uu[(ca->jx[ia]+1)*fdm->nzpad+ca->jz[ia]]	+= wa * scal * ca->w10[ia];
		uu[(ca->jx[ia]+1)*fdm->nzpad+(ca->jz[ia]+1)]+= wa * scal * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void lint2d_extract(float *uu,
					float *dd,
					lint2d ca,
					fdm2d fdm)
					
/*< extract from wavefield >*/
{
    int ia;

	#ifdef _OPENMP
	#pragma omp parallel for \
	private(ia) \
	shared(dd)
	#endif
    for (ia=0; ia<ca->n; ia++) {
		dd[ia] =
	    	uu[ca->jx[ia]*fdm->nzpad+ca->jz[ia]        ] * 1.0; /*ca->w00[ia] +
	    	uu[ca->jx[ia]*fdm->nzpad+(ca->jz[ia]+1)    ] * ca->w01[ia] +
	    	uu[(ca->jx[ia]+1)*fdm->nzpad+ca->jz[ia]    ] * ca->w10[ia] +
	    	uu[(ca->jx[ia]+1)*fdm->nzpad+(ca->jz[ia]+1)] * ca->w11[ia];*/
    }
}  

/*------------------------------------------------------------*/
void fdbell_init(int n)
/*< init bell taper >*/
{
    int   iz,ix;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell=sf_floatalloc2(2*nbell+1,2*nbell+1);

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    bell[nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix)/s);
	}
    }    
}

/*------------------------------------------------------------*/
void lint2d_bell(float**uu,
		 float *ww,
		 lint2d ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix;
    float wa;

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    
	    for (ia=0;ia<ca->n;ia++) {
		wa = ww[ia] * bell[nbell+ix][nbell+iz];

		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w00[ia];
		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w01[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w10[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w11[ia];
	    }

	}
    }
}

/*------------------------------------------------------------*/
abcone2d abcone2d_make(	float	dt,
		       			float	*vv,
		       			bool	free, 
		       			fdm2d   fdm)
		       			
/*< init 2D ABC >*/

/*This absorbing boundary condition follows the work done in
Absorbing Boundary Conditions , by Robert Clayton and Bjorn Engquist

Which can be found at:

http://sepwww.stanford.edu/public/docs/sep11/11_12_abs.html
*/
{
    abcone2d abc;
    int iz,ix;
    float d;
    
    abc = (abcone2d) sf_alloc(1,sizeof(*abc));

    abc->free = free;

    abc->bzl = sf_floatalloc(fdm->nxpad);
    abc->bzh = sf_floatalloc(fdm->nxpad);
    abc->bxl = sf_floatalloc(fdm->nzpad);
    abc->bxh = sf_floatalloc(fdm->nzpad);

	#ifdef _OPENMP
	#pragma omp parallel for \
    private(ix,d)
	#endif
    for (ix=0; ix<fdm->nxpad; ix++) {
		
		d = vv[ix*fdm->nzpad+fdm->nop]*dt/fdm->dz;
		abc->bzl[ix] = (1-d)/(1+d);
		
		d = vv[ix*fdm->nzpad+fdm->nzpad-fdm->nop-1]*dt/fdm->dz;
		abc->bzh[ix] = (1-d)/(1+d);
    }

	#ifdef _OPENMP
	#pragma omp parallel for \
    private(iz,d)
	#endif
    for (iz=0; iz<fdm->nzpad; iz++) {
		
		d = vv[fdm->nop*fdm->nzpad+iz]*dt/fdm->dx;
		abc->bxl[iz] = (1-d)/(1+d);
		
		d = vv[(fdm->nxpad-fdm->nop-1)*fdm->nzpad+iz]*dt/fdm->dx;
		abc->bxh[iz] = (1-d)/(1+d);
    }

    return abc;
}

/*------------------------------------------------------------*/
void abcone2d_apply(float *a,
					float *b,
					abcone2d abc,
					fdm2d    fdm)
					
/*< apply 2D ABC >*/
{
    int iz,ix,iop;

	#ifdef _OPENMP
	#pragma omp for \
    private(iz,ix,iop)
	#endif
    for(ix=0; ix<fdm->nxpad; ix++) {
		for(iop=0; iop<fdm->nop; iop++) {
			
			/* top BC */
			if(!abc->free) { /* not free surface, apply ABC */
				iz = fdm->nop-iop;
				a[ix*fdm->nzpad+iz] = b[ix*fdm->nzpad+(iz+1)] + (b[ix*fdm->nzpad+iz] - a[ix*fdm->nzpad+(iz+1)])*abc->bzl[ix];
			}
			
			/* bottom BC */
			iz = fdm->nzpad-fdm->nop+iop-1;
			a[ix*fdm->nzpad+iz] = b[ix*fdm->nzpad+(iz-1)] + (b[ix*fdm->nzpad+iz] - a[ix*fdm->nzpad+(iz-1)])*abc->bzh[ix];
		}
	}

	#ifdef _OPENMP
	#pragma omp for \
    private(iz,ix,iop)
	#endif
    for(iz=0; iz<fdm->nzpad ;iz++) {
		for(iop=0; iop<fdm->nop; iop++) {
			
			/* left BC */
			ix = fdm->nop-iop;
			a[ix*fdm->nzpad+iz] = b[(ix+1)*fdm->nzpad+iz] + (b[ix*fdm->nzpad+iz] - a[(ix+1)*fdm->nzpad+iz])*abc->bxl[iz];
			
			/* right BC */
			ix = fdm->nxpad-fdm->nop+iop-1;
			a[ix*fdm->nzpad+iz] = b[(ix-1)*fdm->nzpad+iz] + (b[ix*fdm->nzpad+iz] - a[(ix-1)*fdm->nzpad+iz])*abc->bxh[iz];
		}
	}
}

/*------------------------------------------------------------*/
sponge sponge_make(fdm2d fdm)
/*< init boundary sponge >*/

/* Sponge boundary conditions multiply incoming wavefields
by smaller coefficients to attenuate the wavefield over time and space.

The sponge coefficients need to deviate from 1 very gradually to ensure
that there are no induced reflections caused by large impedance 
contrasts */
{
    sponge spo;
    int   ib,nbp;
    float sb,fb;
    
    nbp = fdm->nb+fdm->nop;
    
    spo = (sponge)sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(nbp);
    sb = 4.0*nbp;
    for(ib=0; ib<nbp; ib++) {
		fb = ib/(sqrt(2.0)*sb);
		spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

/*------------------------------------------------------------*/
void sponge2d_apply(float *a,
					float *b, 
					sponge   spo,
					fdm2d    fdm)
					
/*< apply boundary sponge >*/
{
    int iz,ix,ib;
    float w;

	/* Top and left sponges are nb+NOP samples. Bottom and right are nb+NOP-1 samples. */

	#ifdef _OPENMP
	#pragma omp for \
	private(ib,ix,w)
	#endif
	for(ix=0; ix<fdm->nxpad; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
		for(ib=0; ib<fdm->nb+fdm->nop; ib++) {
			w = spo->w[fdm->nb+fdm->nop-ib-1];
		    a[ix*fdm->nzpad+ib] *= w; /*    top sponge */
		    b[ix*fdm->nzpad+ib] *= w; /*    top sponge */
		}
	}
	
	#ifdef _OPENMP
	#pragma omp for \
    private(ib,ix,w)
	#endif
	for(ix=0; ix<fdm->nxpad; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
		for(ib=0; ib<fdm->nb+fdm->nop-1; ib++) {
			w = spo->w[fdm->nb+fdm->nop-ib-2];
		    a[ix*fdm->nzpad+fdm->nzpad-ib-1] *= w; /* bottom sponge */
		    b[ix*fdm->nzpad+fdm->nzpad-ib-1] *= w; /* bottom sponge */
		}
	}	

	#ifdef _OPENMP
	#pragma omp for \
    private(ib,iz,w)
	#endif
    for(ib=0; ib<fdm->nb+fdm->nop; ib++) {
		w = spo->w[fdm->nb+fdm->nop-ib-1];
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
		for(iz=0; iz<fdm->nzpad; iz++) {
		    a[ib*fdm->nzpad+iz] *= w; /*   left sponge */
		    b[ib*fdm->nzpad+iz] *= w; /*   left sponge */
		}

    }
    
	#ifdef _OPENMP
	#pragma omp for \
    private(ib,iz,w)
	#endif
    for(ib=0; ib<fdm->nb+fdm->nop-1; ib++) {
		w = spo->w[fdm->nb+fdm->nop-ib-2];
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
		for(iz=0; iz<fdm->nzpad; iz++) {
		    a[(fdm->nxpad-ib-1)*fdm->nzpad+iz] *= w; /*  right sponge */
		    b[(fdm->nxpad-ib-1)*fdm->nzpad+iz] *= w; /*  right sponge */
		}

    }    
    
}


