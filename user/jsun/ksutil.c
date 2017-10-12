#include <rsf.h>
#include "ksutil.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

/* wall clock */
#include <sys/time.h>

/* fft */
#include "fftomp.h"

#ifndef _ksutil_h

/*------------------------------------------------------------*/
/* book keeping of staggered grid FD                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Muir's derivative operator */
/*#define C1 +0.598144 */
/*  1225/ 1024      /2 */
/*#define C2 -0.039876 */
/* -1225/(1024*  15)/2 */
/*#define C3 +0.004785 */
/*  1225/(1024* 125)/2 */
/*#define C4 -0.000348 */
/* -1225/(1024*1715)/2 */

/* forward FD derivative stencils */
/*#define Fx(a,ix,iy,iz,s) (C4*(a[iy  ][ix+4][iz  ] - a[iy  ][ix-3][iz  ]) + \*/
/*		            C3*(a[iy  ][ix+3][iz  ] - a[iy  ][ix-2][iz  ]) + \*/
/*		            C2*(a[iy  ][ix+2][iz  ] - a[iy  ][ix-1][iz  ]) + \*/
/*                          C1*(a[iy  ][ix+1][iz  ] - a[iy  ][ix  ][iz  ]) )*s*/
/*#define Fy(a,ix,iy,iz,s) (C4*(a[iy+4][ix  ][iz  ] - a[iy-3][ix  ][iz  ]) + \*/
/*		            C3*(a[iy+3][ix  ][iz  ] - a[iy-2][ix  ][iz  ]) + \*/
/*		            C2*(a[iy+2][ix  ][iz  ] - a[iy-1][ix  ][iz  ]) + \*/
/*                          C1*(a[iy+1][ix  ][iz  ] - a[iy  ][ix  ][iz  ]) )*s*/
/*#define Fz(a,ix,iy,iz,s) (C4*(a[iy  ][ix  ][iz+4] - a[iy  ][ix  ][iz-3]) + \*/
/*		            C3*(a[iy  ][ix  ][iz+3] - a[iy  ][ix  ][iz-2]) + \*/
/*		            C2*(a[iy  ][ix  ][iz+2] - a[iy  ][ix  ][iz-1]) + \*/
/*                          C1*(a[iy  ][ix  ][iz+1] - a[iy  ][ix  ][iz  ]) )*s*/

/* backward FD derivative stencils */
/*#define Bx(a,ix,iy,iz,s) (C4*(a[iy  ][ix+3][iz  ] - a[iy  ][ix-4][iz  ]) + \*/
/*		            C3*(a[iy  ][ix+2][iz  ] - a[iy  ][ix-3][iz  ]) + \*/
/*		            C2*(a[iy  ][ix+1][iz  ] - a[iy  ][ix-2][iz  ]) + \*/
/*                          C1*(a[iy  ][ix  ][iz  ] - a[iy  ][ix-1][iz  ]) )*s*/
/*#define By(a,ix,iy,iz,s) (C4*(a[iy+3][ix  ][iz  ] - a[iy-4][ix  ][iz  ]) + \*/
/*		            C3*(a[iy+2][ix  ][iz  ] - a[iy-3][ix  ][iz  ]) + \*/
/*		            C2*(a[iy+1][ix  ][iz  ] - a[iy-2][ix  ][iz  ]) + \*/
/*                          C1*(a[iy  ][ix  ][iz  ] - a[iy-1][ix  ][iz  ]) )*s*/
/*#define Bz(a,ix,iy,iz,s) (C4*(a[iy  ][ix  ][iz+3] - a[iy  ][ix  ][iz-4]) + \*/
/*		            C3*(a[iy  ][ix  ][iz+2] - a[iy  ][ix  ][iz-3]) + \*/
/*		            C2*(a[iy  ][ix  ][iz+1] - a[iy  ][ix  ][iz-2]) + \*/
/*                          C1*(a[iy  ][ix  ][iz  ] - a[iy  ][ix  ][iz-1]) )*s*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* CENTERED derivatives                                       */
/*------------------------------------------------------------*/
#define NOP 4 /* derivative operator half-size */
#define C1 4.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 4.0f/105.0f
#define C4 -1.0f/280.0f
#define Dx(a,ix,iy,iz,s) (C4*(a[iy][ix+4][iz] - a[iy][ix-4][iz]) +	\
			  C3*(a[iy][ix+3][iz] - a[iy][ix-3][iz]) +	\
			  C2*(a[iy][ix+2][iz] - a[iy][ix-2][iz]) +	\
			  C1*(a[iy][ix+1][iz] - a[iy][ix-1][iz])  )*s
#define Dy(a,ix,iy,iz,s) (C4*(a[iy+4][ix][iz] - a[iy-4][ix][iz]) +	\
			  C3*(a[iy+3][ix][iz] - a[iy-3][ix][iz]) +	\
			  C2*(a[iy+2][ix][iz] - a[iy-2][ix][iz]) +	\
			  C1*(a[iy+1][ix][iz] - a[iy-1][ix][iz])  )*s
#define Dz(a,ix,iy,iz,s) (C4*(a[iy][ix][iz+4] - a[iy][ix][iz-4]) +	\
			  C3*(a[iy][ix][iz+3] - a[iy][ix][iz-3]) +	\
			  C2*(a[iy][ix][iz+2] - a[iy][ix][iz-2]) +	\
			  C1*(a[iy][ix][iz+1] - a[iy][ix][iz-1])  )*s
/*^*/

#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define XZ 4
#define YZ 5
/*^*/

typedef struct fdm2 *fdm2d;
/*^*/

typedef struct fdm3 *fdm3d;
/*^*/

typedef struct lcoef2 *lint2d;
/*^*/

typedef struct lcoef3 *lint3d;
/*^*/

typedef struct abc2 *abcone2d;
/*^*/

typedef struct abc3 *abcone3d;
/*^*/

typedef struct spon *sponge;
/*^*/

typedef struct ofg *ofg2d;
/*^*/

typedef struct dft3 *dft3d;
/*^*/

typedef struct lps3 *lps3d;
/*^*/

typedef struct ksp3 *ksp3d;
/*^*/

typedef struct vksp3 *vksp3d;
/*^*/

typedef struct clr3 *clr3d;
/*^*/

typedef struct dat3 *dat3d;
/*^*/

typedef struct mut3 *mut3d;
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
};
/*^*/

struct fdm3{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    int   ny,nypad;
    float oz,ozpad;
    float ox,oxpad;
    float oy,oypad;
    float dz;
    float dx;
    float dy;
    bool verb;
    bool free;
    int ompchunk;
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

struct lcoef3{
    int n;
    float *w000;
    float *w001;
    float *w010;
    float *w011;
    float *w100;
    float *w101;
    float *w110;
    float *w111;
    int *jz;
    int *jx;
    int *jy;
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

struct abc3{
    bool free;
    float**bzl;
    float**bzh;
    float**bxl;
    float**bxh;
    float**byl;
    float**byh;
};
/*^*/

struct spon{
    float *w;
};
/*^*/

struct ofg{
    float **tt;
};
/*^*/

struct dft3{
    int nkz;
    int nkx;
    int nky;
    float okz;
    float okx;
    float oky;
    float dkz;
    float dkx;
    float dky;
};
/*^*/

struct lps3{
    float ***flt;
};
/*^*/

struct ksp3{
    float***kz;
    float***kx;
    float***ky;
};
/*^*/

struct vksp3{
    float***kpa;
    float***ksa;
    float***kpb;
    float***ksb;
};
/*^*/

struct clr3{
    int *n2_start;
    int *n2_local;
    int **map;
    int n2_sum;
    int n2_max;
};
/*^*/

struct dat3{
    int   ns,nr,nc;
    int   nt,nx,ny;
    float ot,ox,oy;
    float dt,dx,dy;
};
/*^*/

struct mut3{
    float **wt;
};
/*^*/

#endif

static float** bell;
static float***bell3d;
static int    nbell;

/*------------------------------------------------------------*/
fdm2d fdutil_init(bool verb_,
		  bool free_,
		  sf_axis az_,
		  sf_axis ax_,
		  int     nb_,
		  int ompchunk_)
/*< init fdm utilities >*/
{
    fdm2d fdm;
    fdm = (fdm2d) sf_alloc(1,sizeof(*fdm));

    fdm->free=free_;
    fdm->verb=verb_;

    fdm->nb=nb_;

    fdm->nz=sf_n(az_);
    fdm->nx=sf_n(ax_);

    fdm->dz=sf_d(az_);
    fdm->dx=sf_d(ax_);

    fdm->oz=sf_o(az_);
    fdm->ox=sf_o(ax_);

    fdm->nzpad=sf_n(az_)+2*fdm->nb;
    fdm->nxpad=sf_n(ax_)+2*fdm->nb;
	
    fdm->ozpad=sf_o(az_)-fdm->nb*fdm->dz;
    fdm->oxpad=sf_o(ax_)-fdm->nb*fdm->dx;

    fdm->ompchunk=ompchunk_;

    return fdm;
}

/*------------------------------------------------------------*/
fdm3d fdutil3d_init(bool verb_, 
		    bool free_,
		    sf_axis az_, 
		    sf_axis ax_, 
		    sf_axis ay_, 
		    int     nb_,
		    int ompchunk_) 
/*< init fdm utilities >*/
{ 
    fdm3d fdm;
    fdm = (fdm3d) sf_alloc(1,sizeof(*fdm));

    fdm->free=free_;
    fdm->verb=verb_;

    fdm->nb=nb_;
    

    fdm->nz=sf_n(az_);
    fdm->nx=sf_n(ax_);
    fdm->ny=sf_n(ay_);

    fdm->dz=sf_d(az_);
    fdm->dx=sf_d(ax_);
    fdm->dy=sf_d(ay_);

    fdm->oz=sf_o(az_);
    fdm->ox=sf_o(ax_);
    fdm->oy=sf_o(ay_);

    fdm->nzpad=sf_n(az_)+2*fdm->nb;
    fdm->nxpad=sf_n(ax_)+2*fdm->nb;
    if(sf_n(ay_)>1) fdm->nypad=sf_n(ay_)+2*fdm->nb;
    else fdm->nypad=1;
	
    fdm->ozpad=sf_o(az_)-fdm->nb*fdm->dz;
    fdm->oxpad=sf_o(ax_)-fdm->nb*fdm->dx;
    if(sf_n(ay_)>1) fdm->oypad=sf_o(ay_)-fdm->nb*fdm->dy;
    else fdm->oypad=sf_o(ay_);

    fdm->ompchunk=ompchunk_;

    return fdm;
}

/*------------------------------------------------------------*/
ofg2d offgrid_init(fdm2d fdm)
/*< init off-grid interpolation >*/
{
    ofg2d ofg;
    ofg = (ofg2d) sf_alloc(1,sizeof(*ofg));
    
    ofg->tt = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    
    return ofg;
}

/*------------------------------------------------------------*/
void offgridfor(float **ti,
		ofg2d  ofg,
		fdm2d  fdm)
/*< forward off-grid interpolation (in place) >*/
{
    int iz,ix;

    /* zero output */
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    ofg->tt[ix][iz]=0;
	}
    }

    for     (ix=0;ix<fdm->nx-1;ix++) {
	for (iz=0;iz<fdm->nz-1;iz++) {
	    ofg->tt[ix][iz] =
		ti[ix  ][iz  ] + 
		ti[ix+1][iz  ] + 
		ti[ix  ][iz+1] + 
		ti[ix+1][iz+1];
	}
    }   

    /* copy to input array */
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    ti[ix][iz] = 0.25 * ofg->tt[ix][iz];
	}
    }
}

/*------------------------------------------------------------*/
void offgridadj(float **ti,
		ofg2d  ofg,
		fdm2d  fdm)
/*< adjoint off-grid interpolation (in place) >*/
{
    int iz,ix;

    /* zero output */
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    ofg->tt[ix][iz]=0;
	}
    }
    
    for     (ix=0;ix<fdm->nx-1;ix++) {
	for (iz=0;iz<fdm->nz-1;iz++) {
	    ofg->tt[ix  ][iz  ] += ti[ix][iz];
	    ofg->tt[ix+1][iz  ] += ti[ix][iz];
	    ofg->tt[ix  ][iz+1] += ti[ix][iz];
	    ofg->tt[ix+1][iz+1] += ti[ix][iz];
	}
    }   

    /* copy to input array */
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    ti[ix][iz] = 0.25 * ofg->tt[ix][iz];
	}
    }
}

/*------------------------------------------------------------*/
void expand(float** a, 
	    float** b, 
	    fdm2d fdm)
/*< expand domain >*/
{
    int iz,ix;

    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    b[fdm->nb+ix][fdm->nb+iz] = a[ix][iz];
	}
    }

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
void expand3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm)
/*< expand domain >*/
{
    int iz,ix,iy;

    for         (iy=0;iy<fdm->ny;iy++) {
	for     (ix=0;ix<fdm->nx;ix++) {
	    for (iz=0;iz<fdm->nz;iz++) {
		b[fdm->nb+iy][fdm->nb+ix][fdm->nb+iz] = a[iy][ix][iz];
	    }
	}
    }

    for         (iy=0; iy<fdm->nypad; iy++) {
	for     (ix=0; ix<fdm->nxpad; ix++) {
	    for (iz=0; iz<fdm->nb;    iz++) {
		b[iy][ix][           iz  ] = b[iy][ix][           fdm->nb  ];
		b[iy][ix][fdm->nzpad-iz-1] = b[iy][ix][fdm->nzpad-fdm->nb-1];
	    }
	}
    }


    for         (iy=0; iy<fdm->nypad; iy++) {
	for     (ix=0; ix<fdm->nb;    ix++) {
	    for (iz=0; iz<fdm->nzpad; iz++) {
		b[iy][           ix  ][iz] = b[iy][           fdm->nb  ][iz];
		b[iy][fdm->nxpad-ix-1][iz] = b[iy][fdm->nxpad-fdm->nb-1][iz];
	    }
	}
    }

    for         (iy=0; iy<fdm->nb;    iy++) {
        for     (ix=0; ix<fdm->nxpad; ix++) {
            for (iz=0; iz<fdm->nzpad; iz++) {
                b[           iy  ][ix][iz] = b[           fdm->nb  ][ix][iz];
                b[fdm->nypad-iy-1][ix][iz] = b[fdm->nypad-fdm->nb-1][ix][iz];
            }
        }
    }

}

/*------------------------------------------------------------*/
void expand2d(float** a, 
              float** b, 
	      fdm3d fdm)
/*< expand domain >*/
{
    int iz,ix;

    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    b[fdm->nb+ix][fdm->nb+iz] = a[ix][iz];
	}
    }

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
void cut2d(float**  a,
	   float**  b,
	   fdm2d  fdm,
	   sf_axis cz, 
	   sf_axis cx)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix, nx,nz, jx,jz;
    int fz,fx;

    fz = (floor)((sf_o(cz)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(cx)-fdm->oxpad)/fdm->dx);

    nx = sf_n(cx);
    nz = sf_n(cz);

    jz = floor(sf_d(cz)/fdm->dz);
    jx = floor(sf_d(cx)/fdm->dx);
 
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iz)				\
    shared(a,b,nx,nz,fx,fz,jz,jx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    b[ix][iz] = a[fx+ix*jx][fz+iz*jz];
	}
    }
}

/*------------------------------------------------------------*/
void cut3d(float*** a,
	   float*** b,
	   fdm3d  fdm,
	   sf_axis cz, 
	   sf_axis cx,
	   sf_axis cy)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix,iy, nz,nx,ny, jz,jx,jy;
    int fz,fx,fy;

    fz = (floor)((sf_o(cz)-fdm->ozpad)/fdm->dz +0.0001f);
    fx = (floor)((sf_o(cx)-fdm->oxpad)/fdm->dx +0.0001f);
    fy = (floor)((sf_o(cy)-fdm->oypad)/fdm->dy +0.0001f);

    nz = sf_n(cz);
    nx = sf_n(cx);
    ny = sf_n(cy);

    jz = floor(sf_d(cz)/fdm->dz);
    jx = floor(sf_d(cx)/fdm->dx);
    jy = floor(sf_d(cy)/fdm->dy);
    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iy,iz)				\
    shared(a,b,nx,ny,nz,fx,fy,fz,jx,jy,jz)
#endif
    for         (iy=0;iy<ny;iy++) {
	for     (ix=0;ix<nx;ix++) {
	    for (iz=0;iz<nz;iz++) {
		b[iy][ix][iz] = a[fy+iy*jy][fx+ix*jx][fz+iz*jz];
	    }
	}
    }
}

/*------------------------------------------------------------*/
void cut3d_complex(sf_complex*** a,
	           sf_complex*** b,
                   fdm3d  fdm,
                   sf_axis cz, 
                   sf_axis cx,
                   sf_axis cy)
/*< cut a rectangular wavefield subset and convert to float >*/
{
    int iz,ix,iy, nz,nx,ny, jz,jx,jy;
    int fz,fx,fy;

    fz = (floor)((sf_o(cz)-fdm->ozpad)/fdm->dz +0.0001f);
    fx = (floor)((sf_o(cx)-fdm->oxpad)/fdm->dx +0.0001f);
    fy = (floor)((sf_o(cy)-fdm->oypad)/fdm->dy +0.0001f);

    nz = sf_n(cz);
    nx = sf_n(cx);
    ny = sf_n(cy);

    jz = floor(sf_d(cz)/fdm->dz);
    jx = floor(sf_d(cx)/fdm->dx);
    jy = floor(sf_d(cy)/fdm->dy);
    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iy,iz)				\
    shared(a,b,nx,ny,nz,fx,fy,fz,jx,jy,jz)
#endif
    for         (iy=0;iy<ny;iy++) {
	for     (ix=0;ix<nx;ix++) {
	    for (iz=0;iz<nz;iz++) {
		b[iy][ix][iz] = a[fy+iy*jy][fx+ix*jx][fz+iz*jz];
	    }
	}
    }
}

/*------------------------------------------------------------*/
void cut3d_slice(float** a,
                 float** b,
                 fdm3d  fdm,
                 sf_axis cz,
                 sf_axis cx)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix, nz,nx, jz,jx;
    int fz,fx;
    
    fz = (floor)((sf_o(cz)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(cx)-fdm->oxpad)/fdm->dx);
    
    nz = sf_n(cz);
    nx = sf_n(cx);
    
    jz = floor(sf_d(cz)/fdm->dz);
    jx = floor(sf_d(cx)/fdm->dx);
    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iz)				\
    shared(a,b,nx,nz,fx,fz,jx,jz)
#endif
    for     (ix=0;ix<nx;ix++) {
        for (iz=0;iz<nz;iz++) {
            b[ix][iz] = a[fx+ix*jx][fz+iz*jz];
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
lint2d lint2d_make(int    na, 
		   pt2d*  aa, 
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

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->ozpad && 
	   aa[ia].z <  fdm->ozpad + (fdm->nzpad-1)*fdm->dz &&
	   aa[ia].x >= fdm->oxpad && 
	   aa[ia].x <  fdm->oxpad + (fdm->nxpad-1)*fdm->dx   ) {
	    
	    ca->jz[ia] = (int)( (aa[ia].z-fdm->ozpad)/fdm->dz);
	    ca->jx[ia] = (int)( (aa[ia].x-fdm->oxpad)/fdm->dx);
	    
	    f1 = (aa[ia].z-fdm->ozpad)/fdm->dz - ca->jz[ia];
	    f2 = (aa[ia].x-fdm->oxpad)/fdm->dx - ca->jx[ia];
	} else {
	    ca->jz[ia] = 0; 
	    ca->jx[ia] = 0;
	    
	    f1 = 1; 
	    f2 = 0;
	}

	ca->w00[ia] = (1-f1)*(1-f2);
	ca->w01[ia] = (  f1)*(1-f2);
	ca->w10[ia] = (1-f1)*(  f2);
	ca->w11[ia] = (  f1)*(  f2);

	/* sf_warning("%g",ca->w00[ia]+ca->w01[ia]+ca->w10[ia]+ca->w11[ia]); */
    }

    return ca;
}

/*------------------------------------------------------------*/
lint3d lint3d_make(int    na, 
		   pt3d*  aa, 
		   fdm3d fdm)
/*< init 3D linear interpolation >*/
{
    lint3d ca;
    int    ia;
    float f1,f2,f3;
    
    ca = (lint3d) sf_alloc(1,sizeof(*ca));

    ca->n = na;

    ca->w000 = sf_floatalloc(na);
    ca->w001 = sf_floatalloc(na);
    ca->w010 = sf_floatalloc(na);
    ca->w011 = sf_floatalloc(na);
    ca->w100 = sf_floatalloc(na);
    ca->w101 = sf_floatalloc(na);
    ca->w110 = sf_floatalloc(na);
    ca->w111 = sf_floatalloc(na);

    ca->jz  = sf_intalloc(na);
    ca->jx  = sf_intalloc(na);
    ca->jy  = sf_intalloc(na);

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->ozpad && 
	   aa[ia].z <  fdm->ozpad + (fdm->nzpad-1)*fdm->dz &&
	   aa[ia].x >= fdm->oxpad && 
	   aa[ia].x <  fdm->oxpad + (fdm->nxpad-1)*fdm->dx &&   
	   aa[ia].y >= fdm->oypad && 
	   aa[ia].y <  fdm->oypad + (fdm->nypad-1)*fdm->dy  ) {
	    
	    ca->jz[ia] = (int)( (aa[ia].z-fdm->ozpad)/fdm->dz);
	    ca->jx[ia] = (int)( (aa[ia].x-fdm->oxpad)/fdm->dx);
	    ca->jy[ia] = (int)( (aa[ia].y-fdm->oypad)/fdm->dy);
	    
	    f1 = (aa[ia].z-fdm->ozpad)/fdm->dz - ca->jz[ia];
	    f2 = (aa[ia].x-fdm->oxpad)/fdm->dx - ca->jx[ia];
	    f3 = (aa[ia].y-fdm->oypad)/fdm->dy - ca->jy[ia];

	} else {
	    ca->jz[ia] = 0; 
	    ca->jx[ia] = 0;
	    ca->jy[ia] = 0;
	    
	    f1 = 1; 
	    f2 = 0;
	    f3 = 0;
	}

	ca->w000[ia] = (1-f3)*(1-f1)*(1-f2);
	ca->w001[ia] = (1-f3)*(  f1)*(1-f2);
	ca->w010[ia] = (1-f3)*(1-f1)*(  f2);
	ca->w011[ia] = (1-f3)*(  f1)*(  f2);

	ca->w100[ia] = (  f3)*(1-f1)*(1-f2);
	ca->w101[ia] = (  f3)*(  f1)*(1-f2);
	ca->w110[ia] = (  f3)*(1-f1)*(  f2);
	ca->w111[ia] = (  f3)*(  f1)*(  f2);
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
#pragma omp parallel for schedule(dynamic,1) private(ia,wa) shared(ca,ww,uu)
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
void lint2d_inject(float**uu,
		   float *dd,
		   lint2d ca)
/*< inject into wavefield >*/
{
    int   ia;
    for(ia=0;ia<ca->n;ia++) {	
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] += dd[ia] * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd[ia] * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd[ia] * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd[ia] * ca->w11[ia];
    }
}

void lint2d_extract(float**uu,
		    float* dd,
		    lint2d ca)
/*< extract from wavefield >*/
{
    int ia;
    for(ia=0;ia<ca->n;ia++) {
	dd[ia] =
	    uu[ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w00[ia] +
	    uu[ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w01[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w10[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w11[ia];
    }
}  

/*------------------------------------------------------------*/
void lint3d_inject(float***uu,
		   float  *dd,
		   lint3d  ca)
/*< inject into wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {	
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd[ia] * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd[ia] * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd[ia] * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd[ia] * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd[ia] * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd[ia] * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd[ia] * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd[ia] * ca->w111[ia];
    }
}

void lint3d_extract(float***uu,
		    float  *dd,
		    lint3d  ca)
/*< extract from wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {
	dd[ia] =
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w000[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w001[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w010[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w011[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w100[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w101[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w110[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w111[ia];
    }
}  

/*------------------------------------------------------------*/
void lint3d_inject_complex(sf_complex***uu,
                           sf_complex  *dd,
                           lint3d  ca)
/*< inject into wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {	
#ifdef SF_HAS_COMPLEX_H
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd[ia] * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd[ia] * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd[ia] * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd[ia] * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd[ia] * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd[ia] * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd[ia] * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd[ia] * ca->w111[ia];
#else
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] = sf_cadd(uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] , sf_crmul(dd[ia],ca->w000[ia]));
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] = sf_cadd(uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] , sf_crmul(dd[ia],ca->w001[ia]));
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] = sf_cadd(uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] , sf_crmul(dd[ia],ca->w010[ia]));
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] = sf_cadd(uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] , sf_crmul(dd[ia],ca->w011[ia]));
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] = sf_cadd(uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] , sf_crmul(dd[ia],ca->w100[ia]));
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] = sf_cadd(uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] , sf_crmul(dd[ia],ca->w101[ia]));
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] = sf_cadd(uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] , sf_crmul(dd[ia],ca->w110[ia]));
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] = sf_cadd(uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] , sf_crmul(dd[ia],ca->w111[ia]));
#endif
    }
}

void lint3d_extract_complex(sf_complex***uu,
                            sf_complex  *dd,
                            lint3d  ca)
/*< extract from wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {
#ifdef SF_HAS_COMPLEX_H
	dd[ia] =
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w000[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w001[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w010[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w011[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w100[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w101[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w110[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w111[ia];
#else
	dd[ia] =
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]  ] , ca->w000[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]+1] , ca->w001[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]  ] , ca->w010[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]+1] , ca->w011[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]  ] , ca->w100[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]+1] , ca->w101[ia]) ,
	    sf_cadd( sf_crmul(uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]  ] , ca->w110[ia]) ,
		     sf_crmul(uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]+1] , ca->w111[ia]))))))));
#endif
    }
}

/*------------------------------------------------------------*/
void lint3d_expl_complex(sf_complex***uz,
                         sf_complex***ux,
                         sf_complex***uy,
                         sf_complex **dd,
                         lint3d  ca)
/*< inject into wavefield >*/
{
    int ia,i,j,k;
    for (ia=0;ia<ca->n;ia++) {	
            for(k=-1;k<=1;k++)
                for(j=-1;j<=1;j++)
                    for(i=-1;i<=1;i++)
                    {
                        if(SF_ABS(i)+SF_ABS(j)+SF_ABS(k)==3)
                        {
#ifdef SF_HAS_COMPLEX_H
                            uz[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] += dd[0][ia] * i;
                            ux[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] += dd[1][ia] * j;
                            uy[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] += dd[2][ia] * k;
#else
                            uz[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] = sf_cadd(uz[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] , sf_crmul(dd[0][ia],(float)i));
                            ux[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] = sf_cadd(ux[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] , sf_crmul(dd[1][ia],(float)j));
                            uy[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] = sf_cadd(uy[ca->jy[ia]+k][ca->jx[ia]+j][ca->jz[ia]+i] , sf_crmul(dd[2][ia],(float)k));


#endif
                        }
                    }
    }
}

/*------------------------------------------------------------*/
void lint2d_inject1(float**uu,
		    float  dd,
		    lint2d ca)
/*< inject into wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] += dd * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void lint3d_inject1(float***uu,
		    float   dd,
		    lint3d  ca)
/*< inject into wavefield >*/
{
    int ia;
    for (ia=0;ia<ca->n;ia++) {
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] += dd * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] += dd * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] += dd * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] += dd * ca->w111[ia];
    }
}

/*------------------------------------------------------------*/
void fdbell_init(int n)
/*< init bell taper >*/
{
    int   iz,ix;
    float s;

    nbell = n;
    s = nbell==0?1:2.0/(nbell*nbell);

    bell=sf_floatalloc2(2*nbell+1,2*nbell+1);

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    bell[nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix)*s);
	}
    }    
}

/*------------------------------------------------------------*/
void fdbell3d_init(int n)
/*< init bell taper >*/
{
    int   iz,ix,iy;
    float s;

    nbell = n;
    s = nbell==0?1:2.0/(nbell*nbell);

    bell3d=sf_floatalloc3(2*nbell+1,2*nbell+1,2*nbell+1);

    for        (iy=-nbell;iy<=nbell;iy++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		bell3d[nbell+iy][nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix+iy*iy)*s);
	    }
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

		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w00[ia];
		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w01[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w10[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w11[ia];
	    }

	}
    }
}

/*------------------------------------------------------------*/
void lint3d_bell(float***uu,
		 float  *ww,
		 lint3d  ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix,iy;
    float wa;

    for        (iy=-nbell;iy<=nbell;iy++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		
		for (ia=0;ia<ca->n;ia++) {
		    wa = ww[ia] * bell3d[nbell+iy][nbell+ix][nbell+iz];
		    
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w000[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w001[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w010[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w011[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w100[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w101[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w110[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w111[ia];
		}
		
	    }
	}
    }
}

/*------------------------------------------------------------*/
void lint3d_bell_complex(sf_complex***uu,
                         sf_complex  *ww,
                         lint3d  ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix,iy;
    sf_complex wa;

    for        (iy=-nbell;iy<=nbell;iy++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		
		for (ia=0;ia<ca->n;ia++) {
#ifdef SF_HAS_COMPLEX_H
		    wa = ww[ia] * bell3d[nbell+iy][nbell+ix][nbell+iz];
		    
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w000[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w001[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w010[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w011[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w100[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w101[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w110[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w111[ia];
#else
		    wa = sf_crmul(ww[ia] , bell3d[nbell+iy][nbell+ix][nbell+iz]);
		    
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] = sf_cadd(uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ],sf_crmul(wa , ca->w000[ia]));
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] = sf_cadd(uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ],sf_crmul(wa , ca->w001[ia]));
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] = sf_cadd(uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ],sf_crmul(wa , ca->w010[ia]));
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] = sf_cadd(uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ],sf_crmul(wa , ca->w011[ia]));
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] = sf_cadd(uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ],sf_crmul(wa , ca->w100[ia]));
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] = sf_cadd(uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ],sf_crmul(wa , ca->w101[ia]));
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] = sf_cadd(uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ],sf_crmul(wa , ca->w110[ia]));
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] = sf_cadd(uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ],sf_crmul(wa , ca->w111[ia]));
#endif
                }

            }
        }
    }
}

/*------------------------------------------------------------*/
void lint2d_bell1(float**uu,
		  float ww,
		  lint2d ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix;
    float wa;

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    
	    for (ia=0;ia<ca->n;ia++) {
		wa = ww * bell[nbell+ix][nbell+iz];

		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w00[ia];
		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w01[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w10[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w11[ia];
	    }

	}
    }
}

/*------------------------------------------------------------*/
void lint3d_bell1(float***uu,
		  float  ww,
		  lint3d  ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix,iy;
    float wa;

    for        (iy=-nbell;iy<=nbell;iy++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		wa = ww * bell3d[nbell+iy][nbell+ix][nbell+iz];

		for (ia=0;ia<ca->n;ia++) {
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w000[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w001[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w010[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w011[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] += wa * ca->w100[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] += wa * ca->w101[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] += wa * ca->w110[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] += wa * ca->w111[ia];
		}
		
	    }
	}
    }
}

/*------------------------------------------------------------*/
abcone2d abcone2d_make(int     nop,
		       float    dt,
		       float**  vv,
		       bool   free, 
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

    for (ix=0;ix<fdm->nxpad;ix++) {
	d = vv[ix][           nop  ] *dt/fdm->dz; abc->bzl[ix] = (1-d)/(1+d);
	d = vv[ix][fdm->nzpad-nop-1] *dt/fdm->dz; abc->bzh[ix] = (1-d)/(1+d);
    }
    for (iz=0;iz<fdm->nzpad;iz++) {
	d = vv[           nop  ][iz] *dt/fdm->dx; abc->bxl[iz] = (1-d)/(1+d);
	d = vv[fdm->nxpad-nop-1][iz] *dt/fdm->dx; abc->bxh[iz] = (1-d)/(1+d);
    }

    return abc;
}

/*------------------------------------------------------------*/
abcone3d abcone3d_make(int     nop,
		       float    dt,
		       float ***vv,
		       bool   free, 
		       fdm3d   fdm)
/*< init 3D ABC >*/

/*This absorbing boundary condition follows the work done in
  Absorbing Boundary Conditions , by Robert Clayton and Bjorn Engquist

  Which can be found at:

  http://sepwww.stanford.edu/public/docs/sep11/11_12_abs.html
*/

{
    abcone3d abc;
    int iz,ix,iy;
    float d;

    abc = (abcone3d) sf_alloc(1,sizeof(*abc));

    abc->free = free;

    /* z */
    abc->bzl = sf_floatalloc2(fdm->nxpad,fdm->nypad);
    abc->bzh = sf_floatalloc2(fdm->nxpad,fdm->nypad);

    /* x */
    abc->bxl = sf_floatalloc2(fdm->nzpad,fdm->nypad);
    abc->bxh = sf_floatalloc2(fdm->nzpad,fdm->nypad);

    /* y */
    abc->byl = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    abc->byh = sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for     (iy=0;iy<fdm->nypad;iy++) {
	for (ix=0;ix<fdm->nxpad;ix++) {
	    d = vv[iy][ix][           nop  ] *dt/fdm->dz; abc->bzl[iy][ix] = (1-d)/(1+d);
	    d = vv[iy][ix][fdm->nzpad-nop-1] *dt/fdm->dz; abc->bzh[iy][ix] = (1-d)/(1+d);
	}
    }

    for     (iy=0;iy<fdm->nypad;iy++) {
	for (iz=0;iz<fdm->nzpad;iz++) {
	    d = vv[iy][           nop  ][iz] *dt/fdm->dx; abc->bxl[iy][iz] = (1-d)/(1+d);
	    d = vv[iy][fdm->nxpad-nop-1][iz] *dt/fdm->dx; abc->bxh[iy][iz] = (1-d)/(1+d);
	}
    }
    
    for     (ix=0;ix<fdm->nxpad;ix++) {
	for (iz=0;iz<fdm->nzpad;iz++) {
	    d = vv[           nop  ][ix][iz] *dt/fdm->dy; abc->byl[ix][iz] = (1-d)/(1+d);
	    d = vv[fdm->nypad-nop-1][ix][iz] *dt/fdm->dy; abc->byh[ix][iz] = (1-d)/(1+d);
	}
    }

    return abc;
}


/*------------------------------------------------------------*/
void abcone2d_apply(float**   uo,
		    float**   um,
		    int      nop,
		    abcone2d abc,
		    fdm2d    fdm)
/*< apply 2D ABC >*/
{
    int iz,ix,iop;

/*#ifdef _OPENMP*/
/*#pragma omp parallel for		\*/
/*    schedule(dynamic)			\*/
/*    private(iz,ix,iop)			\*/
/*    shared(fdm,abc,nop,uo,um)*/
/*#endif*/
	for(ix=0;ix<fdm->nxpad;ix++) {
	    for(iop=0;iop<nop;iop++) {

		/* top BC */
		if(!abc->free) { /* not free surface, apply ABC */
		    iz = nop-iop;
		    uo      [ix][iz  ] 
			= um[ix][iz+1] 
			+(um[ix][iz  ]
			  - uo[ix][iz+1]) * abc->bzl[ix];
		}

		/* bottom BC */
		iz = fdm->nzpad-nop+iop-1;
		uo      [ix][iz  ] 
		    = um[ix][iz-1]
		    +(um[ix][iz  ]
		      - uo[ix][iz-1]) * abc->bzh[ix];
	    }
	}

/*#ifdef _OPENMP*/
/*#pragma omp parallel for		\*/
/*    schedule(dynamic)			\*/
/*    private(iz,ix,iop)			\*/
/*    shared(fdm,nop,uo,um,abc)*/
/*#endif*/
	    for(iz=0;iz<fdm->nzpad;iz++) {
		for(iop=0;iop<nop;iop++) {

		    /* left BC */
		    ix = nop-iop;
		    uo      [ix  ][iz] 
			= um[ix+1][iz] 
			+(um[ix  ][iz]
			  - uo[ix+1][iz]) * abc->bxl[iz];

		    /* right BC */
		    ix = fdm->nxpad-nop+iop-1;
		    uo      [ix  ][iz] 
			= um[ix-1][iz]
			+(um[ix  ][iz]
			  - uo[ix-1][iz]) * abc->bxh[iz];
		}
	    }
}

/*------------------------------------------------------------*/
void abcone3d_apply(float  ***uo,
		    float  ***um,
		    int      nop,
		    abcone3d abc,
		    fdm3d    fdm)
/*< apply 3D ABC >*/
{
    int iz,ix,iy,iop;

/*#ifdef _OPENMP*/
/*#pragma omp parallel for			\*/
/*    schedule(dynamic)				\*/
/*    private(iz,ix,iy,iop)			\*/
/*    shared(fdm,abc,nop,uo,um)*/
/*#endif*/
	for    (iy=0;iy<fdm->nypad;iy++) {
	    for(ix=0;ix<fdm->nxpad;ix++) {
		for(iop=0;iop<nop;iop++) {
		
		    /* z min */
		    if(!abc->free) { /* not free surface, apply ABC */
			iz = nop-iop;
			uo      [iy][ix][iz  ] 
			    = um[iy][ix][iz+1] 
			    +(um[iy][ix][iz  ]
			      - uo[iy][ix][iz+1]) * abc->bzl[iy][ix];
		    }
		
		    /* z max */
		    iz = fdm->nzpad-nop+iop-1;
		    uo      [iy][ix][iz  ] 
			= um[iy][ix][iz-1]
			+(um[iy][ix][iz  ]
			  - uo[iy][ix][iz-1]) * abc->bzh[iy][ix];
		}
	    }
	}
    
/*#ifdef _OPENMP*/
/*#pragma omp parallel for			\*/
/*    schedule(dynamic)				\*/
/*    private(iz,ix,iy,iop)			\*/
/*    shared(fdm,abc,nop,uo,um)*/
/*#endif*/
	    for    (iy=0;iy<fdm->nypad;iy++) {
		for(iz=0;iz<fdm->nzpad;iz++) {
		    for(iop=0;iop<nop;iop++) {
		
			/* x min */
			ix = nop-iop;
			uo      [iy][ix  ][iz] 
			    = um[iy][ix+1][iz] 
			    +(um[iy][ix  ][iz]
			      - uo[iy][ix+1][iz]) * abc->bxl[iy][iz];
		
			/* x max */
			ix = fdm->nxpad-nop+iop-1;
			uo      [iy][ix  ][iz] 
			    = um[iy][ix-1][iz]
			    +(um[iy][ix  ][iz]
			      - uo[iy][ix-1][iz]) * abc->bxh[iy][iz];
		    }
		}
	    }

/*#ifdef _OPENMP*/
/*#pragma omp parallel for			\*/
/*    schedule(dynamic)				\*/
/*    private(iz,ix,iy,iop)			\*/
/*    shared(fdm,abc,nop,uo,um)*/
/*#endif*/
		for    (ix=0;ix<fdm->nxpad;ix++) {
		    for(iz=0;iz<fdm->nzpad;iz++) {
			for(iop=0;iop<nop;iop++) {
		
			    /* y min */
			    iy = nop-iop;
			    uo      [iy  ][ix][iz] 
				= um[iy+1][ix][iz] 
				+(um[iy  ][ix][iz]
				  - uo[iy+1][ix][iz]) * abc->byl[ix][iz];
		
			    /* y max */
			    iy = fdm->nypad-nop+iop-1;
			    uo      [iy  ][ix][iz] 
				= um[iy-1][ix][iz]
				+(um[iy  ][ix][iz]
				  - uo[iy-1][ix][iz]) * abc->byh[ix][iz];
			}
		    }
		}

}


/*------------------------------------------------------------*/
sponge sponge_make(int nb)
/*< init boundary sponge >*/
    
/* Sponge boundary conditions multiply incoming wavefields
   by smaller coefficients to attenuate the wavefield over time and space.

   The sponge coefficients need to deviate from 1 very gradually to ensure
   that there are no induced reflections caused by large impedance 
   contrasts */
{
    sponge spo;
    int   ib;
    float sb,fb;
    
    spo = (sponge) sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(nb);
    sb = 4.0*nb;               
    for(ib=0; ib<nb; ib++) {
	fb = ib/(sqrt(2.0)*sb);
	spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

sponge sponge_make2(int nb, float cb)
/*< init boundary sponge >*/
    
/* Sponge boundary conditions multiply incoming wavefields
   by smaller coefficients to attenuate the wavefield over time and space.

   The sponge coefficients need to deviate from 1 very gradually to ensure
   that there are no induced reflections caused by large impedance 
   contrasts */
{
    sponge spo;
    int   ib;
    float sb,fb;
    
    spo = (sponge) sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(nb);
    sb = 4.0*nb;               
    for(ib=0; ib<nb; ib++) {
	fb = ib/(sqrt(2.0)*sb)*cb;
	spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

/*------------------------------------------------------------*/
void sponge2d_apply(float**   uu,
		    sponge   spo,
		    fdm2d    fdm)
/*< apply boundary sponge >*/
{
    int iz,ix,ib,ibz,ibx;
    float w;

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,ix,ibz,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibz = fdm->nzpad-ib-1;
	for(ix=0; ix<fdm->nxpad; ix++) {
	    uu[ix][ib ] *= w; /*    top sponge */
	    uu[ix][ibz] *= w; /* bottom sponge */
	}
    }

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ibx,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibx = fdm->nxpad-ib-1;
	for(iz=0; iz<fdm->nzpad; iz++) {
	    uu[ib ][iz] *= w; /*   left sponge */
	    uu[ibx][iz] *= w; /*  right sponge */
	}
    }
}

/*------------------------------------------------------------*/
void sponge3d_apply(float  ***uu,
		    sponge   spo,
		    fdm3d    fdm)
/*< apply boundary sponge >*/
{
    int iz,ix,iy,ib,ibz,ibx,iby;
    float w;

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,ix,iy,ibz,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibz = fdm->nzpad-ib-1;
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(ix=0; ix<fdm->nxpad; ix++) {
		uu[iy][ix][ib ] *= w; /* z min */
		uu[iy][ix][ibz] *= w; /* z max */
	    }
	}
    }

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,iy,ibx,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibx = fdm->nxpad-ib-1;
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		uu[iy][ib ][iz] *= w; /* x min */
		uu[iy][ibx][iz] *= w; /* x max */
	    }
	}
    }

    if (fdm->nypad>1) {
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ix,iby,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];
	
	iby = fdm->nypad-ib-1;
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		uu[ib ][ix][iz] *= w; /* y min */
		uu[iby][ix][iz] *= w; /* y max */
	    }
	}

    }
    }
}

void sponge3d_apply_complex(sf_complex  ***uu,
		            sponge   spo,
		            fdm3d    fdm)
/*< apply boundary sponge for complex-valued wavefield >*/
{
    int iz,ix,iy,ib,ibz,ibx,iby;
    float w;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ib,ix,iy,ibz,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibz = fdm->nzpad-ib-1;
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(ix=0; ix<fdm->nxpad; ix++) {
#ifdef SF_HAS_COMPLEX_H
		uu[iy][ix][ib ] *= w; /* z min */
		uu[iy][ix][ibz] *= w; /* z max */
#else
		uu[iy][ix][ib ] = sf_crmul(uu[iy][ix][ib ],w); /* z min */
		uu[iy][ix][ibz] = sf_crmul(uu[iy][ix][ibz],w); /* z max */
#endif
	    }
	}
    }

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ib,iz,iy,ibx,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibx = fdm->nxpad-ib-1;
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
#ifdef SF_HAS_COMPLEX_H
		uu[iy][ib ][iz] *= w; /* x min */
		uu[iy][ibx][iz] *= w; /* x max */
#else
		uu[iy][ib ][iz] = sf_crmul(uu[iy][ib ][iz],w); /* x min */
		uu[iy][ibx][iz] = sf_crmul(uu[iy][ibx][iz],w); /* x max */
#endif
	    }
	}
    }

    if (fdm->nypad>1) {
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ib,iz,ix,iby,w)			\
    shared(fdm,spo,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];
	
	iby = fdm->nypad-ib-1;
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
#ifdef SF_HAS_COMPLEX_H
		uu[ib ][ix][iz] *= w; /* y min */
		uu[iby][ix][iz] *= w; /* y max */
#else
		uu[ib ][ix][iz] = sf_crmul(uu[ib ][ix][iz],w); /* y min */
		uu[iby][ix][iz] = sf_crmul(uu[iby][ix][iz],w); /* y max */
#endif
	    }
	}

    }
    }
}


bool cfl_generic(
    float vpmin, float vpmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals, char *wave)
/*< cfl check for both 2d and 3d acoustic fdcode >*/
{
    int dim = 3;
    float dmin;
    float tdp;
    float wplength;
    bool passed;

    if (dy < 0)dim = 2;
    
    if (dim == 2) dmin = SF_MIN(dx,dz);
    else dmin = SF_MIN(dx,SF_MIN(dy,dz));
    
    tdp = dt * vpmax *sqrt(2); /* maximum distance */
    if (dmin > tdp) sf_warning("CFL: Stability check ... %s-wave... PASSED", wave);
    else {
        sf_error("CFL: Stability check ... FAILED ... minimum grid sampling: %f !> %f", dmin, tdp);
    }
    wplength = safety*vpmin / fmax;
    
    if (dim == 2) passed = wplength > intervals*sqrt(dx*dx+dz*dz);
    else passed = wplength > intervals*sqrt(dx*dx+dy*dy+dz*dz);
    
    if (passed) {
        sf_warning("CFL: Accuracy check ... %s-wave ... PASSED", wave);
    }
    else {
        if (dim == 2) {
            sf_warning("CFL: fmax must be less than < %f",
                       safety*vpmin/intervals/sqrt(dx*dx+dz*dz));
        }
        else {
            sf_warning("CFL: fmax must be less than < %f",
                       safety*vpmin/intervals/sqrt(dx*dx+dy*dy+dz*dz));
            sf_error("CFL: Accuracy check ... %s-wave ... FAILED",wave);
        }
    }
    return true;
}

bool cfl_elastic(
    float vpmin, float vpmax,
    float vsmin, float vsmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals)
/*< cfl check for both 2d and 3d elastic fdcode >*/
{
    bool ppass = cfl_generic(vpmin,vpmax,
                             dx,dy,dz,dt,
                             fmax,safety,intervals,"P");
    bool spass = cfl_generic(vsmin,vsmax,
                             dx,dy,dz,dt,
                             fmax,safety,intervals,"S");
    
    return ppass && spass;
}

bool cfl_acoustic(
    float vpmin, float vpmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals)
/*< cfl check for acoustic wave equation >*/
{
    bool ppass = cfl_generic(vpmin,vpmax,
                             dx,dy,dz,dt,
                             fmax,safety,intervals,"P");
    return ppass;
}

/*------------------------------------------------------------*/
dft3d dft3d_init(int pad1,
                 bool rtoc,
                 bool padio,
                 fdm3d fdm)
/*< init 3d fft >*/
{
    dft3d dft;

    dft = (dft3d) sf_alloc(1,sizeof(*dft));
    fft_init(pad1,fdm->nzpad,fdm->nxpad,fdm->nypad,&dft->nkz,&dft->nkx,&dft->nky,rtoc,padio);

    dft->dkz = 1./(dft->nkz*fdm->dz); dft->okz = -0.5/fdm->dz;
    dft->dkx = 1./(dft->nkx*fdm->dx); dft->okx = -0.5/fdm->dx;
    if(dft->nky>1) {
        dft->dky = 1./(dft->nky*fdm->dy); dft->oky = -0.5/fdm->dy;
    } else {
        dft->dky = 0.; dft->oky = 0.;
    }
    if(fdm->verb) sf_warning("DFT: nkz=%d, dkz=%g, okz=%g, nkx=%d, dkx=%g, okx=%g, nky=%d, dky=%g, oky=%g",dft->nkz,dft->dkz,dft->okz,dft->nkx,dft->dkx,dft->okx,dft->nky,dft->dky,dft->oky);

    return dft;
}

/*------------------------------------------------------------*/
void dft3d_finalize()
/*< finalize 3d fft >*/
{
    fft_finalize();
}

/*------------------------------------------------------------*/
lps3d lps3d_init(float pcut /* pcut/2 is the width of tapered region w.r.t. 1 */,
                 float fcut /* cutoff frequency */,
                 float vmax /* maximum p-wave velocity */,
                 dft3d dft)
/*< initialize low-pass filter coefficients >*/
{
    lps3d lps;
    int ikz,ikx,iky;
    float kz,kx,ky,kk,kmax,kbnd;

    lps = (lps3d) sf_alloc(1,sizeof(*lps));

    lps->flt = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);

    kmax = sqrtf(powf(2.*SF_PI*dft->okz,2.) + powf(2.*SF_PI*dft->okx,2.) + powf(2.*SF_PI*dft->oky,2.));
    kbnd = 2.*SF_PI*fcut/vmax;

    sf_warning("Tukey taper: pcut = %f, fcut = %f, vmax = %f.",pcut,fcut,vmax);
    if (kbnd > kmax) {
        sf_warning("cutoff wavenumber %f larger than maximum wavenumber %f! Setting kbnd = kmax...",kbnd,kmax);
	kbnd = kmax;
    }

    /* construct the pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,kz,kx,ky,kk)		\
    shared(dft,lps,kmax,kbnd,pcut)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        ky = 2.*SF_PI*(dft->oky + iky*dft->dky);
        for (ikx=0; ikx<dft->nkx; ikx++) {
            kx = 2.*SF_PI*(dft->okx+ikx*dft->dkx);
            for (ikz=0; ikz<dft->nkz; ikz++) {
                kz = 2.*SF_PI*(dft->okz+ikz*dft->dkz);
                kk = sqrtf(kz*kz + kx*kx + ky*ky);

                if (pcut==0) {
                    lps->flt[iky][ikx][ikz] = 1.;
                } else {
                    if (kk > kbnd)
                        lps->flt[iky][ikx][ikz] = 0.; 
                    else if (kk >= kbnd*(1.-0.5*pcut))
                        lps->flt[iky][ikx][ikz] = 0.5*(1.+cosf(SF_PI*(2.*kk/(pcut*kbnd)-2./pcut+1.)));
                    else
                        lps->flt[iky][ikx][ikz] = 1.;
                }
            }
        }
    }

    return lps;
}

/*------------------------------------------------------------*/
static sf_complex *wavek =NULL;
static sf_complex *wavekz=NULL;
static sf_complex *wavekx=NULL;
static sf_complex *waveky=NULL;

ksp3d ksp3d_make(dft3d dft)
/*< init k-space arrays for pseudo-spectral method >*/
{
    int ikz,ikx,iky;
    int nk;
    float kz,kx,ky;
    ksp3d ksp;

    ksp = (ksp3d) sf_alloc(1,sizeof(*ksp));
    nk = dft->nky*dft->nkx*dft->nkz;
    wavek  = sf_complexalloc(nk);
    wavekz = sf_complexalloc(nk);
    wavekx = sf_complexalloc(nk);
    waveky = sf_complexalloc(nk);
    ksp->kz = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);
    ksp->kx = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);
    ksp->ky = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);

    /* construct the pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,kz,kx,ky)		\
    shared(dft,ksp)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        ky = 2.*SF_PI*(dft->oky + iky*dft->dky);
        for (ikx=0; ikx<dft->nkx; ikx++) {
            kx = 2.*SF_PI*(dft->okx+ikx*dft->dkx);
            for (ikz=0; ikz<dft->nkz; ikz++) {
                kz = 2.*SF_PI*(dft->okz+ikz*dft->dkz);
                ksp->kz[iky][ikx][ikz] = kz; 
                ksp->kx[iky][ikx][ikz] = kx; 
                ksp->ky[iky][ikx][ikz] = ky; 
            }
        }
    }

    return ksp;
}

/*------------------------------------------------------------*/
void ksp3d_apply(float *wavedx,
                 float *wavedy,
                 float *wavedz,
                 float *wave,
                 dft3d dft,
                 ksp3d ksp)
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/
{
    int ikz,ikx,iky,ik;

    fft(wave,wavek);

    /* apply pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,ksp,wavek,wavekz,wavekx,waveky)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        for (ikx=0; ikx<dft->nkx; ikx++) {
            for (ikz=0; ikz<dft->nkz; ikz++) {
                ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                wavekz[ik] = wavek[ik]*sf_cmplx(0.,ksp->kz[iky][ikx][ikz]); 
                wavekx[ik] = wavek[ik]*sf_cmplx(0.,ksp->kx[iky][ikx][ikz]); 
                waveky[ik] = wavek[ik]*sf_cmplx(0.,ksp->ky[iky][ikx][ikz]); 
#else
                wavekz[ik] = sf_cmul(wavek[ik],sf_cmplx(0.,ksp->kz[iky][ikx][ikz])); 
                wavekx[ik] = sf_cmul(wavek[ik],sf_cmplx(0.,ksp->kx[iky][ikx][ikz])); 
                waveky[ik] = sf_cmul(wavek[ik],sf_cmplx(0.,ksp->ky[iky][ikx][ikz])); 
#endif
            }
        }
    }

    ifft(wavedz,wavekz);
    ifft(wavedx,wavekx);
    ifft(wavedy,waveky);
}

/*------------------------------------------------------------*/
void ksp3d_apply1(float *wavedx,
                  float *wavedy,
                  float *wavedz,
                  float *wavex,
                  float *wavey,
                  float *wavez,
                  dft3d dft,
                  ksp3d ksp)
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/
{
    int ikz,ikx,iky,ik;

    fft(wavez,wavekz);
    fft(wavex,wavekx);
    fft(wavey,waveky);

    /* apply pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,ksp,wavekz,wavekx,waveky)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        for (ikx=0; ikx<dft->nkx; ikx++) {
            for (ikz=0; ikz<dft->nkz; ikz++) {
                ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                wavekz[ik] = wavekz[ik]*sf_cmplx(0.,ksp->kz[iky][ikx][ikz]); 
                wavekx[ik] = wavekx[ik]*sf_cmplx(0.,ksp->kx[iky][ikx][ikz]); 
                waveky[ik] = waveky[ik]*sf_cmplx(0.,ksp->ky[iky][ikx][ikz]); 
#else
                wavekz[ik] = sf_cmul(wavekz[ik],sf_cmplx(0.,ksp->kz[iky][ikx][ikz])); 
                wavekx[ik] = sf_cmul(wavekx[ik],sf_cmplx(0.,ksp->kx[iky][ikx][ikz])); 
                waveky[ik] = sf_cmul(waveky[ik],sf_cmplx(0.,ksp->ky[iky][ikx][ikz])); 
#endif
            }
        }
    }

    ifft(wavedz,wavekz);
    ifft(wavedx,wavekx);
    ifft(wavedy,waveky);
}

/*------------------------------------------------------------*/
void ksp3d_apply2(float *wavedxx,
                  float *wavedyy,
                  float *wavedzz,
                  float *wavedxy,
                  float *wavedyz,
                  float *wavedzx,
                  float *wave,
                  dft3d dft,
                  ksp3d ksp)
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/
{
    int ikz,ikx,iky,ik;

    fft(wave,wavek);

    /* apply pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,ksp,wavek,wavekz,wavekx,waveky)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        for (ikx=0; ikx<dft->nkx; ikx++) {
            for (ikz=0; ikz<dft->nkz; ikz++) {
                ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                wavekz[ik] = -wavek[ik]*ksp->kz[iky][ikx][ikz]*ksp->kz[iky][ikx][ikz]; 
                wavekx[ik] = -wavek[ik]*ksp->kx[iky][ikx][ikz]*ksp->kx[iky][ikx][ikz]; 
                waveky[ik] = -wavek[ik]*ksp->ky[iky][ikx][ikz]*ksp->ky[iky][ikx][ikz]; 
#else
                wavekz[ik] = sf_crmul(wavek[ik],-ksp->kz[iky][ikx][ikz]*ksp->kz[iky][ikx][ikz]); 
                wavekx[ik] = sf_crmul(wavek[ik],-ksp->kx[iky][ikx][ikz]*ksp->kx[iky][ikx][ikz]); 
                waveky[ik] = sf_crmul(wavek[ik],-ksp->ky[iky][ikx][ikz]*ksp->ky[iky][ikx][ikz]); 
#endif
            }
        }
    }

    ifft(wavedzz,wavekz);
    ifft(wavedxx,wavekx);
    ifft(wavedyy,waveky);

    /* apply pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,ksp,wavek,wavekz,wavekx,waveky)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        for (ikx=0; ikx<dft->nkx; ikx++) {
            for (ikz=0; ikz<dft->nkz; ikz++) {
                ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                wavekz[ik] = -wavek[ik]*ksp->kx[iky][ikx][ikz]*ksp->ky[iky][ikx][ikz]; 
                wavekx[ik] = -wavek[ik]*ksp->ky[iky][ikx][ikz]*ksp->kz[iky][ikx][ikz]; 
                waveky[ik] = -wavek[ik]*ksp->kz[iky][ikx][ikz]*ksp->kx[iky][ikx][ikz]; 
#else
                wavekz[ik] = sf_crmul(wavek[ik],-ksp->kx[iky][ikx][ikz]*ksp->ky[iky][ikx][ikz]); 
                wavekx[ik] = sf_crmul(wavek[ik],-ksp->ky[iky][ikx][ikz]*ksp->kz[iky][ikx][ikz]); 
                waveky[ik] = sf_crmul(wavek[ik],-ksp->kz[iky][ikx][ikz]*ksp->kx[iky][ikx][ikz]); 
#endif
            }
        }
    }

    ifft(wavedxy,wavekz);
    ifft(wavedyz,wavekx);
    ifft(wavedzx,waveky);

}

/*------------------------------------------------------------*/
void ksp3d_finalize()
/*< free static memory allocated for ksp >*/
{
    free(wavek);
    free(wavekz);
    free(wavekx);
    free(waveky);
}

/*------------------------------------------------------------*/
static sf_complex *waveka =NULL;
static sf_complex *wavekb =NULL;

vksp3d vksp3d_make(float gpavg,
                   float gsavg,
                   dft3d dft,
                   lps3d lps)
/*< init k-space arrays for fractional-Laplacian pseudo-spectral method >*/
{
    int ikz,ikx,iky;
    int nk;
    float kz,kx,ky,k2,k2min;
    vksp3d vksp;

    vksp = (vksp3d) sf_alloc(1,sizeof(*vksp));
    nk = dft->nky*dft->nkx*dft->nkz;
    k2min = dft->dkz*dft->dkz + dft->dkx*dft->dkx + dft->dky*dft->dky; /* protection against zero */

    waveka  = sf_complexalloc(nk);
    wavekb  = sf_complexalloc(nk);
    vksp->kpa = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);
    vksp->ksa = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);
    vksp->kpb = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);
    vksp->ksb = sf_floatalloc3(dft->nkz,dft->nkx,dft->nky);

    /* construct the pseudo-spectral op */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,kz,kx,ky,k2)		\
    shared(dft,lps,vksp,gpavg,gsavg)
#endif
    for (iky=0; iky<dft->nky; iky++) {
        ky = 2.*SF_PI*(dft->oky + iky*dft->dky);
        for (ikx=0; ikx<dft->nkx; ikx++) {
            kx = 2.*SF_PI*(dft->okx+ikx*dft->dkx);
            for (ikz=0; ikz<dft->nkz; ikz++) {
                kz = 2.*SF_PI*(dft->okz+ikz*dft->dkz);
                k2 = kz*kz + kx*kx + ky*ky;

                vksp->kpb[iky][ikx][ikz] = powf(k2,gpavg);
                vksp->ksb[iky][ikx][ikz] = powf(k2,gsavg);

                if (k2<k2min) k2=k2min;
                vksp->kpa[iky][ikx][ikz] = powf(k2,gpavg-0.5);
                vksp->ksa[iky][ikx][ikz] = powf(k2,gsavg-0.5);

                if (NULL!=lps) { 
                    vksp->kpa[iky][ikx][ikz] *= lps->flt[iky][ikx][ikz];
                    vksp->ksa[iky][ikx][ikz] *= lps->flt[iky][ikx][ikz];
                }
            }
        }
    }

    return vksp;
}

/*------------------------------------------------------------*/
void vksp3d_apply(float *wavea,
                  float *waveb,
                  int mode /* 1 -> p mode; 2 -> s mode */,
                  dft3d dft,
                  vksp3d vksp)
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/
{
    int ikz,ikx,iky,ik;

    fft(wavea,waveka);
    fft(waveb,wavekb);

    /* apply pseudo-spectral op */
    if (mode == 1) { /* p mode */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,vksp,waveka,wavekb)
#endif
        for (iky=0; iky<dft->nky; iky++) {
            for (ikx=0; ikx<dft->nkx; ikx++) {
                for (ikz=0; ikz<dft->nkz; ikz++) {
                    ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                    waveka[ik] = waveka[ik]*vksp->kpa[iky][ikx][ikz]; 
                    wavekb[ik] = wavekb[ik]*vksp->kpb[iky][ikx][ikz]; 
#else
                    waveka[ik] = sf_crmul(waveka[ik],vksp->kpa[iky][ikx][ikz]); 
                    wavekb[ik] = sf_crmul(wavekb[ik],vksp->kpb[iky][ikx][ikz]); 
#endif
                }
            }
        }
    } else { /* s mode */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ikz,ikx,iky,ik)                     \
    shared(dft,vksp,waveka,wavekb)
#endif
        for (iky=0; iky<dft->nky; iky++) {
            for (ikx=0; ikx<dft->nkx; ikx++) {
                for (ikz=0; ikz<dft->nkz; ikz++) {
                    ik = (iky*dft->nkx+ikx)*dft->nkz+ikz; 
#ifdef SF_HAS_COMPLEX_H
                    waveka[ik] = waveka[ik]*vksp->ksa[iky][ikx][ikz]; 
                    wavekb[ik] = wavekb[ik]*vksp->ksb[iky][ikx][ikz]; 
#else
                    waveka[ik] = sf_crmul(waveka[ik],vksp->ksa[iky][ikx][ikz]); 
                    wavekb[ik] = sf_crmul(wavekb[ik],vksp->ksb[iky][ikx][ikz]); 
#endif
                }
            }
        }
    }

    ifft(wavea,waveka);
    ifft(waveb,wavekb);
}

/*------------------------------------------------------------*/
void vksp3d_finalize()
/*< free static memory allocated for fft and ksp >*/
{
    free(waveka);
    free(wavekb);
}

/*------------------------------------------------------------*/
static sf_complex *cwave   =NULL;
static sf_complex **cwaven2=NULL;
static sf_complex **waven2 =NULL;
static sf_complex ***waves =NULL;

clr3d clr3d_make(int *n2s,
                 fdm3d fdm)
/*< prepare lowrank arrays for rite method >*/
{
    int i,n2_sum,n2_max;
    clr3d clr;

    clr = (clr3d) sf_alloc(1,sizeof(*clr));

    clr->n2_start = sf_intalloc(6);
    clr->n2_local = sf_intalloc(6);
    clr->map = sf_intalloc2(3,3);

    n2_sum = 0;
    n2_max = 0;
    for (i=0; i<6; i++) {
        clr->n2_start[i] = n2_sum; 
        clr->n2_local[i] = n2s[i]; 
        n2_sum += n2s[i];
        if(n2_max<n2s[i]) n2_max = n2s[i];
    }
    clr->n2_sum=n2_sum;
    clr->n2_max=n2_max;

    if(fdm->verb) {
        sf_warning("n2_start=%d,%d,%d,%d,%d,%d",clr->n2_start[0],clr->n2_start[1],clr->n2_start[2],clr->n2_start[3],clr->n2_start[4],clr->n2_start[5]);
        sf_warning("n2_local=%d,%d,%d,%d,%d,%d",clr->n2_local[0],clr->n2_local[1],clr->n2_local[2],clr->n2_local[3],clr->n2_local[4],clr->n2_local[5]);
        sf_warning("n2_sum=%d, n2_max=%d",clr->n2_sum,clr->n2_max);
    }

    clr->map[0][0] = 0; clr->map[0][1] = 3; clr->map[0][2] = 4;
    clr->map[1][0] = 3; clr->map[1][1] = 1; clr->map[1][2] = 5;
    clr->map[2][0] = 4; clr->map[2][1] = 5; clr->map[2][2] = 2;

    return clr;
}

clr3d clr3d_make2(int *n2s,
                  fdm3d fdm)
/*< prepare lowrank arrays for rite method >*/
{
    int i,n2_sum,n2_max;
    clr3d clr;

    clr = (clr3d) sf_alloc(1,sizeof(*clr));

    clr->n2_start = sf_intalloc(9);
    clr->n2_local = sf_intalloc(9);
    clr->map = sf_intalloc2(3,3);

    n2_sum = 0;
    n2_max = 0;
    for (i=0; i<9; i++) {
        clr->n2_start[i] = n2_sum; 
        clr->n2_local[i] = n2s[i]; 
        n2_sum += n2s[i];
        if(n2_max<n2s[i]) n2_max = n2s[i];
    }
    clr->n2_sum=n2_sum;
    clr->n2_max=n2_max;

    if(fdm->verb) {
        sf_warning("n2_start=%d,%d,%d,%d,%d,%d,%d,%d,%d",clr->n2_start[0],clr->n2_start[1],clr->n2_start[2],clr->n2_start[3],clr->n2_start[4],clr->n2_start[5],clr->n2_start[6],clr->n2_start[7],clr->n2_start[8]);
        sf_warning("n2_local=%d,%d,%d,%d,%d,%d,%d,%d,%d",clr->n2_local[0],clr->n2_local[1],clr->n2_local[2],clr->n2_local[3],clr->n2_local[4],clr->n2_local[5],clr->n2_local[6],clr->n2_local[7],clr->n2_local[8]);
        sf_warning("n2_sum=%d, n2_max=%d",clr->n2_sum,clr->n2_max);
    }

    clr->map[0][0] = 0; clr->map[0][1] = 1; clr->map[0][2] = 2;
    clr->map[1][0] = 3; clr->map[1][1] = 4; clr->map[1][2] = 5;
    clr->map[2][0] = 6; clr->map[2][1] = 7; clr->map[2][2] = 8;

    return clr;
}

void clr3d_init(fdm3d fdm,
                dft3d dft,
                clr3d clr)
/*< initialize computation array >*/
{
    int nxyz,nk,ik,i,j;

    nxyz = fdm->nypad *fdm->nxpad *fdm->nzpad ;
    nk   = dft->nky*dft->nkx*dft->nkz;
    cwave  = sf_complexalloc(nk);
    cwaven2= sf_complexalloc2(nk,  clr->n2_max);
    waven2 = sf_complexalloc2(nxyz,clr->n2_max);
    waves  = sf_complexalloc3(nxyz,3,3);

    for (i=0; i<3; i++)     {
        for (j=0; j<3; j++) {
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ik)                                 \
    shared(nxyz,i,j,waves)
#endif
            for (ik=0; ik<nxyz; ik++)
                waves[i][j][ik] = sf_cmplx(0.,0.);
        }
    }
}

/*------------------------------------------------------------*/
void clr3d_apply(sf_complex **uo,
                 sf_complex **ui,
                 sf_complex **lt,
                 sf_complex **rt,
                 fdm3d fdm,
                 dft3d dft,
                 clr3d clr)
/*< apply lowrank matrices for time stepping (can be in-place) >*/
{
    int nxyz,nk,ik,ir,ic,im,in,iz,ix,iy,i;
    sf_complex c;

    nxyz = fdm->nypad *fdm->nxpad *fdm->nzpad ;
    nk   = dft->nky*dft->nkx*dft->nkz;

    for (ic=0; ic<3; ic++) {

        fft(ui[ic],cwave);

        for (ir=0; ir<3; ir++) {

#ifdef _OPENMP
#pragma omp parallel for              \
    private(im,in,ik)                 \
    shared(clr,cwaven2,cwave,rt,nk)
#endif
            for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++) {
                in = clr->n2_start[clr->map[ir][ic]]+im;
                for (ik=0; ik<nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                    cwaven2[im][ik] = cwave[ik]*rt[in][ik];
#else
                    cwaven2[im][ik] = sf_cmul(cwave[ik],rt[in][ik]);
#endif
                }
            }
            for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++)
                ifft(waven2[im],cwaven2[im]);

#ifdef _OPENMP
#pragma omp parallel for			\
    private(iy,ix,iz,i,c,im,in)                 \
    shared(fdm,clr,lt,waven2,waves,ir,ic)
#endif
            for (iy=0; iy<fdm->nypad; iy++)         {
                for (ix=0; ix<fdm->nxpad; ix++)     {
                    for (iz=0; iz<fdm->nzpad; iz++) {
                        i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz; /* flattened coordinate */
                        c = sf_cmplx(0.,0.);
                        for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++) {
                            in = clr->n2_start[clr->map[ir][ic]]+im;
#ifdef SF_HAS_COMPLEX_H
                            c += lt[in][i]*waven2[im][i];
#else
                            c = sf_cadd(c,sf_cmul(lt[in][i],waven2[im][i]));
#endif
                        }
                        waves[ir][ic][i] = c;
                    }
                }
            }

        } /* ir loop */

    } /* ic loop */

    /* linear combination forms the output vector */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ik)                                 \
    shared(nxyz,waves,uo)
#endif
    for (ik=0; ik<nxyz; ik++) {
#ifdef SF_HAS_COMPLEX_H
        uo[0][ik] = waves[0][0][ik] + waves[0][1][ik] + waves[0][2][ik];
        uo[1][ik] = waves[1][0][ik] + waves[1][1][ik] + waves[1][2][ik];
        uo[2][ik] = waves[2][0][ik] + waves[2][1][ik] + waves[2][2][ik];
#else
        uo[0][ik] = sf_cadd(waves[0][0][ik],sf_cadd(waves[0][1][ik],waves[0][2][ik]));
        uo[1][ik] = sf_cadd(waves[1][0][ik],sf_cadd(waves[1][1][ik],waves[1][2][ik]));
        uo[2][ik] = sf_cadd(waves[2][0][ik],sf_cadd(waves[2][1][ik],waves[2][2][ik]));
#endif
    }

}

/*------------------------------------------------------------*/
void clr3d_apply2(sf_complex **u2,
                  sf_complex **u1,
                  sf_complex **u0,
                  sf_complex **lt,
                  sf_complex **rt,
                  fdm3d fdm,
                  dft3d dft,
                  clr3d clr)
/*< apply lowrank matrices for time stepping (can be in-place) -- two-step version >*/
{
    int nxyz,nk,ik,ir,ic,im,in,iz,ix,iy,i;
    sf_complex c;

    nxyz = fdm->nypad *fdm->nxpad *fdm->nzpad ;
    nk   = dft->nky*dft->nkx*dft->nkz;

    for (ic=0; ic<3; ic++) {

        fft(u1[ic],cwave);

        for (ir=0; ir<3; ir++) {

#ifdef _OPENMP
#pragma omp parallel for              \
    private(im,in,ik)                 \
    shared(clr,cwaven2,cwave,rt,nk)
#endif
            for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++) {
                in = clr->n2_start[clr->map[ir][ic]]+im;
                for (ik=0; ik<nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                    cwaven2[im][ik] = cwave[ik]*rt[in][ik];
#else
                    cwaven2[im][ik] = sf_cmul(cwave[ik],rt[in][ik]);
#endif
                }
            }
            for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++)
                ifft(waven2[im],cwaven2[im]);

#ifdef _OPENMP
#pragma omp parallel for			\
    private(iy,ix,iz,i,c,im,in)                 \
    shared(fdm,clr,lt,waven2,waves,ir,ic)
#endif
            for (iy=0; iy<fdm->nypad; iy++)         {
                for (ix=0; ix<fdm->nxpad; ix++)     {
                    for (iz=0; iz<fdm->nzpad; iz++) {
                        i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz; /* flattened coordinate */
                        c = sf_cmplx(0.,0.);
                        for (im=0; im<clr->n2_local[clr->map[ir][ic]]; im++) {
                            in = clr->n2_start[clr->map[ir][ic]]+im;
#ifdef SF_HAS_COMPLEX_H
                            c += lt[in][i]*waven2[im][i];
#else
                            c = sf_cadd(c,sf_cmul(lt[in][i],waven2[im][i]));
#endif
                        }
                        waves[ir][ic][i] = c;
                    }
                }
            }

        } /* ir loop */

    } /* ic loop */

    /* linear combination forms the output vector */
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ik)                                 \
    shared(nxyz,waves,u2,u1,u0)
#endif
    for (ik=0; ik<nxyz; ik++) {
#ifdef SF_HAS_COMPLEX_H
        u2[0][ik] = 2.*u1[0][ik] - u0[0][ik] + waves[0][0][ik] + waves[0][1][ik] + waves[0][2][ik];
        u2[1][ik] = 2.*u1[1][ik] - u0[1][ik] + waves[1][0][ik] + waves[1][1][ik] + waves[1][2][ik];
        u2[2][ik] = 2.*u1[2][ik] - u0[2][ik] + waves[2][0][ik] + waves[2][1][ik] + waves[2][2][ik];
#else
        u2[0][ik] = sf_cadd(sf_csub(sf_crmul(u1[0][ik],2.),u0[0][ik]),sf_cadd(waves[0][0][ik],sf_cadd(waves[0][1][ik],waves[0][2][ik])));
        u2[1][ik] = sf_cadd(sf_csub(sf_crmul(u1[1][ik],2.),u0[1][ik]),sf_cadd(waves[1][0][ik],sf_cadd(waves[1][1][ik],waves[1][2][ik])));
        u2[2][ik] = sf_cadd(sf_csub(sf_crmul(u1[2][ik],2.),u0[2][ik]),sf_cadd(waves[2][0][ik],sf_cadd(waves[2][1][ik],waves[2][2][ik])));
#endif
    }

}

/*------------------------------------------------------------*/
void clr3d_apply_dbg(sf_complex **uo,
                     sf_complex **ui,
                     sf_complex **lt,
                     sf_complex **rt,
                     fdm3d fdm,
                     dft3d dft,
                     clr3d clr)
/*< apply lowrank matrices for time stepping (can be in-place) >*/
{
    int nxyz,nk,ik,ir,ic,im,in,iz,ix,iy,i;
    sf_complex c;

    nxyz = fdm->nypad *fdm->nxpad *fdm->nzpad ;
    nk   = dft->nky*dft->nkx*dft->nkz;

    for (ic=0; ic<1; ic++) {

        fft(ui[ic],cwave);

        for (ir=0; ir<1; ir++) {

            for (im=0; im<clr->n2_local[clr->map[0][0]]; im++) {
                in = clr->n2_start[clr->map[0][0]]+im;
                for (ik=0; ik<nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                    cwaven2[im][ik] = cwave[ik]*rt[in][ik];
#else
                    cwaven2[im][ik] = sf_cmul(cwave[ik],rt[in][ik]);
#endif
                }
            }
            for (im=0; im<clr->n2_local[clr->map[0][0]]; im++)
                ifft(waven2[im],cwaven2[im]);

            for (iy=0; iy<fdm->nypad; iy++)         {
                for (ix=0; ix<fdm->nxpad; ix++)     {
                    for (iz=0; iz<fdm->nzpad; iz++) {
                        i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz; /* flattened coordinate */
                        c = sf_cmplx(0.,0.);
                        for (im=0; im<clr->n2_local[clr->map[0][0]]; im++) {
                            in = clr->n2_start[clr->map[0][0]]+im;
#ifdef SF_HAS_COMPLEX_H
                            c += lt[in][i]*waven2[im][i];
#else
                            c = sf_cadd(c,sf_cmul(lt[in][i],waven2[im][i]));
#endif
                        }
                        waves[ir][ic][i] = c;
                    }
                }
            }

        } /* ir loop */

    } /* ic loop */

    /* linear combination forms the output vector */
    for (ic=0; ic<1; ic++) {
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ik)                                 \
    shared(nk,ir,waves,uo)
#endif
        for (ik=0; ik<nxyz; ik++) {
#ifdef SF_HAS_COMPLEX_H
            uo[ic][ik] = waves[0][ic][ik];
#else
            uo[ic][ik] = waves[0][ic][ik];
#endif
        }
    }

}

/*------------------------------------------------------------*/
void clr3d_finalize()
/*< free static memory allocated for clr >*/
{
    free(cwave);
    free(*cwaven2); free(cwaven2);
    free(*waven2); free(waven2);
    free(**waves); free(*waves); free(waves);
}

/*------------------------------------------------------------*/
double walltime( double *t0 )
/*< report runtime >*/
{
    double mic, time;
    double mega = 0.000001;
    struct timeval tp;
    struct timezone tzp;
    static long base_sec = 0;
    static long base_usec = 0;

    (void) gettimeofday(&tp,&tzp);
    if (base_sec == 0)
    {
        base_sec = tp.tv_sec;
        base_usec = tp.tv_usec;
    }

    time = (double) (tp.tv_sec - base_sec);
    mic = (double) (tp.tv_usec - base_usec);
    time = (time + mic * mega) - *t0;
    return(time);
}

/*------------------------------------------------------------*/
dat3d dat3d_init(int ns_,
                 int nr_,
                 int nc_,
		 sf_axis at_,
		 sf_axis ax_,
                 sf_axis ay_)
/*< init 3D data >*/
{
    dat3d dat;
    dat = (dat3d) sf_alloc(1,sizeof(*dat));

    dat->ns = ns_;
    dat->nr = nr_;
    dat->nc = nc_;

    dat->nt = sf_n(at_);
    dat->nx = sf_n(ax_);
    dat->ny = sf_n(ay_);
    
    dat->dt = sf_d(at_);
    dat->dx = sf_d(ax_);
    dat->dy = sf_d(ay_);

    dat->ot = sf_o(at_);
    dat->ox = sf_o(ax_);
    dat->oy = sf_o(ay_);

    return dat;
}

/*------------------------------------------------------------*/
mut3d mut3d_make(float t0  /* source delay */,
                 float velw /* water velocity */,
                 float eps /* decay parameter */,
                 pt3d* ss,
                 pt3d* rr,
                 dat3d dat)
/*< init 3D mutting, currently only handles single source >*/
{
    mut3d mut;
    int   it,ir;
    float dist,t1,tt;
    
    mut = (mut3d) sf_alloc(1,sizeof(*mut));
    mut->wt = sf_floatalloc2(dat->nt,dat->nr);

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ir,it,dist,t1,tt)                   \
    shared(t0,velw,eps,ss,rr,dat,mut)
#endif
    for (ir=0;ir<dat->nr;ir++) {
	
        dist = powf(rr[ir].z-ss[0].z,2.) + powf(rr[ir].x-ss[0].x,2.) + powf(rr[ir].y-ss[0].y,2.);
        t1 = t0 + sqrtf(dist)/velw;
        for (it=0; it<dat->nt; it++) {
            tt = it*dat->dt;
            if (tt > t1)
                mut->wt[ir][it] = 1.;
            else
                mut->wt[ir][it] = expf(eps*(tt-t1));
        }

    }

    return mut;
}

/*------------------------------------------------------------*/
void mut3d_apply(float ***dd,
                 dat3d dat,
                 mut3d mut)
/*< apply in-place muting >*/
{
    int ir,ic,it;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ir,ic,it)                           \
    shared(dd,mut,dat)
#endif
    for         (it=0;it<dat->nt;it++) {
        for     (ic=0;ic<dat->nc;ic++) {
            for (ir=0;ir<dat->nr;ir++) {
                dd[it][ic][ir] *= mut->wt[ir][it];
            }
        }
    }

}

