#include <rsf.h>
#include "fdlib.h"
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif 


#ifndef _fdlib_h

typedef struct fdm2 *fdm2d;
/*^*/

typedef struct fdm3 *fdm3d;
/*^*/

typedef struct bcoef2d *bell2d;
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

struct bcoef2d{
	int nbell;
	float **bell;
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



#endif

static float** bell;
static float***bell3d;
static int    nbell;

/*------------------------------------------------------------*/
/* stencils and stencil parameters*/

#define NOP 4 /* derivative operator half-size */
/* Muir's derivative operator */
#define C1 +0.598144  //  1225/ 1024      /2 
#define C2 -0.039876  // -1225/(1024*  15)/2
#define C3 +0.004785  //  1225/(1024* 125)/2 
#define C4 -0.000348  // -1225/(1024*1715)/2 

/***************/
/* 2D STENCILS */
/***************/
/*  2D forward FD derivative stencils */
#define Fz(a,ix,iz,s) (2*C4*(a[ix  ][iz+4] - a[ix  ][iz-3]) +     \
                       2*C3*(a[ix  ][iz+3] - a[ix  ][iz-2]) +     \
                       2*C2*(a[ix  ][iz+2] - a[ix  ][iz-1]) +     \
                       2*C1*(a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Fx(a,ix,iz,s) (2*C4*(a[ix+4][iz  ] - a[ix-3][iz  ]) +     \
                       2*C3*(a[ix+3][iz  ] - a[ix-2][iz  ]) +     \
                       2*C2*(a[ix+2][iz  ] - a[ix-1][iz  ]) +     \
                       2*C1*(a[ix+1][iz  ] - a[ix  ][iz  ])  )*s

/* 2D backward FD derivative stencils */
#define Bz(a,ix,iz,s) (2*C4*(a[ix  ][iz+3] - a[ix  ][iz-4]) +     \
                       2*C3*(a[ix  ][iz+2] - a[ix  ][iz-3]) +     \
                       2*C2*(a[ix  ][iz+1] - a[ix  ][iz-2]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix  ][iz-1])  )*s
#define Bx(a,ix,iz,s) (2*C4*(a[ix+3][iz  ] - a[ix-4][iz  ]) +     \
                       2*C3*(a[ix+2][iz  ] - a[ix-3][iz  ]) +     \
                       2*C2*(a[ix+1][iz  ] - a[ix-2][iz  ]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix-1][iz  ])  )*s

/***************/
/* 3D STENCILS */
/***************/
/*  3D forward FD derivative stencils */
#define Fz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+4] - a[iy ][ix  ][iz-3]) +     \
                       2*C3*( a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-2]) +     \
                       2*C2*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-1]) +     \
                       2*C1*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz  ])  )*s
#define Fx3(a,iy,ix,iz,s) (2*C4*(a[iy][ix+4][iz  ] - a[iy ][ix-3][iz  ]) +     \
                       2*C3*( a[iy][ix+3][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       2*C2*( a[iy][ix+2][iz  ] - a[iy ][ix-1][iz  ]) +     \
                       2*C1*( a[iy][ix+1][iz  ] - a[iy ][ix  ][iz  ])  )*s
#define Fy3(a,iy,ix,iz,s) (2*C4*(a[iy+4][ix][iz  ] - a[iy-3][ix ][iz  ]) +     \
                       2*C3*( a[iy+3][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       2*C2*( a[iy+2][ix][iz  ] - a[iy-1][ix ][iz  ]) +     \
                       2*C1*( a[iy+1][ix][iz  ] - a[iy  ][ix ][iz  ])  )*s

/* 3D backward FD derivative stencils */
#define Bz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-4]) +     \
                       2*C3*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-3]) +     \
                       2*C2*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz-2]) +     \
                       2*C1*( a[iy][ix  ][iz  ] - a[iy ][ix  ][iz-1])  )*s
#define Bx3(a,iy,ix,iz,s) (2*C4*(a[iy][ix+3][iz  ] - a[iy ][ix-4][iz  ]) +     \
                       2*C3*( a[iy][ix+2][iz  ] - a[iy ][ix-3][iz  ]) +     \
                       2*C2*( a[iy][ix+1][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       2*C1*( a[iy][ix  ][iz  ] - a[iy ][ix-1][iz  ])  )*s
#define By3(a,iy,ix,iz,s) (2*C4*(a[iy+3][ix][iz  ] - a[iy-4][ix ][iz  ]) +     \
                       2*C3*( a[iy+2][ix][iz  ] - a[iy-3][ix ][iz  ]) +     \
                       2*C2*( a[iy+1][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       2*C1*( a[iy  ][ix][iz  ] - a[iy-1][ix ][iz  ])  )*s
/*------------------------------------------------------------*/


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
    fdm->nypad=sf_n(ay_)+2*fdm->nb;
	
    fdm->ozpad=sf_o(az_)-fdm->nb*fdm->dz;
    fdm->oxpad=sf_o(ax_)-fdm->nb*fdm->dx;
    fdm->oypad=sf_o(ay_)-fdm->nb*fdm->dy;

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
	int nb;
	
	nb = fdm->nb;
	
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    b[fdm->nb+ix][fdm->nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<nb;    iz++) {
	    b[ix][           iz  ] = b[ix][           fdm->nb  ];
	    b[ix][nb+fdm->nz-1+iz] = b[ix][nb+fdm->nz-1];
	}
    }

    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
	    b[           ix  ][iz] = b[           fdm->nb  ][iz];
	    b[nb+fdm->nx-1+ix][iz] = b[nb+fdm->nx-1][iz];
	}
    }
}

/*------------------------------------------------------------*/
void expand3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm)
/*< expand domain >*/
{
    int iz,ix,i3;

    for         (i3=0;i3<fdm->ny;i3++) {
	for     (ix=0;ix<fdm->nx;ix++) {
	    for (iz=0;iz<fdm->nz;iz++) {
		b[fdm->nb+i3][fdm->nb+ix][fdm->nb+iz] = a[i3][ix][iz];
	    }
	}
    }

    for         (i3=0; i3<fdm->nypad; i3++) {
	for     (ix=0; ix<fdm->nxpad; ix++) {
	    for (iz=0; iz<fdm->nb;    iz++) {
		b[i3][ix][           iz  ] = b[i3][ix][           fdm->nb  ];
		b[i3][ix][fdm->nzpad-iz-1] = b[i3][ix][fdm->nzpad-fdm->nb-1];
	    }
	}
    }


    for         (i3=0; i3<fdm->nypad; i3++) {
	for     (ix=0; ix<fdm->nb;    ix++) {
	    for (iz=0; iz<fdm->nzpad; iz++) {
		b[i3][           ix  ][iz] = b[i3][           fdm->nb  ][iz];
		b[i3][fdm->nxpad-ix-1][iz] = b[i3][fdm->nxpad-fdm->nb-1][iz];
	    }
	}
    }

    for         (i3=0; i3<fdm->nb;    i3++) {
	for     (ix=0; ix<fdm->nxpad; ix++) {
	    for (iz=0; iz<fdm->nzpad; iz++) {
		b[           i3  ][ix][iz] = b[           fdm->nb  ][ix][iz];
		b[fdm->nypad-i3-1][ix][iz] = b[fdm->nypad-fdm->nb-1][ix][iz];
	    }
	}
    }

}


/*------------------------------------------------------------*/
void cut2d(float**  a,
	   float**  b,
	   fdm2d  fdm,
	   sf_axis c1, 
	   sf_axis c2)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix;
    int fz,fx;

    fz = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);

    for     (ix=0;ix<sf_n(c2);ix++) {
	for (iz=0;iz<sf_n(c1);iz++) {
	    b[ix][iz] = a[fx+ix][fz+iz];
	}
    }
}

/*------------------------------------------------------------*/
void cut3d(float*** a,
	   float*** b,
	   fdm3d  fdm,
	   sf_axis c1, 
	   sf_axis c2,
	   sf_axis c3)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix,iy;
    int fz,fx,fy;

    fz = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);
    fy = (floor)((sf_o(c3)-fdm->oypad)/fdm->dy);

    for         (iy=0;iy<sf_n(c3);iy++) {
	for     (ix=0;ix<sf_n(c2);ix++) {
	    for (iz=0;iz<sf_n(c1);iz++) {
		b[iy][ix][iz] = a[fy+iy][fx+ix][fz+iz];
	    }
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
        sf_warning("YOUR SOURCES ARE OUTSIDE OF THE GRID!!!\n");
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
		   float *ww,
		   lint2d ca)
/*< inject into wavefield >*/
{
    int   ia;
    float wa;

#ifdef _OPENMP
#pragma omp parallel for \
    schedule(dynamic,1) \
    private(ia,wa) \
    shared(ca,ww,uu)
#endif
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void lint3d_inject(float***uu,
		   float  *ww,
		   lint3d  ca)
/*< inject into wavefield >*/
{
    int   ia;
    float wa;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ia,wa) shared(ca,ww,uu)
#endif
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w111[ia];
    }
}

/*------------------------------------------------------------*/
void lint2d_inject1(float**uu,
		    float  ww,
		    lint2d ca)
/*< inject into wavefield >*/
{
    int   ia;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,ww,uu)
#endif
    for (ia=0;ia<ca->n;ia++) {

	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void lint3d_inject1(float***uu,
		    float   ww,
		    lint3d  ca)
/*< inject into wavefield >*/
{
    int   ia;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,ww,uu)
#endif
    for (ia=0;ia<ca->n;ia++) {

	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w111[ia];
    }
}

/*------------------------------------------------------------*/
void lint2d_extract(float**uu,
		    float* dd,
		    lint2d ca)
/*< extract from wavefield >*/
{
    int ia;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,dd,uu)
#endif
    for (ia=0;ia<ca->n;ia++) {
	dd[ia] =
	    uu[ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w00[ia] +
	    uu[ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w01[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w10[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w11[ia];
    }
}  

void lint3d_extract(float***uu,
		    float  *dd,
		    lint3d  ca)
/*< extract from wavefield >*/
{
    int ia;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,dd,uu)
#endif
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
bell2d fdbell_init(int n)
/*< init bell taper >*/
{
	int   i1,i2;
	float s;
	bell2d bell;

	bell = sf_alloc(1,sizeof(*bell));

	bell->nbell = n;
	s = 0.5*bell->nbell;

	bell->bell=sf_floatalloc2(2*n+1,2*n+1);

	for		(i2=-n;i2<=n;i2++) {
		for	(i1=-n;i1<=n;i1++) {
			bell->bell[n+i2][n+i1] = exp(-(i1*i1+i2*i2)/s);
		}
	}
	return bell;
}

/*------------------------------------------------------------*/
void fdbell3d_init(int n)
/*< init bell taper >*/
{
    int   iz,ix,i3;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell3d=sf_floatalloc3(2*nbell+1,2*nbell+1,2*nbell+1);

    for        (i3=-nbell;i3<=nbell;i3++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		bell3d[nbell+i3][nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix+i3*i3)/s);
	    }
	}    
    }
}

/*------------------------------------------------------------*/
void lint2d_bell(float**uu,
		 float *ww,
		 bell2d bell,
		 lint2d ca)
/*< apply bell taper >*/
{
	int   ia,iz,ix;
	float wa;
	int nbell;
	float **b;

	nbell = bell->nbell;
	b     = bell->bell;

	for    (ix=-nbell;ix<=nbell;ix++) {
		for(iz=-nbell;iz<=nbell;iz++) {

			for (ia=0;ia<ca->n;ia++) {
				wa = ww[ia] * b[nbell+ix][nbell+iz];
				//fprintf(stderr,"\n ix+ca->jx[ia]=%d,\t iz+ca->jz[ia]=%d,\t ix=%d, ca->jx[ia]=%d\t ",ix+ca->jx[ia],iz+ca->jz[ia],ix,ca->jx[ia]);

				uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w00[ia];
				uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w01[ia];
				uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w10[ia];
				uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w11[ia];
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
		    
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w000[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w001[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w010[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w011[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w100[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w101[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w110[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w111[ia];
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
void free_abcone2d(abcone2d abc)
/*< Free the space allocated for the one-way radiating boundary >*/
{

	free(abc->bzl);
	free(abc->bzh);
	free(abc->bxl);
	free(abc->bxh);
	
	free(abc);
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

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iop)				\
    shared(fdm,nop,uo,um,abc)
#endif
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

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iop)				\
    shared(fdm,nop,uo,um,abc)
#endif
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

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
#endif

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
    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
#endif

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

#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
#endif

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
    sb = nb;  
    for(ib=0; ib<nb; ib++) {
	fb = .3*ib/sb;
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
    schedule(dynamic,1)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ibz = fdm->nzpad-ib-1;
	for(ix=0; ix<fdm->nxpad; ix++) {
//        fprintf(stderr,"Value before: %f\n",uu[ix][ib]);
	    uu[ix][ib ] *= w; /*    top sponge */
	    uu[ix][ibz] *= w; /* bottom sponge */
 //       fprintf(stderr,"Value after: %f\n",uu[ix][ib]);
	}

	ibx = fdm->nxpad-ib-1;
	for(iz=0; iz<fdm->nzpad; iz++) {
	    uu[ib ][iz] *= w; /*   left sponge */
	    uu[ibx][iz] *= w; /*  right sponge */
	}

    }
}

/*------------------------------------------------------------*/
void free_sponge(sponge spo)
/*< deallocate the sponge >*/
{
	free(spo->w);
	free(spo);
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
    schedule(dynamic,1)				\
    private(ib,iz,ix,iy,ibz,ibx,iby,w)		\
    shared(fdm,uu)
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

	ibx = fdm->nxpad-ib-1;
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		uu[iy][ib ][iz] *= w; /* x min */
		uu[iy][ibx][iz] *= w; /* x max */
	    }
	}
	
	iby = fdm->nypad-ib-1;
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		uu[ib ][ix][iz] *= w; /* x min */
		uu[iby][ix][iz] *= w; /* x max */
	    }
	}

    }
} 

float min(float x, float y){
    if (x < y) return x;
    else return y;
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

    if (dy < 0)dim = 2;  
    
    if (dim == 2) dmin = min(dx,dz); 
    else dmin = min(dx,min(dy,dz)); 
    
    float tdp = dt * vpmax *sqrt(2); // maximum distance
    if (dmin > tdp) sf_warning("CFL: Stability check ... %s-wave... PASSED", wave);
    else { 
        sf_error("CFL: Stability check ... FAILED ... minimum grid sampling: %f !> %f", dmin, tdp); 
    } 
    float wplength = safety*vpmin / fmax;
  
    bool passed;
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


