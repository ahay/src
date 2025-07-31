#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
/*#include "fdutil3d.h"*/


#ifndef _fdutil3d_h

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

typedef struct spon *sponge1d;
/*^*/
typedef struct sponge *sponge2d;
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

struct spon{
    float *w;
};
/*^*/
struct sponge{
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
int omp_init(void)
/*< init OMP parameters >*/
{
#ifdef _OPENMP
    int ompath;
#endif
    int ompnth,ompchunk;
    
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OMP data chunk size */

    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
    /* OMP available threads */

#ifdef _OPENMP
#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    return ompnth;
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
    int f1,f2;

    f1 = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    f2 = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);

    for     (ix=0;ix<sf_n(c2);ix++) {
	for (iz=0;iz<sf_n(c1);iz++) {
	    b[ix][iz] = a[f2+ix][f1+iz];
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
    int iz,ix,i3;
    int f1,f2,f3;

    f1 = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    f2 = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);
    f3 = (floor)((sf_o(c3)-fdm->oypad)/fdm->dy);

    for         (i3=0;i3<sf_n(c3);i3++) {
	for     (ix=0;ix<sf_n(c2);ix++) {
	    for (iz=0;iz<sf_n(c1);iz++) {
		b[i3][ix][iz] = a[f3+i3][f2+ix][f1+iz];
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
#pragma omp parallel for schedule(dynamic,1) private(ia,wa) shared(ca,ww,uu)
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
void lint3d_bell(float***uu,
		 float  *ww,
		 lint3d  ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix,i3;
    float wa;

    for        (i3=-nbell;i3<=nbell;i3++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		
		for (ia=0;ia<ca->n;ia++) {
		    wa = ww[ia] * bell3d[nbell+i3][nbell+ix][nbell+iz];
		    
		    uu[ i3+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w000[ia];
		    uu[ i3+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w001[ia];
		    uu[ i3+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w010[ia];
		    uu[ i3+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w011[ia];
		    uu[ i3+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w100[ia];
		    uu[ i3+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w101[ia];
		    uu[ i3+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w110[ia];
		    uu[ i3+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w111[ia];
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
sponge1d sponge_make(int nb)
/*< init boundary sponge >*/
{
    sponge1d spo;
    int   ib;
    float sb,fb;
    
    spo = (sponge1d) sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(nb);

    sb = 4.0*nb;               
    for(ib=0; ib<nb; ib++) {
	fb = ib/(sqrt(2.0)*sb);
	spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

/*------------------------------------------------------------*/
void sponge2d_apply(float**   uu,
		    sponge1d spo,
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
	    uu[ix][ib ] *= w; /*    top sponge */
	    uu[ix][ibz] *= w; /* bottom sponge */
	}

	ibx = fdm->nxpad-ib-1;
	for(iz=0; iz<fdm->nzpad; iz++) {
	    uu[ib ][iz] *= w; /*   left sponge */
	    uu[ibx][iz] *= w; /*  right sponge */
	}

    }
}

/*------------------------------------------------------------*/
void sponge3d_apply(float  ***uu,
		    sponge1d spo,
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
