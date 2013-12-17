#include <rsf.h>
#include "fd3util.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _fd3util_h

typedef struct fdm3 *fdm3d;
/*^*/


typedef struct spon *sponge;
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

struct spon{
    float *w;
};
/*^*/

#endif

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
void expand3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm)
/*< expand domain a to b >*/
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

void window3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm)
/*< window b domain to a >*/
{
    for         (int i3=0;i3<fdm->ny;i3++) {
	for     (int ix=0;ix<fdm->nx;ix++) {
	    for (int iz=0;iz<fdm->nz;iz++) {
		a[i3][ix][iz]=b[fdm->nb+i3][fdm->nb+ix][fdm->nb+iz];
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
		uu[ib ][ix][iz] *= w; /* y min */
		uu[iby][ix][iz] *= w; /* y max */
	    }
	}

    }
}

static	float c0, c11, c12, c21, c22, c31, c32;

void fd3_init(fdm3d fdm)
/*< initialize >*/
{
	float t;

	t = 1.0/(fdm->dz*fdm->dz);
	c11 = 4.0*t/3.0;
	c12=  -t/12.0;

	t = 1.0/(fdm->dx*fdm->dx);
	c21 = 4.0*t/3.0;
	c22=  -t/12.0;

	t = 1.0/(fdm->dy*fdm->dy);
	c31 = 4.0*t/3.0;
	c32=  -t/12.0;

	c0  = -2.0 * (c11 + c12 + c21 + c22 +c31 + c32);
}

void step_forward(float ***u0, float ***u1, float ***vv, fdm3d fdm)
/*< step forward >*/
{
#ifdef _OPENMP
#pragma omp for	schedule(dynamic,fdm->ompchunk)	
#endif
	for(int iy=2; iy<fdm->nypad-2; iy++) 
	for(int ix=2; ix<fdm->nxpad-2; ix++)
	for(int iz=2; iz<fdm->nzpad-2; iz++) 
	{			    
		/* 4th order Laplacian operator */
		float ua=c0 * u1[iy  ][ix  ][iz  ] + 
			c11*(u1[iy  ][ix  ][iz-1] + u1[iy  ][ix  ][iz+1]) +
			c12*(u1[iy  ][ix  ][iz-2] + u1[iy  ][ix  ][iz+2]) +
			c21*(u1[iy  ][ix-1][iz  ] + u1[iy  ][ix+1][iz  ]) +
			c22*(u1[iy  ][ix-2][iz  ] + u1[iy  ][ix+2][iz  ]) +
			c31*(u1[iy-1][ix  ][iz  ] + u1[iy+1][ix  ][iz  ]) +
			c32*(u1[iy-2][ix  ][iz  ] + u1[iy+2][ix  ][iz  ]) ;
		u0[iy][ix][iz] = 2*u1[iy][ix][iz]-u0[iy][ix][iz]+ua * vv[iy][ix][iz];
	}
}


void add_source(float ***u1, float *wlt, int **Szxy, int ns, fdm3d fdm, bool add)
/*< add source >*/
{
	for (int is=0; is<ns; is++)
	{
		int sz=Szxy[is][0]+fdm->nb;
		int sx=Szxy[is][1]+fdm->nb;
		int sy=Szxy[is][2]+fdm->nb;

		if (add) u1[sy][sx][sz] += wlt[is];
		else 	u1[sy][sx][sz] -= wlt[is];
	}
}

void free_surface(float ***u0, float ***u1, fdm3d fdm)
/*< handle free surface at the top >*/
{
	if(fdm->free) /* free surface */
	{ 	    	
		for(int iy=0; iy<fdm->nypad; iy++) 
		for(int ix=0; ix<fdm->nxpad; ix++)
		for(int iz=0; iz<fdm->nb; iz++) 
		{
			u0[iy][ix][iz]=0;
			u1[iy][ix][iz]=0;
		}
	}
}
