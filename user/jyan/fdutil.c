#include <rsf.h>
#include "fdutil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _fdutil_h

typedef struct fdm2 *fdm2d;
/*^*/

typedef struct fdm3 *fdm3d;
/*^*/

typedef struct lcoef2 *lint2d;
/*^*/

typedef struct lcoef3 *lint3d;
/*^*/

typedef struct abcone *abcone2d;
/*^*/

typedef struct sponge *sponge2d;
/*^*/

typedef struct ofg *ofg2d;
/*^*/

struct fdm2{
    int nb;
    int   n1,n1pad;
    int   n2,n2pad;
    float o1,o1pad;
    float o2,o2pad;
    float d1;
    float d2;
    bool verb;
    bool free;
    int ompchunk;
};
/*^*/

struct fdm3{
    int nb;
    int   n1,n1pad;
    int   n2,n2pad;
    int   n3,n3pad;
    float o1,o1pad;
    float o2,o2pad;
    float o3,o3pad;
    float d1;
    float d2;
    float d3;
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
    int *j1;
    int *j2;
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
    int *j1;
    int *j2;
    int *j3;
};
/*^*/

struct abcone{
    bool free;
    float *b1l;
    float *b1h;
    float *b2l;
    float *b2h;
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
		  sf_axis a1_, 
		  sf_axis a2_, 
		  int     nb_,
		  int ompchunk_) 
/*< init fdm utilities >*/
{ 
    fdm2d fdm;
    fdm = (fdm2d) sf_alloc(1,sizeof(*fdm));

    fdm->free=free_;
    fdm->verb=verb_;

    fdm->nb=nb_;

    fdm->n1=sf_n(a1_);
    fdm->n2=sf_n(a2_);

    fdm->d1=sf_d(a1_);
    fdm->d2=sf_d(a2_);

    fdm->o1=sf_o(a1_);
    fdm->o2=sf_o(a2_);

    fdm->n1pad=sf_n(a1_)+2*fdm->nb;
    fdm->n2pad=sf_n(a2_)+2*fdm->nb;
	
    fdm->o1pad=sf_o(a1_)-fdm->nb*fdm->d1;
    fdm->o2pad=sf_o(a2_)-fdm->nb*fdm->d2;

    fdm->ompchunk=ompchunk_;

    return fdm;
}

/*------------------------------------------------------------*/
fdm3d fdutil3d_init(bool verb_, 
		    bool free_,
		    sf_axis a1_, 
		    sf_axis a2_, 
		    sf_axis a3_, 
		    int     nb_,
		    int ompchunk_) 
/*< init fdm utilities >*/
{ 
    fdm3d fdm;
    fdm = (fdm3d) sf_alloc(1,sizeof(*fdm));

    fdm->free=free_;
    fdm->verb=verb_;

    fdm->nb=nb_;

    fdm->n1=sf_n(a1_);
    fdm->n2=sf_n(a2_);
    fdm->n3=sf_n(a3_);

    fdm->d1=sf_d(a1_);
    fdm->d2=sf_d(a2_);
    fdm->d3=sf_d(a3_);

    fdm->o1=sf_o(a1_);
    fdm->o2=sf_o(a2_);
    fdm->o3=sf_o(a3_);

    fdm->n1pad=sf_n(a1_)+2*fdm->nb;
    fdm->n2pad=sf_n(a2_)+2*fdm->nb;
    fdm->n3pad=sf_n(a3_)+2*fdm->nb;
	
    fdm->o1pad=sf_o(a1_)-fdm->nb*fdm->d1;
    fdm->o2pad=sf_o(a2_)-fdm->nb*fdm->d2;
    fdm->o3pad=sf_o(a3_)-fdm->nb*fdm->d3;

    fdm->ompchunk=ompchunk_;

    return fdm;
}

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
ofg2d offgrid_init(fdm2d fdm)
/*< init off-grid interpolation >*/
{
    ofg2d ofg;
    ofg = (ofg2d) sf_alloc(1,sizeof(*ofg));
    
    ofg->tt = sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    
    return ofg;
}

/*------------------------------------------------------------*/
void offgridfor(float **ti,
		ofg2d  ofg,
		fdm2d  fdm)
/*< forward off-grid interpolation (in place) >*/
{
    int i1,i2;

    /* zero output */
    for     (i2=0;i2<fdm->n2;i2++) {
	for (i1=0;i1<fdm->n1;i1++) {
	    ofg->tt[i2][i1]=0;
	}
    }

    for     (i2=0;i2<fdm->n2-1;i2++) {
	for (i1=0;i1<fdm->n1-1;i1++) {
	    ofg->tt[i2][i1] =
		ti[i2  ][i1  ] + 
		ti[i2+1][i1  ] + 
		ti[i2  ][i1+1] + 
		ti[i2+1][i1+1];
	}
    }   

    /* copy to input array */
    for     (i2=0;i2<fdm->n2;i2++) {
	for (i1=0;i1<fdm->n1;i1++) {
	    ti[i2][i1] = 0.25 * ofg->tt[i2][i1];
	}
    }
}

/*------------------------------------------------------------*/
void offgridadj(float **ti,
		ofg2d  ofg,
		fdm2d  fdm)
/*< adjoint off-grid interpolation (in place) >*/
{
    int i1,i2;

    /* zero output */
    for     (i2=0;i2<fdm->n2;i2++) {
	for (i1=0;i1<fdm->n1;i1++) {
	    ofg->tt[i2][i1]=0;
	}
    }
    
    for     (i2=0;i2<fdm->n2-1;i2++) {
	for (i1=0;i1<fdm->n1-1;i1++) {
	    ofg->tt[i2  ][i1  ] += ti[i2][i1];
	    ofg->tt[i2+1][i1  ] += ti[i2][i1];
	    ofg->tt[i2  ][i1+1] += ti[i2][i1];
	    ofg->tt[i2+1][i1+1] += ti[i2][i1];
	}
    }   

    /* copy to input array */
    for     (i2=0;i2<fdm->n2;i2++) {
	for (i1=0;i1<fdm->n1;i1++) {
	    ti[i2][i1] = 0.25 * ofg->tt[i2][i1];
	}
    }
}

/*------------------------------------------------------------*/
void expand(float** a, 
	    float** b, 
	    fdm2d fdm)
/*< expand domain >*/
{
    int i1,i2;

    for     (i2=0;i2<fdm->n2;i2++) {
	for (i1=0;i1<fdm->n1;i1++) {
	    b[fdm->nb+i2][fdm->nb+i1] = a[i2][i1];
	}
    }

    for     (i2=0; i2<fdm->n2pad; i2++) {
	for (i1=0; i1<fdm->nb;    i1++) {
	    b[i2][           i1  ] = b[i2][           fdm->nb  ];
	    b[i2][fdm->n1pad-i1-1] = b[i2][fdm->n1pad-fdm->nb-1];
	}
    }

    for     (i2=0; i2<fdm->nb;    i2++) {
	for (i1=0; i1<fdm->n1pad; i1++) {
	    b[           i2  ][i1] = b[           fdm->nb  ][i1];
	    b[fdm->n2pad-i2-1][i1] = b[fdm->n2pad-fdm->nb-1][i1];
	}
    }
}

void expand3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm)
/*< expand domain >*/
{
    int i1,i2,i3;

    for         (i3=0;i3<fdm->n3;i3++) {
	for     (i2=0;i2<fdm->n2;i2++) {
	    for (i1=0;i1<fdm->n1;i1++) {
		b[fdm->nb+i3][fdm->nb+i2][fdm->nb+i1] = a[i3][i2][i1];
	    }
	}
    }

    for         (i3=0; i3<fdm->n3pad; i3++) {
	for     (i2=0; i2<fdm->n2pad; i2++) {
	    for (i1=0; i1<fdm->nb;    i1++) {
		b[i3][i2][           i1  ] = b[i3][i2][           fdm->nb  ];
		b[i3][i2][fdm->n1pad-i1-1] = b[i3][i2][fdm->n1pad-fdm->nb-1];
	    }
	}
    }


    for         (i3=0; i3<fdm->n3pad; i3++) {
	for     (i2=0; i2<fdm->nb;    i2++) {
	    for (i1=0; i1<fdm->n1pad; i1++) {
		b[i3][           i2  ][i1] = b[i3][           fdm->nb  ][i1];
		b[i3][fdm->n2pad-i2-1][i1] = b[i3][fdm->n2pad-fdm->nb-1][i1];
	    }
	}
    }

    for         (i3=0; i3<fdm->nb;    i3++) {
	for     (i2=0; i2<fdm->n2pad; i2++) {
	    for (i1=0; i1<fdm->n1pad; i1++) {
		b[           i3  ][i2][i1] = b[           fdm->nb  ][i2][i1];
		b[fdm->n3pad-i3-1][i2][i1] = b[fdm->n3pad-fdm->nb-1][i2][i1];
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
    int i1,i2;
    int f1,f2;

    f1 = (floor)((sf_o(c1)-fdm->o1pad)/fdm->d1);
    f2 = (floor)((sf_o(c2)-fdm->o2pad)/fdm->d2);

    for     (i2=0;i2<sf_n(c2);i2++) {
	for (i1=0;i1<sf_n(c1);i1++) {
	    b[i2][i1] = a[f2+i2][f1+i1];
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
    int i1,i2,i3;
    int f1,f2,f3;

    f1 = (floor)((sf_o(c1)-fdm->o1pad)/fdm->d1);
    f2 = (floor)((sf_o(c2)-fdm->o2pad)/fdm->d2);
    f3 = (floor)((sf_o(c3)-fdm->o3pad)/fdm->d3);

    for         (i3=0;i3<sf_n(c3);i3++) {
	for     (i2=0;i2<sf_n(c2);i2++) {
	    for (i1=0;i1<sf_n(c1);i1++) {
		b[i3][i2][i1] = a[f3+i3][f2+i2][f1+i1];
	    }
	}
    }
}

/*------------------------------------------------------------*/
void bfill(float** b, 
	   fdm2d fdm)
/*< fill boundaries >*/
{
    int i1,i2;
    
    for     (i2=0; i2<fdm->n2pad; i2++) {
	for (i1=0; i1<fdm->nb;    i1++) {
	    b[i2][           i1  ] = b[i2][           fdm->nb  ];
	    b[i2][fdm->n1pad-i1-1] = b[i2][fdm->n1pad-fdm->nb-1];
	}
    }
    
    for     (i2=0; i2<fdm->nb;    i2++) {
	for (i1=0; i1<fdm->n1pad; i1++) {
	    b[           i2  ][i1] = b[           fdm->nb  ][i1];
	    b[fdm->n2pad-i2-1][i1] = b[fdm->n2pad-fdm->nb-1][i1];
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

    ca->j1  = sf_intalloc(na);
    ca->j2  = sf_intalloc(na);

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->o1pad && 
	   aa[ia].z <  fdm->o1pad + (fdm->n1pad-1)*fdm->d1 &&
	   aa[ia].x >= fdm->o2pad && 
	   aa[ia].x <  fdm->o2pad + (fdm->n2pad-1)*fdm->d2   ) {
	    
	    ca->j1[ia] = (int)( (aa[ia].z-fdm->o1pad)/fdm->d1);
	    ca->j2[ia] = (int)( (aa[ia].x-fdm->o2pad)/fdm->d2);
	    
	    f1 = (aa[ia].z-fdm->o1pad)/fdm->d1 - ca->j1[ia];
	    f2 = (aa[ia].x-fdm->o2pad)/fdm->d2 - ca->j2[ia];
	} else {
	    ca->j1[ia] = 0; 
	    ca->j2[ia] = 0;
	    
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

    ca->j1  = sf_intalloc(na);
    ca->j2  = sf_intalloc(na);
    ca->j3  = sf_intalloc(na);

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->o1pad && 
	   aa[ia].z <  fdm->o1pad + (fdm->n1pad-1)*fdm->d1 &&
	   aa[ia].x >= fdm->o2pad && 
	   aa[ia].x <  fdm->o2pad + (fdm->n2pad-1)*fdm->d2 &&   
	   aa[ia].y >= fdm->o3pad && 
	   aa[ia].y <  fdm->o3pad + (fdm->n3pad-1)*fdm->d3  ) {
	    
	    ca->j1[ia] = (int)( (aa[ia].z-fdm->o1pad)/fdm->d1);
	    ca->j2[ia] = (int)( (aa[ia].x-fdm->o2pad)/fdm->d2);
	    ca->j3[ia] = (int)( (aa[ia].y-fdm->o3pad)/fdm->d3);
	    
	    f1 = (aa[ia].z-fdm->o1pad)/fdm->d1 - ca->j1[ia];
	    f2 = (aa[ia].x-fdm->o2pad)/fdm->d2 - ca->j2[ia];
	    f3 = (aa[ia].y-fdm->o3pad)/fdm->d3 - ca->j3[ia];

	} else {
	    ca->j1[ia] = 0; 
	    ca->j2[ia] = 0;
	    ca->j3[ia] = 0;
	    
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
	
	uu[ ca->j2[ia]   ][ ca->j1[ia]   ] = wa;
	uu[ ca->j2[ia]   ][ ca->j1[ia]+1 ] = wa;
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]   ] = wa;
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]+1 ] = wa;
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
	
	uu[ ca->j2[ia]   ][ ca->j1[ia]   ] -= wa * ca->w00[ia];
	uu[ ca->j2[ia]   ][ ca->j1[ia]+1 ] -= wa * ca->w01[ia];
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]   ] -= wa * ca->w10[ia];
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]+1 ] -= wa * ca->w11[ia];
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

	uu[ ca->j2[ia]   ][ ca->j1[ia]   ] -= ww * ca->w00[ia];
	uu[ ca->j2[ia]   ][ ca->j1[ia]+1 ] -= ww * ca->w01[ia];
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]   ] -= ww * ca->w10[ia];
	uu[ ca->j2[ia]+1 ][ ca->j1[ia]+1 ] -= ww * ca->w11[ia];
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
	    uu[ ca->j2[ia]  ][ ca->j1[ia]  ] * ca->w00[ia] +
	    uu[ ca->j2[ia]  ][ ca->j1[ia]+1] * ca->w01[ia] +
	    uu[ ca->j2[ia]+1][ ca->j1[ia]  ] * ca->w10[ia] +
	    uu[ ca->j2[ia]+1][ ca->j1[ia]+1] * ca->w11[ia];
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
	    uu[ ca->j3[ia]  ][ ca->j2[ia]  ][ ca->j1[ia]  ] * ca->w000[ia] +
	    uu[ ca->j3[ia]  ][ ca->j2[ia]  ][ ca->j1[ia]+1] * ca->w001[ia] +
	    uu[ ca->j3[ia]  ][ ca->j2[ia]+1][ ca->j1[ia]  ] * ca->w010[ia] +
	    uu[ ca->j3[ia]  ][ ca->j2[ia]+1][ ca->j1[ia]+1] * ca->w011[ia] +
	    uu[ ca->j3[ia]+1][ ca->j2[ia]  ][ ca->j1[ia]  ] * ca->w100[ia] +
	    uu[ ca->j3[ia]+1][ ca->j2[ia]  ][ ca->j1[ia]+1] * ca->w101[ia] +
	    uu[ ca->j3[ia]+1][ ca->j2[ia]+1][ ca->j1[ia]  ] * ca->w110[ia] +
	    uu[ ca->j3[ia]+1][ ca->j2[ia]+1][ ca->j1[ia]+1] * ca->w111[ia];
    }
}  


/*------------------------------------------------------------*/
void fdbell_init(int n)
/*< init bell taper >*/
{
    int   i1,i2;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell=sf_floatalloc2(2*nbell+1,2*nbell+1);

    for    (i2=-nbell;i2<=nbell;i2++) {
	for(i1=-nbell;i1<=nbell;i1++) {
	    bell[nbell+i2][nbell+i1] = exp(-(i1*i1+i2*i2)/s);
	}
    }    
}

/*------------------------------------------------------------*/
void fdbell3d_init(int n)
/*< init bell taper >*/
{
    int   i1,i2,i3;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell3d=sf_floatalloc3(2*nbell+1,2*nbell+1,2*nbell+1);

    for        (i3=-nbell;i3<=nbell;i3++) {
	for    (i2=-nbell;i2<=nbell;i2++) {
	    for(i1=-nbell;i1<=nbell;i1++) {
		bell3d[nbell+i3][nbell+i2][nbell+i1] = exp(-(i1*i1+i2*i2+i3*i3)/s);
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
    int   ia,i1,i2;
    float wa;

    for    (i2=-nbell;i2<=nbell;i2++) {
	for(i1=-nbell;i1<=nbell;i1++) {
	    
	    for (ia=0;ia<ca->n;ia++) {
		wa = ww[ia] * bell[nbell+i2][nbell+i1];

		uu[ i2+ca->j2[ia]   ][ i1+ca->j1[ia]   ] -= wa * ca->w00[ia];
		uu[ i2+ca->j2[ia]   ][ i1+ca->j1[ia]+1 ] -= wa * ca->w01[ia];
		uu[ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]   ] -= wa * ca->w10[ia];
		uu[ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]+1 ] -= wa * ca->w11[ia];
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
    int   ia,i1,i2,i3;
    float wa;

    for        (i3=-nbell;i3<=nbell;i3++) {
	for    (i2=-nbell;i2<=nbell;i2++) {
	    for(i1=-nbell;i1<=nbell;i1++) {
		
		for (ia=0;ia<ca->n;ia++) {
		    wa = ww[ia] * bell3d[nbell+i3][nbell+i2][nbell+i1];
		    
		    uu[ i3+ca->j3[ia]   ][ i2+ca->j2[ia]   ][ i1+ca->j1[ia]   ] -= wa * ca->w000[ia];
		    uu[ i3+ca->j3[ia]   ][ i2+ca->j2[ia]   ][ i1+ca->j1[ia]+1 ] -= wa * ca->w001[ia];
		    uu[ i3+ca->j3[ia]   ][ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]   ] -= wa * ca->w010[ia];
		    uu[ i3+ca->j3[ia]   ][ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]+1 ] -= wa * ca->w011[ia];
		    uu[ i3+ca->j3[ia]+1 ][ i2+ca->j2[ia]   ][ i1+ca->j1[ia]   ] -= wa * ca->w100[ia];
		    uu[ i3+ca->j3[ia]+1 ][ i2+ca->j2[ia]   ][ i1+ca->j1[ia]+1 ] -= wa * ca->w101[ia];
		    uu[ i3+ca->j3[ia]+1 ][ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]   ] -= wa * ca->w110[ia];
		    uu[ i3+ca->j3[ia]+1 ][ i2+ca->j2[ia]+1 ][ i1+ca->j1[ia]+1 ] -= wa * ca->w111[ia];
		}
		
	    }
	}
    }
}

/*------------------------------------------------------------*/
abcone2d abcone2d_make(int     nop,
		       float    dt,
		       float**  vp,
		       bool   free, 
		       fdm2d   fdm)
/*< init 2D ABC >*/
{
    abcone2d abc;
    int i1,i2;
    float d;

    abc = (abcone2d) sf_alloc(1,sizeof(*abc));

    abc->free = free;

    abc->b1l = sf_floatalloc(fdm->n2pad);
    abc->b1h = sf_floatalloc(fdm->n2pad);
    abc->b2l = sf_floatalloc(fdm->n1pad);
    abc->b2h = sf_floatalloc(fdm->n1pad);

    for (i2=0;i2<fdm->n2pad;i2++) {
	d = vp[i2][           nop  ] *dt/fdm->d1; abc->b1l[i2] = (1-d)/(1+d);
	d = vp[i2][fdm->n1pad-nop-1] *dt/fdm->d1; abc->b1h[i2] = (1-d)/(1+d);
    }
    for (i1=0;i1<fdm->n1pad;i1++) {
	d = vp[           nop  ][i1] *dt/fdm->d2; abc->b2l[i1] = (1-d)/(1+d);
	d = vp[fdm->n2pad-nop-1][i1] *dt/fdm->d2; abc->b2h[i1] = (1-d)/(1+d);
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
    int i1,i2,iop;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(i1,i2,iop) shared(fdm,nop,uo,um,abc)
#endif
    for(i2=0;i2<fdm->n2pad;i2++) {
	for(iop=0;iop<nop;iop++) {

	    /* top BC */
	    if(!abc->free) { /* not free surface, apply ABC */
		i1 = nop-iop;
		uo      [i2][i1  ] 
		    = um[i2][i1+1] 
		    +(um[i2][i1  ]
		      - uo[i2][i1+1]) * abc->b1l[i2];
	    }

	    /* bottom BC */
	    i1 = fdm->n1pad-nop+iop-1;
	    uo      [i2][i1  ] 
		= um[i2][i1-1]
		+(um[i2][i1  ]
		- uo[i2][i1-1]) * abc->b1h[i2];
	}
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(i1,i2,iop) shared(fdm,nop,uo,um,abc)
#endif
    for(i1=0;i1<fdm->n1pad;i1++) {
	for(iop=0;iop<nop;iop++) {

	    /* left BC */
	    i2 = nop-iop;
	    uo      [i2  ][i1] 
		= um[i2+1][i1] 
		+(um[i2  ][i1]
		- uo[i2+1][i1]) * abc->b2l[i1];

	    /* right BC */
	    i2 = fdm->n2pad-nop+iop-1;
	    uo      [i2  ][i1] 
		= um[i2-1][i1]
		+(um[i2  ][i1]
		- uo[i2-1][i1]) * abc->b2h[i1];
	}
    }
}

/*------------------------------------------------------------*/
sponge2d sponge2d_make(fdm2d fdm)
/*< init boundary sponge >*/
{
    sponge2d spo;
    int   ib;
    float sb,fb;
    
    spo = (sponge2d) sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(fdm->nb);

    sb = 4.0*fdm->nb;               /*                 sigma */
    for(ib=0; ib<fdm->nb; ib++) {
	fb = ib/(sqrt(2.0)*sb);     /*  x / (sqrt(2) * sigma) */
	spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

/*------------------------------------------------------------*/
void sponge2d_apply(float**   uu,
		    sponge2d spo,
		    fdm2d    fdm)
/*< apply boundary sponge >*/
{
    int i1,i2,ib,ib1,ib2;
    float w;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ib,i1,i2,ib1,ib2,w) shared(fdm,uu)
#endif
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[fdm->nb-ib-1];

	ib1 = fdm->n1pad-ib-1;
	for(i2=0; i2<fdm->n2pad; i2++) {
	    uu[i2][ib ] *= w; /*    top sponge */
	    uu[i2][ib1] *= w; /* bottom sponge */
	}

	ib2 = fdm->n2pad-ib-1;
	for(i1=0; i1<fdm->n1pad; i1++) {
	    uu[ib ][i1] *= w; /*   left sponge */
	    uu[ib2][i1] *= w; /*  right sponge */
	}

    }
}

