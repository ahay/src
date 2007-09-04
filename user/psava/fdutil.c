#include <rsf.h>
#include "fdutil.h"

#ifndef _fdutil_h

typedef struct fdm *fdm2d;
/*^*/

typedef struct lcoef *lint2d;
/*^*/

typedef struct abcone *abcone2d;
/*^*/

typedef struct sponge *sponge2d;
/*^*/

struct fdm{
    int n1,n1pad;
    int n2,n2pad;
    int nb;
    float o1,o1pad;
    float o2,o2pad;
    float d1;
    float d2;
    bool verb;
    bool free;
    int ompchunk;
};
/*^*/

struct lcoef{
    int n;
    float *w00;
    float *w01;
    float *w10;
    float *w11;
    int *j1;
    int *j2;
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

#endif

static float** bell;
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
	   aa[ia].x <  fdm->o2pad + (fdm->n2pad-1)*fdm->d2) {
	    
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
