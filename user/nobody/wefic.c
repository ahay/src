#include <math.h>

#include <rsf.h>
/*^*/

#include "wefic.h"

#include "slice.h"
/*^*/

static bool verb;
static axa am,aw,ah;

static float          *qq; /* image */
static float complex **us; /* wavefield */
static float complex **ur;

static float complex **tt; /* phase shift */
/*static float complex **ts;*/
/*static float complex **tr;*/

/*------------------------------------------------------------*/

void wefic_init(bool verb_,
		axa amx_        /* i-line (data) */,
		axa amy_        /* x-line (data) */,
		axa amz_        /* depth */,
		axa aw_         /* frequency */
)
/*< initialize WE IC >*/
{
    verb=verb_;

    aw = aw_;

    am.n = amx_.n*amy_.n*amz_.n;

    /* allocate wavefield storage */
    us = sf_complexalloc2(am.n,aw.n);
    ur = sf_complexalloc2(am.n,aw.n);
    qq = sf_floatalloc   (am.n);

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;
}

void hhfic_init(axa ah_)
/*< initialize prestack IC >*/
{
    int iw,ih;
    float w;
    float h;

    ah = ah_;

    tt = sf_complexalloc2(aw.n,ah.n);
/*    ts = sf_complexalloc2(aw.n,ah.n);*/
/*    tr = sf_complexalloc2(aw.n,ah.n);*/
    
    for (ih=0; ih<ah.n; ih++) {
	h = ah.o+ih*ah.d;
	for (iw=0; iw<aw.n; iw++) {
	    w = aw.o+iw*aw.d;
	    tt[ih][iw] = cexpf(-2*I*w*h);
	    
/*	    tr[ih][iw] = cexpf(I*w*h);*/
/*	    ts[ih][iw] = conjf(tr[ih][iw]);*/
	}
    }
}

/*------------------------------------------------------------*/

void wefic_close(void)
/*< free allocated storage >*/
{
    free( *us); free( us);
    free( *ur); free( ur);
    ;           free( qq);
}

void hhfic_close(void)
/*< free allocated storage >*/
{ 
    free( *tt); free( tt);
}

/*------------------------------------------------------------*/

void zofic(fslice sdat /* source   data [nw][nz][ny][nx] */,
	   fslice rdat /* receiver data [nw][nz][ny][nx] */,
	   fslice imag /*         image     [nz][ny][nx] */
    )
/*< apply imaging condition >*/
{
    int im,iw;

    fslice_get(sdat,0,us[0]);
    fslice_get(rdat,0,ur[0]);

    for(im=0; im<am.n; im++){
	qq[im] = 0.0;
    }

    for (iw=0; iw<aw.n; iw++) {
	for(im=0; im<am.n; im++){
	    ;         qq    [im] += crealf( 
		conjf(us[iw][im]) 
		*     ur[iw][im]          ); 
	}
    }
    fslice_put(imag, 0,qq);
}

/*------------------------------------------------------------*/

void hhfic(fslice sdat /* source   data [nw][nz][ny][nx] */,
	   fslice rdat /* receiver data [nw][nz][ny][nx] */,
	   fslice imag /*         image     [nz][ny][nx] */
    )
/*< apply imaging condition >*/
{
    int im,ih,iw;

    /* read wavefields */
    fslice_get(sdat,0,us[0]);
    fslice_get(rdat,0,ur[0]);

    for(ih=0; ih<ah.n; ih++){
	if(verb) sf_warning ("ih=%3d of %3d",ih+1,ah.n);
	
	for(im=0; im<am.n; im++){
	    qq[im] = 0.0;	    
	}


	for (iw=0; iw<aw.n; iw++) {
	    for(im=0; im<am.n; im++){
		;         qq    [im] +=  crealf ( 
		    conjf(us[iw][im]) 
		    *     ur[iw][im]*tt[ih][iw] ); 
	    }
	}
/*
	for (iw=0; iw<aw.n; iw++) {
	    for(im=0; im<am.n; im++){
		;         qq    [im] +=  crealf ( 
		    conjf(us[iw][im]*ts[ih][iw]) 
		    *     ur[iw][im]*tr[ih][iw] ); 
	    }
	}
*/
	fslice_put(imag,ih,qq);
    } /* ih */
}

