#include <rsf.h>
#include <math.h>
#include "eigv.h"

#define  KMAP(i,n) (i<n/2.) ? SF_PI*i/(n/2.) : SF_PI*(-n+i)/(n/2.);

/* Joe's taper */
/*#define TAPER(k) (0.5*(1+cosf(k))) */

/* 2nd order */
/*#define TAPER(k) (k!=0 ? sinf(  k)/k : 1)*/

/* 4th order */
/*#define TAPER(k) (k!=0 ?	       \*/
/*		  4./3.  *sin(  k)/k - \*/
/*		  1./6.  *sin(2*k)/k : 1)*/

/* 6th order */
/*#define TAPER(k) (k!=0 ?	       \*/
/*		  3./2.  *sin(  k)/k - \*/
/*		  3./10. *sin(2*k)/k + \*/
/*		  1./60. *sin(3*k)/k : 1)*/

/* 8th order */
#define TAPER(k) (k!=0 ?	       \
		  8./5.  *sin(  k)/k - \
		  2./5.  *sin(2*k)/k +		\
		  8./105.*sin(3*k)/k -		\
		  1./140.*sin(4*k)/k : 1)

/*#define TAPER(k) 1.*/
/*------------------------------------------------------------*/

typedef struct wfs *wfs2d;
/*^*/

struct wfs{
    sf_complex ***temp;
    float        *eigv; /* eigenvalues */
    float        *polr; /* eigenvector */ 
    float       **ctfl; /* Christoffel matrix */
    float **UPz, **UPx;
    sf_fft3d   ftz,   ftx;
};
/*^*/


/*------------------------------------------------------------*/
wfs2d wfsep_init(sf_axis  ax,
		 sf_axis  az)
/*< init wavefield separator >*/
{
    wfs2d wfs;
    wfs = (wfs2d) sf_alloc(1,sizeof(*wfs));

    wfs->temp = sf_complexalloc3(sf_n(az),sf_n(ax),1);
    wfs->eigv = sf_floatalloc (2);           /* eigenvalues */
    wfs->polr = sf_floatalloc (2);           /* eigenvector */
    wfs->ctfl = sf_floatalloc2(2,2);         /* Christoffel matrix */

    wfs->UPz  = sf_floatalloc2(sf_n(az),sf_n(ax));
    wfs->UPx  = sf_floatalloc2(sf_n(az),sf_n(ax));

    wfs->ftz = sf_fft3a1_init(sf_n(az),sf_n(ax),1);
    wfs->ftx = sf_fft3a2_init(sf_n(az),sf_n(ax),1);

    return wfs;
}

/*------------------------------------------------------------*/
void wfsep_close(wfs2d wfs)
/*< close wavefield separator >*/
{
    sf_fft3a1_close(wfs->ftz);
    sf_fft3a2_close(wfs->ftx);
}

/*------------------------------------------------------------*/
void wfsep(float **zdel,
	   float **xdel,
	   sf_axis  ax,
	   sf_axis  az,
	   float   c11,
	   float   c33,
	   float   c55,
	   float   c13,
	   wfs2d   wfs
    )
/*< wavefield separator >*/
{
    
    float a11,a12,a22;
    float uz,ux,k;
    
    int    jz,  jx;
    int    nz,  nx;
    float  kz,  kx;

    nz = sf_n(az);
    nx = sf_n(ax);

    for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
	    wfs->UPz[jx][jz]=kz;
	    wfs->UPx[jx][jz]=kx;
	}
    }

    /* eigenvalues and eigenvectors */
    for    (jx=1;jx<nx;jx++){ kx = KMAP(jx,nx); 
	for(jz=1;jz<nz;jz++){ kz = KMAP(jz,nz);
	    k=sqrt(kz*kz+kx*kx);

	    a11=  c11*kx*kx+
		  c55*kz*kz;
	    
	    a12= (c13+c55)*kx*kz;

	    a22=  c55*kx*kx+
		  c33*kz*kz;

	    wfs->ctfl[0][0] = a11;
	    wfs->ctfl[0][1] = a12;
	    wfs->ctfl[1][0] = a12;
	    wfs->ctfl[1][1] = a22;

	    eigval2(wfs->ctfl,wfs->eigv);
	    eigvec2(wfs->ctfl,wfs->eigv[0],wfs->polr); /* use largest eigenvalue (qP mode) */

	    ux=wfs->polr[0]*k;
	    uz=wfs->polr[1]*k;

	    /* get the closest direction to k */
	    if(ux*kx + uz*kz <0) {
		wfs->UPz[jx][jz] = -uz;
		wfs->UPx[jx][jz] = -ux;
	    } else {
		wfs->UPz[jx][jz] =  uz;
		wfs->UPx[jx][jz] =  ux;
	    }

	}
    }

    /* boundaries */
    for(jx=0;jx<nx;jx++) {
	wfs->UPz[jx][0] = 0.5*(wfs->UPz[jx][1] + wfs->UPz[jx][nz-1]);
    }
    for(jz=0;jz<nz;jz++) {
	wfs->UPx[0][jz] = 0.5*(wfs->UPx[1][jz] + wfs->UPx[nx-1][jz]);
    }
    
    /*------------------------------------------------------------*/
    /* 
     * Z derivative
     */
    for    (jx=0;jx<nx;jx++){
	for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
	    wfs->temp[0][jx][jz]=sf_cmplx( 0, wfs->UPz[jx][jz] * TAPER(kz) );
	}
    }   

    sf_cnt3a2(wfs->temp,wfs->ftx);
    sf_cnt3a1(wfs->temp,wfs->ftz);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs->temp,wfs->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs->temp,wfs->ftz);
    for(    jx=0;jx<nx;jx++){
	for(jz=0;jz<nz;jz++){
	    zdel[jx][jz]=crealf(wfs->temp[0][jx][jz]);
	}
    }
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* 
     * X derivative
     */
    for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	for(jz=0;jz<nz;jz++){
	    wfs->temp[0][jx][jz]=sf_cmplx( 0, wfs->UPx[jx][jz] * TAPER(kx) );
	}
    }   

    sf_cnt3a2(wfs->temp,wfs->ftx);
    sf_cnt3a1(wfs->temp,wfs->ftz);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs->temp,wfs->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs->temp,wfs->ftz);
    for(    jx=0;jx<nx;jx++){
	for(jz=0;jz<nz;jz++){
	    xdel[jx][jz]=crealf(wfs->temp[0][jx][jz]);
	}
    }
    /*------------------------------------------------------------*/

}
