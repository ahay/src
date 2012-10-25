#include <rsf.h>
#include <math.h>
/*^*/
#include "eigen2x2.h"
/*^*/

/*------------------------------------------------------------*/

#define  KMAP(i,n) (i<n/2.) ? SF_PI*i/(n/2.) : SF_PI*(-n+i)/(n/2.);
#define  KMAPK(i,n) (  ( i-n/2 )/(n/2.)*SF_PI  );
/*#define  KMAP(i,n) (  ( i-n/2 )/(n/2.)*SF_PI  );*/

/* Joe's taper */
/*#define TAPER(k) (0.5*(1+cosf(k))) */

/* 2nd order */
#define TAPER2(k) (k!=0 ? sin(  k)/k : 1)
/*#define TAPER(k) (sin(  k)/k )*/
/* 4th order */
#define TAPER4(k) (k!=0 ?	       \
		  4./3.  *sin(  k)/k - \
		  1./6.  *sin(2*k)/k : 1)

/* 6th order */
#define TAPER6(k) (k!=0 ?	       \
		  3./2.  *sin(  k)/k - \
		  3./10. *sin(2*k)/k + \
		  1./60. *sin(3*k)/k : 1)

/* 8th order */
#define TAPER8(k) (k!=0 ?			\
		  8./5.  *sin(  k)/(k) -	\
		  2./5.  *sin(2*k)/(k) +	\
		  8./105.*sin(3*k)/(k) -	\
		  1./140.*sin(4*k)/(k) : 1)

/*#define TAPER(k) 1.*/

/* Dave's 2nd order */
#define TAPERD(k1,k2) (k1!=0 ? sinf(k1/2)*cosf(k2/2)/k1 :  2*cosf(k2/2) )

#define TAPERG(kmag,sig) exp(-(kmag*kmag)/2/sig/sig)/2.


#define TAPERS(k) (k!=0 ? sin(  k)/k : 1)

/*#define filter(k) sin(k)/k*/
#define TAPER(k) (k!=0 ?			\
		  (k<SF_PI ? filter(k) :0)	\
		  :1)

#define filter(k)	  8./5.  *sin(  k)/k -			\
					  2./5.  *sin(2*k)/k +	\
					  8./105.*sin(3*k)/k -  \
					  1./140.*sin(4*k)/k 
#define TAPERK(k1,k2) (k1*k1+k2*k2 != 0 ?				\
		       (sqrt(k1*k1+k2*k2)<SF_PI  ? filter(k1) :0.)	\
		       :1)

/*Butterworth 3rd order low pass filter*/
#define filterB(k)    (1+cos(k))*(1+cos(k))*(1+cos(k))/(5+3*cos(2*k))
#define filterBreal(k)   4*pow(cos(k/2),4) * (-1+3*cos(k))         /(5+3*cos(2*k))
#define filterBimag(k)   4*pow(cos(k/2),3) * ( 1+3*cos(k))*sin(k/2)/(5+3*cos(2*k))

#define TAPERB(kmag)  (kmag<SF_PI  ? filterB(kmag)  :0 )

#define TAPERBR(kmag) (kmag<SF_PI  ? filterBreal(kmag)  :0 )

#define TAPERBI(kmag) (kmag<SF_PI  ? filterBimag(kmag)  :0 ) 
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/


typedef struct wfs *wfs2d;
/*^*/
struct wfs{
    sf_complex ***temp;
    sf_fft3d   ftz,   ftx;
};/*<test>*/


/*------------------------------------------------------------*/
wfs2d wfsep_init(sf_axis  ax,
		 sf_axis  az)
/*< init wavefield separator >*/
{
    wfs2d wfs;
    wfs = (wfs2d) sf_alloc(1,sizeof(*wfs));

    wfs->temp = sf_complexalloc3(sf_n(az),sf_n(ax),1);


    wfs->ftz=sf_fft3a1_init(sf_n(az),sf_n(ax),1);
    wfs->ftx=sf_fft3a2_init(sf_n(az),sf_n(ax),1);

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

void fftshift1d(float *x, int n)
/*< test >*/
{
    float tmp;
    int i,n2;
    n2 = n / 2;    /* half of vector length */

    for (i = 0; i < n2; i++)
    {
	tmp     = x[i];
	x[i]    = x[i+n2];
	x[i+n2] = tmp;
    }   

}


void fftshift2d(float **x, int nx, int nz)
/*< test >*/
{
    
    int m2, n2;
    int i, k;
   
    float tmp13, tmp24;

    m2 = nx / 2;    /* half of row dimension */
    n2 = nz / 2;    /* half of column dimension */

/* interchange entries in 4 quadrants, 1 <--> 3 and 2 <--> 4 */

    for (i = 0; i < m2; i++)
    {
	for (k = 0; k < n2; k++)
	{
	    tmp13         = x[i][k];
	    x[i][k]       = x[i+m2][k+n2];
	    x[i+m2][k+n2] = tmp13;

	    tmp24         = x[i+m2][k];
	    x[i+m2][k]    = x[i][k+n2];
	    x[i][k+n2]    = tmp24;
	}
    }



}


/*------------------------------------------------------------*/
void wfsep(float **zdel,
	   float **xdel,
	   sf_axis  ax,
	   sf_axis  az,
	   float   c11,
	   float   c13,
	   float   c15,
	   float   c33,
	   float   c35,
	   float   c55,
	   wfs2d   wfs
	   
    )
/*< test >*/
{
  
    float a11,a12,a22;
    float uz,ux,k;
    
    int    jz,  jx;
    int    nz,  nx;
    float  kz,  kx;


    /*   float tmp; */

    extern float **tpvalue;
    extern char *tapertype; 
/*    sf_warning("tapertype= %s",tapertype);*/
    float v[2][2],A[2][2],d[2];
    extern char *domain;
       
    nz = sf_n(az);
    nx = sf_n(ax);
/*    sf_warning("nx=%d nz=%d",nx,nz);*/

    for    (jx=0;jx<nx;jx++){ kx = KMAPK(jx,nx);
	for(jz=0;jz<nz;jz++){ kz = KMAPK(jz,nz);

	    zdel[jx][jz]=kz;
	    xdel[jx][jz]=kx;
	}
    }

    for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);

	    k=sqrt(kz*kz+kx*kx);
		
	    a11=  c11*kx*kx+
		2*c15*kx*kz +
		c55*kz*kz;
		
	    a12= c15     *kx*kx +
		(c13+c55)*kx*kz +
		c35     *kz*kz;

	    a22=  c55*kx*kx+
		2*c35*kx*kz +
		c33*kz*kz;

	
	    A[0][0] = a11;
	    A[0][1] = a12;
	    A[1][0] = a12;
	    A[1][1] = a22;

	    solveSymmetric22(A,v,d);

	    ux=v[0][0]*k;
	    uz=v[0][1]*k;
	    
	    /* get the closest direction to k */
	    if(ux*kx + uz*kz <0) {
		uz = -uz;
		ux = -ux;
	    }
	  

	    xdel[jx][jz]=ux*tpvalue[jx][jz];
	    if(tapertype[0]=='s')
		zdel[jx][jz]=uz*tpvalue[jz][jx];
	    else
		zdel[jx][jz]=uz*tpvalue[jx][jz];
	    
	}
    }

    if (domain[0]=='k'){
	fftshift2d(xdel,nx,nz);
	fftshift2d(zdel,nx,nz);
    }
    else{
	/*------------------------------------------------------------*/
	/* 
	 * Z derivative
	 */
	for    (jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		/* tmp=sqrt(kx*kx+kz*kz); */
		wfs->temp[0][jx][jz]=sf_cmplx( 0, zdel[jx][jz]);
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
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		/* tmp=sqrt(kx*kx+kz*kz); */
		wfs->temp[0][jx][jz]=sf_cmplx( 0, xdel[jx][jz]);
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
    }



}

