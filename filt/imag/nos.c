/* Narrow-offset space-domain shift */

#include <math.h>
#include <rsf.h>
#include "fft2.h"
#include "nos.h"

#define X2K(a,b,p) b.n=a.n+p; \
                   b.d=2.0*SF_PI/(b.n*a.d); \
                   b.o=(1==b.n)?0:-SF_PI/a.d;
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);

#define KOOP(a) for(ih=0;ih<bh.n;ih++){ \
                for(im=0;im<bm.n;im++){ {a} }}
#define LOOP(a) for(ih=0;ih<ah.n;ih++){ \
                for(im=0;im<am.n;im++){ {a} }}

/*------------------------------------------------------------*/
static axa am,ah,aw,az;
static axa bm,bh;

static float complex **wk,**uk;

static float          *km;

static kiss_fft_cfg forw1, invs1;
static float        fftscale;
/*------------------------------------------------------------*/

void nos_init(axa am_,
	      axa ah_,
	      axa aw_,
	      axa az_)
/*< initialize >*/
{
    int   im;
    int   jm;

    am = am_;
    ah = ah_;
    aw = aw_;
    az = az_;

    X2K(am,bm,0);
    X2K(ah,bh,0);
    fft2_init(bm.n,bh.n);

    wk = sf_complexalloc2(bm.n,bh.n);
    uk = sf_complexalloc2(bm.n,bh.n);

    /* precompute wavenumbers */
    km = sf_floatalloc(bm.n);
    
    for (im=0; im<bm.n; im++) {
	jm = KMAP(im,bm.n);
	km[im] = bm.o + jm*bm.d;
    }

    nosfft_init();
} 

void nos(float            w,
	 float           sz,
	 float complex **wx)
/*< extrapolate >*/
{
    int ih,im;
    int jh;

    float ho,hh,dh;
    float         h_;     /* H     */
    float complex w_;     /* W     */
    float         k_;     /* K, K' */
    float complex kp;
    float complex s_, sp; /* S, S' */
    float complex ss, ds;
    float complex a,b,c;

    float complex cc;
    float complex w2;
    float         s2;

    int jn;
/*------------------------------------------------------------*/

    w2 = 0.0001*aw.d - I*w;
    w2 = w2*w2;
    s2 = sz*sz;

    KOOP( wk[ih][im] = 0.; );

    LOOP( wk[ih][im] = 
	  wx[ih][im]; );

/*    sf_warning("%f %f",crealf(wk[0][20]),cimagf(wk[0][20]));*/
    nosfft(false,wk);
    
/*------------------------------------------------------------*/
/* begin K-domain computation */

    w_ = 4*w2*s2;

    /* loop over km */
    for(im=0; im<bm.n; im++) {
	k_ = km[im]*km[im];  /* K  */
	
	/* output (h) loop */
	for(ih=0; ih<ah.n; ih++) {
	    hh = ah.o + ih*ah.d;
	    
	    uk[ih][im] = 0;

	    /* input (ho) loop */
	    for(jh=0; jh<ah.n; jh++) {
		ho = ah.o + jh*ah.d;
		dh = hh - ho;

		h_ = (dh*dh)/(az.d*az.d);    /* H  */
		kp = k_ * ((1+h_)/w_);       /* K' = K (1+H)/W */
		
		/*   Newton iterations for S */		
		sp = 0.5*(1+csqrtf(1+8*h_*kp)); /* So' */
		ss = 0;
		ds = cabs(sp-ss);
		jn=0;
		while(cabs(ds)>0.001*cabs(sp) && jn<10 ) {
		    a = h_*(1+h_)*kp*kp;
		    b = 2*h_*kp + (2-3*sp)*sp;
		    c = 4*h_*kp + (3-4*sp)*sp;
		    
		    sp = (a + b * sp*sp)/(c * sp);

		    ds = cabs(sp-ss);
		    ss = sp;

/*		    sf_warning("%d (%f,%f)",jn,crealf(sp),cimagf(sp));*/

		    jn++;
		} /* Newton iterations for S */
		s_ = sp / ((1+h_)/w_);
		
		cc = csqrtf(s_ + k_);
		sf_warning("%d %f %f",im,crealf(cc),cimagf(cc));
		uk[ih][im] += wk[jh][im] * cexpf(-cc*az.d);
	    } /* ho */
	    uk[ih][im] /= ah.n;
	} /* hh */

    } /* km */    

/* end K-domain computation */
/*------------------------------------------------------------*/
    nosfft(true,uk);
    sf_warning("%f %f",crealf(uk[0][20]),cimagf(uk[0][20]));

    LOOP( wx[ih][im] = 
	  uk[ih][im]; );
}


/*------------------------------------------------------------*/

void nos2(float            w,
	 float           sz,
	 float complex **wx)
/*< extrapolate >*/
{
    int ih,im;
    int jh;

    float ho,hh,dh;
    float h_;     /* H     */
    float w_;     /* W     */
    float k_;     /* K, K' */
    float kp;
    float s_, sp; /* S, S' */
    float ss, ds;
    float a,b,c;

    float complex cc;
    float w2;
    float s2;

    int jn;
/*------------------------------------------------------------*/

    w2 =  w* w;
    s2 = sz*sz;

    KOOP( wk[ih][im] = 0.; );

    LOOP( wk[ih][im] = 
	  wx[ih][im]; );

/*    sf_warning("%f %f",crealf(wk[0][20]),cimagf(wk[0][20]));*/
    nosfft(false,wk);
    
/*------------------------------------------------------------*/
/* begin K-domain computation */

    w_ = 4*w2*s2;

    /* loop over km */
    for(im=0; im<bm.n; im++) {
	k_ = km[im]*km[im];  /* K  */
	
	/* output (h) loop */
	for(ih=0; ih<ah.n; ih++) {
	    hh = ah.o + ih*ah.d;
	    
	    uk[ih][im] = 0;

	    /* input (ho) loop */
	    for(jh=0; jh<ah.n; jh++) {
		ho = ah.o + jh*ah.d;
		dh = hh - ho;

		h_ = (dh*dh)/(az.d*az.d);    /* H  */
		kp = k_ * ((1+h_)/w_);       /* K' = K (1+H)/W */
		
		/*   Newton iterations for S */		
		sp = 0.5*(1+sqrtf(1+8*h_*kp)); /* So' */
		ss = 0;
		ds = abs(sp-ss);
		jn=0;
		while(abs(ds)>0.001*abs(sp) && jn<10 ) {
		    a = h_*(1+h_)*kp*kp;
		    b = 2*h_*kp + (2-3*sp)*sp;
		    c = 4*h_*kp + (3-4*sp)*sp;
		    
		    sp = (a + b * sp*sp)/(c * sp);

		    ds = abs(sp-ss);
		    ss = sp;

/*		    sf_warning("%d (%f,%f)",jn,crealf(sp),cimagf(sp));*/

		    jn++;
		} /* Newton iterations for S */
		s_ = sp / ((1+h_)/w_);
		
		if(s_>k_) {
		    cc = -I*sqrtf(s_ - k_);
		} else {
		    cc=0;
		}
		
		sf_warning("%d %f %f",im,crealf(cc),cimagf(cc));
		uk[ih][im] += wk[jh][im] * cexpf(-cc*az.d);
	    } /* ho */
	    uk[ih][im] /= ah.n;
	} /* hh */

    } /* km */    

/* end K-domain computation */
/*------------------------------------------------------------*/
    nosfft(true,uk);
    sf_warning("%f %f",crealf(uk[0][20]),cimagf(uk[0][20]));

    LOOP( wx[ih][im] = 
	  uk[ih][im]; );
}


/*------------------------------------------------------------*/

void nosfft_init()
/*< initialize fft >*/
{
    forw1 = kiss_fft_alloc(am.n,0,NULL,NULL);
    invs1 = kiss_fft_alloc(am.n,1,NULL,NULL);
    
    if (NULL == forw1 || NULL == invs1)
	sf_error("%s: KISS FFT allocation error",__FILE__);
    
    fftscale = 1./am.n;
}

void nosfft_close(void)
/*< Free allocated storage >*/
{
    free (forw1);
    free (invs1);
}

void nosfft(bool inv,
	    complex float **pp /* [n1][n2] */) 
/*< Apply 1-D FFT >*/
{
    int im,ih;

    if (inv) {
	/* IFT 1 */
	for (ih=0; ih < ah.n; ih++) {
	    kiss_fft(invs1,
		     (const kiss_fft_cpx *) pp[ih], 
		     (      kiss_fft_cpx *) pp[ih]);
	}
	/* scaling */
	for (ih=0; ih < ah.n; ih++) {
	    for (im=0; im < am.n; im++) {
		pp[ih][im] *= fftscale;
	    }
	}
    } else {
	/* FFT 1 */
	for (ih=0; ih < ah.n; ih++) {
	    kiss_fft(forw1,
		     (const kiss_fft_cpx *) pp[ih], 
		     (      kiss_fft_cpx *) pp[ih]);
	}
    }
}
