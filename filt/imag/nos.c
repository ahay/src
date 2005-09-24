/* Narrow-offset space-domain shift */

#include <math.h>
#include <rsf.h>
#include "nos.h"
#include "fft1.h"

#define X2K(a,b,p) b.n=a.n+p; \
                   b.d=2.0*SF_PI/(b.n*a.d); \
                   b.o=(1==b.n)?0:-SF_PI/a.d;
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);

#define KOOP(a) for(ih=0;ih<ah.n;ih++){ \
                for(im=0;im<bm.n;im++){ {a} }}
#define LOOP(a) for(ih=0;ih<ah.n;ih++){ \
                for(im=0;im<am.n;im++){ {a} }}

/*------------------------------------------------------------*/
static axa am,ah,aw,az;
static axa bm;

static float complex **wk,**uk;
static float          *km;
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

    wk = sf_complexalloc2(bm.n,ah.n);
    uk = sf_complexalloc2(bm.n,ah.n);

    /* precompute wavenumbers */
    km = sf_floatalloc(bm.n);
    
    for (im=0; im<bm.n; im++) {
	jm = KMAP(im,bm.n);
	km[im] = bm.o + jm*bm.d;
	km[im] *= km[im];
    }

    fft1_init(bm.n,ah.n,1); /* "1": FFT on axis 1 */
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
/*    float complex kp;     */
/*    float complex s_, sp; */ /* S, S' */
/*    float complex ss, ds; */
/*    float complex a,b,c; */

    float complex cc;
    float complex w2;
    float         s2, aa, fold;

/*    int jn; */
/*------------------------------------------------------------*/

    w2 = 0.1*aw.d - I*w;
    w2 = w2*w2;
    s2 = sz*sz;

    KOOP( wk[ih][im] = 0.; );

    LOOP( wk[ih][im] = 
	  wx[ih][im]; );

    fft1a1(false,wk);
    /*------------------------------------------------------------*/
    /* begin K-domain computation */

    w_ = 4*w2*s2;

    /* loop over km */
    for(im=0; im<bm.n; im++) {
	k_ = km[im];  /* K  */
	
	/* output (h) loop */
	for(ih=0; ih<ah.n; ih++) {
	    hh = ah.o + ih*ah.d;
	    
	    uk[ih][im] = 0;
	    fold = 0.;

	    /* input (ho) loop */
	    for(jh=0; jh<ah.n; jh++) {
		ho = ah.o + jh*ah.d;
		dh = hh - ho;

		h_ = (dh*dh)/(az.d*az.d);    /* H  */

#ifdef fix_later
		kp = - k_ * ((1+h_)/w_);       /* K' = K (1+H)/W */
		
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

		    jn++;
		} /* Newton iterations for S */
		s_ = sp / ((1+h_)/w_);
		

		cc = csqrtf(s_ + k_);
#endif		

/*		cc = -I*w*sz*sqrtf(1.+h_)*az.d;*/
/*		aa = sqrtf(cabsf(cc))/(1+h_);*/
/*		uk[ih][im] += wk[jh][im] * cexpf(-cc) * aa;*/
/*		fold += aa;*/

		cc = -I*2*w*sz*sqrtf(1.+h_)*az.d;
		aa = sqrtf(cabsf(cc)) / (1+h_);

		uk[ih][im] += wk[jh][im] * cexpf(-cc - I*SF_PI/4.) * aa;
		fold += aa;

	    } /* ho */	
	    
	    if (fold > 0.) uk[ih][im] *= sqrtf(2.*w*sz*az.d)/fold;

	} /* hh */	
    } /* km */    

    /* end K-domain computation */
    /*------------------------------------------------------------*/
    fft1a1(false,uk);

    LOOP( wx[ih][im] = 
	  uk[ih][im]; );
}
