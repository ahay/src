/* K-domain phase shift */

#include <math.h>
#include <rsf.h>
#include "fft2.h"
#include "phs.h"

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

static float complex **wk;

static float          **ks;
static float          **kr;
/*------------------------------------------------------------*/

void phs_init(axa am_,
	      axa ah_,
	      axa aw_,
	      axa az_)
/*< initialize >*/
{
    int   im,ih;
    int   jm,jh;
    float km,kh,k;

    am = am_;
    ah = ah_;
    aw = aw_;
    az = az_;

    X2K(am,bm,0);
    X2K(ah,bh,0);
    fft2_init(bm.n,bh.n);

    wk = sf_complexalloc2(bm.n,bh.n);

    /* precompute wavenumbers */
    ks= sf_floatalloc2(bm.n,bh.n);
    kr= sf_floatalloc2(bm.n,bh.n);
    
    for (im=0; im<bm.n; im++) {
	jm = KMAP(im,bm.n);
	km = bm.o + jm*bm.d;
	
	for (ih=0; ih<bh.n; ih++) {
	    jh = KMAP(ih,bh.n);
	    kh = bh.o + jh*bh.d;
	    
	    k = 0.5*(km-kh);
	    ks[ih][im] = k*k;  /* ks^2 */
	    
	    k = 0.5*(km+kh);
	    kr[ih][im] = k*k;  /* kr^2 */
	}
    }
} 

void phs(float            w,
	 float           sz,
	 float complex **wx)
/*< extrapolate >*/
{
    int   ih,im;

    float complex cs,cr,cc;
    float complex w2;
    float         s2;
    w2 = 0.1*aw.d - I*w;
    w2 = w2*w2;
    s2 = sz*sz;

    KOOP( wk[ih][im] = 0.; );

    LOOP( wk[ih][im] = 
	  wx[ih][im]; );
    fft2(false,wk);
    
/*------------------------------------------------------------*/
    KOOP(
	cs  = csqrtf(w2*s2 + ks[ih][im]);
	cr  = csqrtf(w2*s2 + kr[ih][im]);

	cc = cs + cr;
	wk[ih][im] *= cexpf(-cc*az.d);
	);
/*------------------------------------------------------------*/
    fft2(true,wk);
    LOOP( wx[ih][im] = 
	  wk[ih][im]; );
}
