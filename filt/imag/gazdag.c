#include <math.h>

#include <rsf.h>

#include "gazdag.h"

static float eps, dt, dz, *vt, dw, fw;
static int nt, nz;
static float complex *pp;
static bool depth;
static kiss_fft_cfg forw, invs;

void gazdag_init (float eps1, int nt1, float dt1, 
		  int nz1, float dz1, float *vt1, bool depth1)
{
    eps = eps1; 
    nt = nt1; dt = dt1; 
    nz = nz1; dz = dz1;
    vt = vt1;
    depth = depth1;

    /* determine frequency sampling */
    dw = 2.0*SF_PI/(nt*dt);
    fw = -SF_PI/dt;
    
    forw = kiss_fft_alloc(nt,0,NULL,NULL);
    invs = kiss_fft_alloc(nt,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);
    
    /* allocate workspace */
    pp = sf_complexalloc (nt);
}

void gazdag_close ()
{    
    /* free workspace */
    free (pp);	
}

void gazdag (bool inv, float k2, float complex *p, float complex *q)
{
    int it,iz,iw;
    float complex cshift, w2;
	
    if (inv) { /* modeling */
	for (iw=0; iw<nt; iw++) {
	    pp[iw] = q[nz-1];
	}

	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nt; iw++) {
		w2 = eps*dw + I*(fw + iw*dw);

		if (depth) {
		    w2 = w2*w2 * vt[iz] + k2;
		} else {
		    w2 = w2*w2 + vt[iz] * k2;
		}
	
		cshift = cexpf(-csqrtf(w2)*dz);
		pp[iw] = pp[iw]*cshift + q[iz];
	    }
	}

	kiss_fft(forw,(const kiss_fft_cpx *) pp, 
		 (kiss_fft_cpx *) pp);
	for (it=0; it<nt; it++) {
	    p[it] = (it%2)? -pp[it] : pp[it];
	}
    } else { /* migration */
	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it++) {
	    /* scale accumulated image just as we would for an FFT */
	    pp[it] = (it%2 ? -p[it] : p[it])/nt;
	}
	kiss_fft(invs,(const kiss_fft_cpx *) pp, 
		 (kiss_fft_cpx *) pp);
    
	/* loop over migrated times z */
	for (iz=0; iz<nz; iz++) {
	    /* initialize migrated sample */
	    q[iz] = 0.0;
      
	    /* loop over frequencies w */
	    for (iw=0; iw<nt; iw++) {
		/* accumulate image (summed over frequency) */
		q[iz] += pp[iw];

		w2 = eps*dw + I*(fw + iw*dw);

		if (depth) { /* depth migration */
		    w2 = w2*w2 * vt[iz] + k2;
		} else { /* time migration */
		    w2 = w2*w2 + vt[iz] * k2;
		}
	
		/* extrapolate down one migrated time step */
		cshift = conjf(cexpf(-csqrtf(w2)*dz));
		pp[iw] *= cshift;
	    }
	}
    }
}

/* 	$Id$	 */
