#include <math.h>

#include <rsf.h>

#include "gazdagz.h"

static float eps, dt, dz, *vt, dw, fw;
static int nt, nz, nw;
static float complex *pp;
static bool depth, midpoint;

void gazdagz_init (float eps1, int nt1, float dt1, 
		   int nz1, float dz1, float *vt1, bool depth1, bool midpoint)
{
    eps = eps1; 
    nt = nt1; dt = dt1; 
    nz = nz1; dz = dz1;
    vt = vt1;
    depth = depth1;

    /* determine frequency sampling */
    nw = sf_npfa(nt);
    dw = 2.0*SF_PI/(nw*dt);
    fw = -SF_PI/dt;
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);
}

void gazdagz_close ()
{    
    /* free workspace */
    free (pp);	
}

void gazdagz (bool inv, float k2, float complex *p, float complex *q)
{
    int it,iz,iw;
    float complex cshift1, cshift2, cshift, w2, y;
    float x;

    if (inv) { /* modeling */
	for (iw=0; iw<nw; iw++) {
	    pp[iw] = q[nz-1];
	}

	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w2 = eps*dw + I*(fw + iw*dw);
		
		if (depth) {
		    if (midpoint) {
			cshift = csqrtf(0.5 * w2*w2 * (vt[iz]+vt[iz+1]) + k2);
		    } else {
			cshift1 = csqrtf(w2*w2 * vt[iz] + k2);
			cshift2 = csqrtf(w2*w2 * vt[iz+1] + k2);
			
			cshift = (cshift1 + cshift2 - 
				  1./(1./cshift1+1./cshift2))/1.5;
		    }
		} else {
		    if (midpoint) {
			cshift = csqrtf(w2*w2 + 0.5 * (vt[iz]+vt[iz+1]) * k2);
		    } else {
			cshift1 = csqrtf(w2*w2 + vt[iz] * k2);
			cshift2 = csqrtf(w2*w2 + vt[iz+1] * k2);
			
			x = 1./(1.+vt[iz+1]/vt[iz]);
			cshift = cshift1 + x*(cshift2-cshift1);
			y = x*cshift2/cshift;
			cshift *= (1.-y*(1.-y))/(1.-x*(1.-x));
		    }
		}

		/* cshift = cexpf(-cshift*dz)*csqrt(cshift2/cshift1); */
		cshift = cexpf(-cshift*dz);
		pp[iw] = pp[iw]*cshift + q[iz];
	    }
	}
	
	sf_pfacc(1,nw,pp);
	for (it=0; it<nt; it++) {
	    p[it] = (it%2)? -pp[it] : pp[it];
	}
    } else { /* migration */
	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it++) {
	    /* scale accumulated image just as we would for an FFT */
	    pp[it] = (it%2 ? -p[it] : p[it])/nw;
	}
	for (it=nt; it<nw; it++) {
	    pp[it] = 0.0;
	}
	sf_pfacc(-1,nw,pp);
    
	/* loop over migrated times z */
	for (iz=0; iz<nz-1; iz++) {
	    /* initialize migrated sample */
	    q[iz] = 0.0;
      
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		/* accumulate image (summed over frequency) */
		q[iz] += pp[iw];

		w2 = eps*dw + I*(fw + iw*dw);
	
		if (depth) {
		    /* extrapolate down one depth step */
		    if (midpoint) {
			cshift = csqrtf(0.5 * w2*w2 * (vt[iz]+vt[iz+1]) + k2);
		    } else {
			cshift1 = csqrtf(w2*w2 * vt[iz] + k2);
			cshift2 = csqrtf(w2*w2 * vt[iz+1] + k2);
			
			cshift = (cshift1 + cshift2 - 
				  1./(1./cshift1+1./cshift2))/1.5;
		    }
		} else {
		    /* extrapolate down one migrated time step */
		    if (midpoint) {
			cshift = csqrtf(w2*w2 + 0.5 * (vt[iz]+vt[iz+1]) * k2);
		    } else {
			cshift1 = csqrtf(w2*w2 + vt[iz] * k2);
			cshift2 = csqrtf(w2*w2 + vt[iz+1] * k2);
			
			x = 1./(1.+vt[iz+1]/vt[iz]);
			cshift = cshift1 + x*(cshift2-cshift1);
			y = x*cshift2/cshift;
			cshift *= (1.-y*(1.-y))/(1.-x*(1.-x));
		    }
		}
		    
		/* cshift = conjf(cexpf(-cshift*dz)*csqrt(cshift2/cshift1)); */
		cshift = conjf(cexpf(-cshift*dz));
		pp[iw] *= cshift;
	    }
	}
	
	q[nz-1] = 0.0;
	/* loop over frequencies w */
	for (iw=0; iw<nw; iw++) {
	    /* accumulate image (summed over frequency) */
	    q[nz-1] += pp[iw];
	}
    }
}

/* 	$Id$	 */
