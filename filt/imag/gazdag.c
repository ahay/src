#include <math.h>
#include <float.h>

#include <rsf.h>

#include "gazdag.h"

static float eps, dt, dz, *vt, *gt, dw, fw;
static int nt, nz, nw;
static float complex *pp;

void gazdag_init (float eps1, int nt1, float dt1, 
		  int nz1, float dz1, float *vt1, float *gt1)
{
    nt = nt1; dt = dt1; 
    nz = nz1; dz = dz1;
    vt = vt1; gt = gt1;

    /* determine frequency sampling */
    nw = sf_npfa(nt);
    dw = 2.0*SF_PI/(nw*dt);
    fw = -SF_PI/dt;
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);

    eps = eps1*dw; 
    eps *= eps;
}

void gazdag_close ()
{    
    /* free workspace */
    free (pp);	
}

void gazdag (bool inv, float k2, float complex *p, float complex *q)
{
    int it,iz,iw,ig;
    float complex cshift;
    float w, w2, kz, gz, add, kz2;
	
    if (inv) { /* modeling */
	for (iw=0; iw<nw; iw++) {
	    pp[iw] = q[nz-1];
	}

	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
/*
  w2 = eps*dw + I*(fw + iw*dw);
  w2 = w2*w2 + vt[iz] * k2;
  
  cshift = cexpf(-csqrtf(w2)*dz);
  pp[iw] = pp[iw]*cshift + q[iz];
*/
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
	for (iz=0; iz<nz; iz++) {
	    /* initialize migrated sample */
	    q[iz] = 0.0;
      
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w = fw + iw*dw;

		w2 = w*w;
		kz = w2 - vt[iz] * k2;

		/* accumulate image (summed over frequency) */
		q[iz] += pp[iw];

		if ((nz-1)*(nz-1)*kz <= iz*iz*w2) {
		    pp[iw] = 0.;
		    continue;
		}

		kz /= w2;
	
		if (NULL != gt && 0. != gt[iz]) {
		    gz = gt[iz] * vt[iz] * w2 * dz;		    
		    kz2 = kz;
		    for (ig=0; ig < 5; ig++) {
			add = 0.25 * kz * kz + gz * vt[iz] *
			    (2.*sqrtf(k2 * kz2) + gz);
			if (add < 0.) break;
			kz2 = 0.5 * kz + sqrtf(add);
		    }
		    kz = kz2;
		}

		/* extrapolate down one migrated time step */
		cshift = cexpf(I*sqrtf(kz)*w*dz);	    
		pp[iw] *= cshift;
	    }
	}
    }
}

/* 	$Id: gazdag.c,v 1.6 2003/12/04 16:06:20 fomels Exp $	 */
