#include <math.h>

#include <rsf.h>

#include "gazdag.h"

static float eps, dz, *vt, dw;
static int nz, nw;
static float complex *pp;
static bool depth;
static kiss_fftr_cfg forw, invs;

void gazdag_init (float eps1, int nt, float dt, 
                  int nz1, float dz1, float *vt1, bool depth1)
{
    eps = eps1; 
    nz = nz1; dz = dz1;
    vt = vt1;
    depth = depth1;

    /* determine frequency sampling */
    nw = nt/2+1;
    dw = 2.0*SF_PI/(nt*dt);
    
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
    
    /* allocate workspace */
    pp = sf_complexalloc (nw);
}

void gazdag_close ()
{    
    /* free workspace */
    free (pp);  
    free (forw);
    free (invs);
}

void gazdag (bool inv, float k2, float *p, float *q)
{
    int iz,iw;
    float complex cshift, w2;
        
    if (inv) { /* modeling */
        for (iw=0; iw<nw; iw++) {
            pp[iw] = q[nz-1];
        }

        /* loop over migrated times z */
        for (iz=nz-2; iz>=0; iz--) {
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
                w2 = (eps + I*iw)*dw;

                if (depth) {
                    w2 = w2*w2 * vt[iz] + k2;
                } else {
                    w2 = w2*w2 + vt[iz] * k2;
                }
        
                cshift = cexpf(-csqrtf(w2)*dz);
                pp[iw] = pp[iw]*cshift + q[iz];
            }
        }

        kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    } else { /* migration */
	kiss_fftr(forw, p, (kiss_fft_cpx *) pp);
    
        /* loop over migrated times z */
        for (iz=0; iz<nz; iz++) {
            /* initialize migrated sample */
            q[iz] = 0.0;
      
            /* loop over frequencies w */
            for (iw=0; iw<nw; iw++) {
                /* accumulate image (summed over frequency) */
                q[iz] += crealf(pp[iw]);

                w2 = (eps + I*iw)*dw;

                if (depth) { /* depth migration */
                    w2 = w2*w2 * vt[iz] + k2;
                } else { /* time migration */
                    w2 = w2*w2 + vt[iz] * k2;
                }
        
                /* extrapolate down one migrated time step */
                cshift = cexpf(-csqrtf(w2)*dz);
                pp[iw] *= conjf(cshift);
            }
        }
    }
}

/*      $Id$     */
