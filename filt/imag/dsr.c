#include <math.h>

#include <rsf.h>

void dsr (int inv, float eps, float kx, float kh, 
	  int nw, float dw, float fw, 
	  int nz, float dz, 
	  float *vt, float complex *p, float complex *q)
{
    int iz,iw;
    float vs2, vr2, s, r;
    float complex cshift,w2;

    s = 0.5*(kx-kh);
    r = 0.5*(kx+kh);

    if (inv) { /* modeling */
	for (iw=0; iw<nw; iw++) {
	    p[iw] = q[nz-1];
	}
	
	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    vs2 = vt[iz]*s; vs2 *= vs2;
	    vr2 = vt[iz]*r; vr2 *= vr2;

	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w2 = eps*dw + I*(fw + iw*dw);
		w2 *= w2;
		w2 = csqrtf(w2+vs2)+ csqrtf(w2+vr2);
	
		cshift = conjf(cexpf(-0.5*w2*dz));
		p[iw] = p[iw]*cshift+q[iz];
	    }
	}
    } else { /* migration */
	/* loop over migrated times z */
	for (iz=0; iz<nz; iz++) {
	    vs2 = vt[iz]*s; vs2 *= vs2;
	    vr2 = vt[iz]*r; vr2 *= vr2;
      
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		/* accumulate image (summed over frequency) */
		q[iz] += p[iw];

		w2 = eps*dw+I*(fw + iw*dw);
		w2 *= w2;
		w2 = csqrtf(w2+vs2)+ csqrtf(w2+vr2);
	
		/* extrapolate down one migrated time step */
		cshift = cexpf(-0.5*w2*dz);
		p[iw] *= cshift;
	    }
	}
    }
}


