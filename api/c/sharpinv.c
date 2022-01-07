/* Sharpening inversion added Bregman iteration */
/*
  Copyright (C) 2008 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/* m_k+1 = S[ A*d + (I-A*F)*m_k], m_0 = A*d */

#include "_bool.h"
#include "_solver.h"
#include "c99.h"
#include "komplex.h"
#include "alloc.h"
#include "error.h"
#include "sharpen.h"
#include "weight.h"

void sf_csharpinv(sf_coperator oper /* inverted operator */, 
		  float scale       /* extra operator scaling */,
		  int niter         /* number of outer iterations */,
		  int ncycle        /* number of iterations */,
		  float perc        /* sharpening percentage */,
		  bool verb         /* verbosity flag */,
		  int nq, int np    /* model and data size */,
		  sf_complex *qq    /* model */, 
		  sf_complex *pp    /* data */,
		  bool twhole       /* thresholding flag */)
/*< sharp inversion for complex-valued operators >*/
{
    int iter, i, i1, ti;
    sf_complex *q0, *p0, *p1, *tp=NULL;
    float qdif0=0., pdif0=0., qdif, pdif, pi;

    if (!twhole) {
	sf_sharpen_init(np,perc,0.5);
    } else {
	sf_sharpen_init(nq,perc,0.5);
    }
    q0 = sf_complexalloc(nq);
    p0 = sf_complexalloc(np);
    p1 = sf_complexalloc(np);

    for (i1=0; i1 < np; i1++) {
	p0[i1] = pp[i1];
    }

    for (iter=0; iter < niter; iter++) { /* outer iteration */
	oper(true,false,nq,np,qq,p0);
	
	for (i1=0; i1 < nq; i1++) {
	    q0[i1] = qq[i1];
	}
	for (i1=0; i1 < np; i1++) {
	    p1[i1] = p0[i1];
	}

	for (i=0; i < ncycle; i++) { /* inner iteration */
	    oper(false,false,nq,np,qq,p1);	    
	    for (i1=0; i1 < np; i1++) {
#ifdef SF_HAS_COMPLEX_H
		p1[i1] *= (-scale);
#else
		p1[i1] = sf_crmul(pp[i1],-scale);
#endif
	    }
	    oper(true,true,nq,np,qq,p1);

	    for (i1=0; i1 < nq; i1++) {
#ifdef SF_HAS_COMPLEX_H
		qq[i1] += q0[i1];
#else
		qq[i1] = sf_cadd(qq[i1],q0[i1]);
#endif
	    }

	    if (!twhole) {
		tp = sf_complexalloc(np);
		for (ti=0; ti < (nq/np); ti++) {
		    for (i1=0; i1 < np; i1++) {
			tp[i1] = qq[ti*np+i1];
		    }
		    sf_csharpen(tp);
		    sf_cweight_apply(np,tp);
		    for (i1=0; i1 < np; i1++) {
			qq[ti*np+i1] = tp[i1];
		    }
		}
	    } else {
		sf_csharpen(qq);
		sf_cweight_apply(nq,qq);
	    }

	    if (verb) {		  	    
		qdif = 0.;
		for (i1=0; i1 < nq; i1++) {
		    qdif += cabsf(qq[i1]);
		}

		if (0==i) {
		    qdif0 = qdif;
		    qdif=1.;
		} else {
		    qdif /= qdif0;
		}

		sf_warning("inner iteration %d mnorm: %f",i,qdif);
	    }
	} /* inner iteration */
	
	oper(false,false,nq,np,qq,p1);
	for (i1=0; i1 < np; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    p1[i1] *= (-scale);
#else
	    p1[i1] = sf_crmul(pp[i1],-scale);
#endif	
	}

	for (i1=0; i1 < np; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    p0[i1] += pp[i1] + p1[i1];
#else
	    p0[i1] = sf_cadd(p0[i1],sf_cadd(pp[i1],p1[i1]));
#endif
	}

	if (verb) {
	    pdif = 0.;
	    for (i1=0; i1 < np; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pi = cabs(pp[i1]+p1[i1]);
#else
		pi = sf_cabsf(sf_cadd(pp[i1],p1[i1]));
#endif
		pdif += pi*pi;
	    }
	    if (0==iter) {
		pdif0 = pdif;
		pdif=1.;
	    } else {
		pdif /= pdif0;
	    }
	    sf_warning("outer iteration %d dres: %f",iter,pdif);
	}
    } /* outer iteration */

    free(q0);
    free(p0);
    free(p1);
    if (!twhole) free(tp);
    sf_sharpen_close();
}

void sf_sharpinv(sf_operator oper  /* inverted operator */, 
		 float scale       /* extra operator scaling */,
		 int niter         /* number of outer iterations */,
		 int ncycle        /* number of iterations */,
		 float perc        /* sharpening percentage */,
		 bool verb         /* verbosity flag */,
		 int nq, int np    /* model and data size */,
		 float *qq         /* model */, 
		 float *pp         /* data */,
                 bool twhole       /* thresholding flag */)
/*< sharp inversion for real-valued operators >*/
{
    int iter, i, i1, ti;
    float *q0, *p0, *p1, *tp=NULL;
    float qdif0=0., pdif0=0., qdif, pdif, pi;

    if (!twhole) {
	sf_sharpen_init(np,perc,0.5);
    } else {
	sf_sharpen_init(nq,perc,0.5);
    }
    q0 = sf_floatalloc(nq);
    p0 = sf_floatalloc(np);
    p1 = sf_floatalloc(np);

    for (i1=0; i1 < np; i1++) {
	p0[i1] = pp[i1];
    }

    for (iter=0; iter < niter; iter++) { /* outer iteration */
	oper(true,false,nq,np,qq,p0);

	for (i1=0; i1 < nq; i1++) {
	    q0[i1] = qq[i1];
	}
 	for (i1=0; i1 < np; i1++) {
	    p1[i1] = p0[i1];
	}
 
	for (i=0; i < ncycle; i++) { /* inner iteration */
	    oper(false,false,nq,np,qq,p1);
	    for (i1=0; i1 < np; i1++) {
		p1[i1] *= (-scale);
	    }
	    oper(true,true,nq,np,qq,p1);

	    for (i1=0; i1 < nq; i1++) {
		qq[i1] += q0[i1];
	    }
	    if (!twhole) {
		tp = sf_floatalloc(np);
		for (ti=0; ti < (nq/np); ti++) {
		    for (i1=0; i1 < np; i1++) {
			tp[i1] = qq[ti*np+i1];
		    }
		    sf_sharpen(tp);
		    sf_weight_apply(np,tp);
		    for (i1=0; i1 < np; i1++) {
			qq[ti*np+i1] = tp[i1];
		    }
		}
	    } else {
		sf_sharpen(qq);
		sf_weight_apply(nq,qq);
	    }

	    if (verb) {		  
		qdif = 0.;
		for (i1=0; i1 < nq; i1++) {
		    qdif += fabsf(qq[i1]);
		}

		if (0==i) {
		    qdif0 = qdif;
		    qdif=1.;
		} else {
		    qdif /= qdif0;
		}

		sf_warning("inner iteration %d mnorm: %f",i,qdif);
	    }
	} /* inner iteration */

	oper(false,false,nq,np,qq,p1);
	for (i1=0; i1 < np; i1++) {
	    p1[i1] *= (-scale);
	}
	for (i1=0; i1 < np; i1++) {
	    p0[i1] += pp[i1] + p1[i1];
	}

	if (verb) {
	    pdif = 0.;
	    for (i1=0; i1 < np; i1++) {
		pi = fabsf(pp[i1]+p1[i1]);
		pdif += pi*pi;
	    }
	    if (0==iter) {
		pdif0 = pdif;
		pdif=1.;
	    } else {
		pdif /= pdif0;
	    }
	    sf_warning("outer iteration %d dres: %f",iter,pdif);
	}
    } /* outer iteration */


    free(q0);
    free(p0);
    free(p1);
    if(!twhole) free(tp);
    sf_sharpen_close();
}
