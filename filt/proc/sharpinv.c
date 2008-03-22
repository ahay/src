/* Sharpening inversion */
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

#include <rsf.h>

void csharpinv(sf_coperator oper /* inverted operator */, 
	       float scale       /* extra operator scale */,
	       int ncycle        /* number of iterations */,
	       float perc        /* sharpening percentage */,
	       bool verb         /* verbosity flag */,
	       int nq, int np    /* model and data size */,
	       sf_complex *qq    /* model */, 
	       sf_complex *pp    /* data */)
/*< sharp inversion for complex-valued operators >*/
{
    int i, i1;
    sf_complex *q0, *p0=NULL;
    float qdif0=0., pdif0=0., qdif, pdif, pi;

    sf_sharpen_init(nq,perc);
    q0 = sf_complexalloc(nq);

    for (i1=0; i1 < nq; i1++) {
	q0[i1] = qq[i1];
    }

    if (verb) {
	p0 = sf_complexalloc(np);
	for (i1=0; i1 < np; i1++) {
	    p0[i1] = pp[i1];
	}
    }

    for (i=0; i < ncycle; i++) {
	oper(false,false,nq,np,qq,pp);
	for (i1=0; i1 < np; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i1] *= (-scale);
#else
	    pp[i1] = sf_crmul(pp[i1],-scale);
#endif
	}
	oper(true,true,nq,np,qq,pp);

	for (i1=0; i1 < nq; i1++) {
	    qq[i1] += q0[i1];
	}
	sf_csharpen(qq);
	sf_cweight_apply(nq,qq);

	if (verb) {		  
	    pdif = 0.;
	    for (i1=0; i1 < np; i1++) {
		pi = cabsf(pp[i1]-p0[i1]);
		pdif += pi*pi;
	    }
	    qdif = 0.;
	    for (i1=0; i1 < nq; i1++) {
		qdif += cabsf(qq[i1]);
	    }

	    if (0==i) {
		pdif0 = pdif;
		qdif0 = qdif;
		pdif=1.;
		qdif=1.;
	    } else {
		pdif /= pdif0;
		qdif /= qdif0;
	    }

	    sf_warning("iteration %d dres: %f mnorm: %f",i,pdif,qdif);
	}
    }

    free(q0);
    if (NULL != p0) free(p0);
}


