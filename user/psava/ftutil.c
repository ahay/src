/* 3-D FFT encapsulated */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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

#include <rsf.h>
/*^*/

#include "ftutil.h"

#ifndef _ftutil_h

typedef struct fft *fft3d;
/*^*/

struct fft{
    int n1;
    int n2;
    int n3;
    kiss_fft_cfg forw; /* forward transform */
    kiss_fft_cfg invs; /* inverse transform */
    float scale;
    kiss_fft_cpx *trace;
};
/*^*/

typedef struct sft *sft3d;
/*^*/

struct sft{
    float o;
    float d;
    sf_complex *www;
};
/*^*/

#endif

/*------------------------------------------------------------*/
fft3d fft3a1_init(int n1_, int n2_, int n3_)
/*< initialize FFT on axis 1 >*/
{
    fft3d fft;
    fft = (fft3d) sf_alloc(1,sizeof(*fft));

    fft->n1 = n1_;
    fft->n2 = n2_;
    fft->n3 = n3_;

    fft->forw = kiss_fft_alloc(fft->n1,0,NULL,NULL);
    fft->invs = kiss_fft_alloc(fft->n1,1,NULL,NULL);

    if (NULL == fft->forw || NULL == fft->invs)
	sf_error("%s: KISS FFT allocation error",__FILE__);

    fft->scale = 1./sqrtf(fft->n1);

    return fft;
}

/*------------------------------------------------------------*/
fft3d fft3a2_init(int n1_, int n2_, int n3_)
/*< initialize FFT on axis 2 >*/
{
    fft3d fft;
    fft = (fft3d) sf_alloc(1,sizeof(*fft));

    fft->n1 = n1_; 
    fft->n2 = n2_;
    fft->n3 = n3_;

    fft->forw = kiss_fft_alloc(fft->n2,0,NULL,NULL);
    fft->invs = kiss_fft_alloc(fft->n2,1,NULL,NULL);

    if (NULL == fft->forw || NULL == fft->invs)
	sf_error("%s: KISS FFT allocation error",__FILE__);

    fft->trace = (kiss_fft_cpx*) sf_complexalloc(fft->n2);

    fft->scale = 1./sqrtf(fft->n2);
    
    return fft;
}

/*------------------------------------------------------------*/
fft3d fft3a3_init(int n1_, int n2_, int n3_)
/*< initialize FFT on axis 3 >*/
{    
    fft3d fft;
    fft = (fft3d) sf_alloc(1,sizeof(*fft));

    fft->n1 = n1_; 
    fft->n2 = n2_;
    fft->n3 = n3_;

    fft->forw = kiss_fft_alloc(fft->n3,0,NULL,NULL);
    fft->invs = kiss_fft_alloc(fft->n3,1,NULL,NULL);
    
    if (NULL == fft->forw || NULL == fft->invs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);
    
    fft->trace = (kiss_fft_cpx*) sf_complexalloc(fft->n3);

    fft->scale = 1./sqrtf(fft->n3);

    return fft;
}

/*------------------------------------------------------------*/
void fft3a1_close(fft3d fft)
/*< free allocated storage for FFT on axis 1 >*/
{
    free (fft->forw);
    free (fft->invs);
}

/*------------------------------------------------------------*/
void fft3a2_close(fft3d fft)
/*< free allocated storage for FFT on axis 2 >*/
{
    free (fft->trace);

    free (fft->forw);
    free (fft->invs);
}

/*------------------------------------------------------------*/
void fft3a3_close(fft3d fft)
/*< free allocated storage for FFT on axis 3 >*/
{
    free (fft->trace);

    free (fft->forw);
    free (fft->invs);
}

/*------------------------------------------------------------*/
void fft3a1(bool inv           /* inverse/forward flag */, 
	    kiss_fft_cpx ***pp /* [n1][n2][n3] */,
	    fft3d fft) 
/*< apply FFT on axis 1 >*/
{
    int i1, i2, i3;
    
    if (inv) {

	/* IFT 1 */
	for    (i3=0; i3 < fft->n3; i3++) {
	    for(i2=0; i2 < fft->n2; i2++) {
		kiss_fft(fft->invs,pp[i3][i2],pp[i3][i2]);
	    }
	}

	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}
	
    } else {
	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}

	/* FFT 1 */
	for    (i3=0; i3 < fft->n3; i3++) {
	    for(i2=0; i2 < fft->n2; i2++) {
		kiss_fft(fft->forw,pp[i3][i2],pp[i3][i2]);
	    }
	}

    }
}

/*------------------------------------------------------------*/
void fft3a2(bool inv           /* inverse/forward flag */, 
	    kiss_fft_cpx ***pp /* [n1][n2][n3] */,
	    fft3d fft) 
/*< apply FFT on axis 2 >*/
{
    int i1, i2, i3;
    
    if (inv) {

	/* IFT 2 */
	for    (i3=0; i3 < fft->n3; i3++) {
	    for(i1=0; i1 < fft->n1; i1++) {
		kiss_fft_stride(fft->invs,pp[i3][0]+i1,fft->trace,fft->n1);
		for(i2=0; i2 < fft->n2; i2++) {
		    pp[i3][i2][i1] = fft->trace[i2];
		}
	    }
	}
	
	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}
	
    } else {
	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}

	/* FFT 2 */
	for    (i3=0; i3 < fft->n3; i3++) {
	    for(i1=0; i1 < fft->n1; i1++) {
		kiss_fft_stride(fft->forw,pp[i3][0]+i1,fft->trace,fft->n1);
		for(i2=0; i2 < fft->n2; i2++) {
		    pp[i3][i2][i1] = fft->trace[i2];
		}
	    }
	}

    }
}

/*------------------------------------------------------------*/
void fft3a3(bool inv           /* inverse/forward flag */, 
	    kiss_fft_cpx ***pp /* [n1][n2][n3] */,
	    fft3d fft) 
/*< apply FFT on axis 3 >*/
{
    int i1, i2, i3;
    
    if (inv) {

	/* IFT 3 */
	for    (i2=0; i2 < fft->n2; i2++) {
	    for(i1=0; i1 < fft->n1; i1++) {
		kiss_fft_stride(fft->invs,pp[0][0]+i1+i2*fft->n1,fft->trace,fft->n1*fft->n2);
		for(i3=0; i3 < fft->n3; i3++) {
		    pp[i3][i2][i1] = fft->trace[i3];
		}
	    }
	}

	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}
	
    } else {
	/* scaling */
	for        (i3=0; i3 < fft->n3; i3++) {
	    for    (i2=0; i2 < fft->n2; i2++) {
		for(i1=0; i1 < fft->n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->scale);
		}
	    }
	}

	/* FFT 3 */
	for    (i2=0; i2 < fft->n2; i2++) {
	    for(i1=0; i1 < fft->n1; i1++) {
		kiss_fft_stride(fft->forw,pp[0][0]+i1+i2*fft->n1,fft->trace,fft->n1*fft->n2);
		for(i3=0; i3 < fft->n3; i3++) {
		    pp[i3][i2][i1] = fft->trace[i3];
		}
	    }
	}

    }
}

/*------------------------------------------------------------*/
void cnt3a1(sf_complex ***pp,
	    fft3d fft)
/*< apply centering on axis 1 >*/
{
    int i1,i2,i3;

    for        (i3=0; i3<fft->n3; i3++) {
	for    (i2=0; i2<fft->n2; i2++) {
	    for(i1=1; i1<fft->n1; i1+=2){
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] = - pp[i3][i2][i1];
#else
		pp[i3][i2][i1] = sf_cneg(pp[i3][i2][i1]);
#endif
	    }
	}
	
    }
}

/*------------------------------------------------------------*/
void cnt3a2(sf_complex ***pp,
	    fft3d fft)
/*< apply centering on axis 2 >*/
{
    int i1,i2,i3;

    for        (i3=0; i3<fft->n3; i3++) {
	for    (i2=1; i2<fft->n2; i2+=2){
	    for(i1=0; i1<fft->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] = - pp[i3][i2][i1];
#else
		pp[i3][i2][i1] = sf_cneg(pp[i3][i2][i1]);
#endif
	    }
	}
    }
}

/*------------------------------------------------------------*/
void cnt3a3(sf_complex ***pp,
	    fft3d fft)
/*< apply centering on axis 3>*/
{
    int i1,i2,i3;

    for        (i3=0; i3<fft->n3; i3+=2){
	for    (i2=1; i2<fft->n2; i2++) {
	    for(i1=0; i1<fft->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] = - pp[i3][i2][i1];
#else
		pp[i3][i2][i1] = sf_cneg(pp[i3][i2][i1]);
#endif
	    }
	}
    }
}

/*------------------------------------------------------------*/
sft3d sft3_init(int n,
		float o,
		float d)
/*< origin shift (assumes no centering) >*/
{
    int k,i;
    float w,s;
    
    sft3d sft;
    sft = (sft3d) sf_alloc(1,sizeof(*sft));
    
    w=2.0*SF_PI/(n*d) * o;
    k=n/2;
    
    sft->www = sf_complexalloc(n);
    for(i=0; i<n; i++) { sft->www[i]=sf_cmplx(1.0,0.0); }
    
    for(i=0; i<k; i++) {
	s = w * i;
	sft->www[i]   = sf_cmplx(cosf(s),sinf(s));
	
	s = w * (-k-1+i);
	s = w * (-k+i);
	sft->www[k+i] = sf_cmplx(cosf(s),sinf(s));
    }
    
    return sft;
}

/*------------------------------------------------------------*/
void sft3_close(sft3d sft)
/*< close shift >*/
{
    free(sft->www);
}

/*------------------------------------------------------------*/
void sft3a3(sf_complex ***pp,
	    sft3d sft,
	    fft3d fft)
/*< apply shift on axis 3 >*/
{
    int i1,i2,i3;

    for        (i3=0; i3<fft->n3; i3++) {
	for    (i2=0; i2<fft->n2; i2++) {
	    for(i1=0; i1<fft->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] *= sft->www[i3];
#else
		pp[i3][i2][i1] = sf_cmul(pp[i3][i2][i1],sft->www[i3]);
#endif
	    }
	}    
    }
}

/*------------------------------------------------------------*/
void sft3a2(sf_complex ***pp,
	    sft3d sft,
	    fft3d fft)
/*< apply shift on axis 2 >*/
{
    int i1,i2,i3;

    for        (i3=0; i3<fft->n3; i3++) {
	for    (i2=0; i2<fft->n2; i2++) {
	    for(i1=0; i1<fft->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] *= sft->www[i2];
#else
		pp[i3][i2][i1] = sf_cmul(pp[i3][i2][i1],sft->www[i2]);
#endif
	    }
	}    
    }
}

/*------------------------------------------------------------*/
void sft3a1(sf_complex ***pp,
	    sft3d sft,
	    fft3d fft)
/*< apply shift on axis 1 >*/
{
    int i1,i2,i3;
    
    for        (i3=0; i3<fft->n3; i3++) {
	for    (i2=0; i2<fft->n2; i2++) {
	    for(i1=0; i1<fft->n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] *= sft->www[i1];
#else
		pp[i3][i2][i1] = sf_cmul(pp[i3][i2][i1],sft->www[i1]);
#endif
	    }
	}    
    }
}
    
