/* Seislet transform using offset continuation */
/*
  Copyright (C) 2009 University of Texas at Austin
   
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

#include "oclet.h"
#include "ocpredict.h"

static int h, nx, nh;
static float dh, h0, w, dw;
static bool inv, unit;
static sf_complex **t, *t1, *t2, *wei;
static void (*transform)(bool);

static void ocpredict_forw(bool adj, sf_complex *tt, int i, int j)
/* Predict forward */
{
    int i2;
    float h, h2;
    for (i2=i; i2 < i+j; i2++) {
	h = h0 + i2*dh;
	h2=h+dh;
	ocpredict_step(adj,true,dw,nx,w,h,h2,tt);
    }
}

static void ocpredict_back(bool adj, sf_complex *tt, int i, int j)    
/* Predict backward */
{
    int i2;
    float h, h2;
    for (i2=i+j; i2 > i; i2--) {
	h2 = h0 + i2*dh;
	h=h2-dh;
	ocpredict_step(adj,false,dw,nx,w,h,h2,tt);
    }
}

static void ochaar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j, i1;

    if (adj) {
	for (j=1; j <= h/2; j *= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		if (inv) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(false,t1,i,j); /* *z0 */
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] -= t1[i1]; 
#else
			t[i+j][i1] = sf_csub(t[i+j][i1],t1[i1]); 
#endif    
			/* d = o - P[e] */
		    }
		    for (i1=0; i1 < nx; i1++) {
			t2[i1] = t[i+j][i1];
		    }
		    ocpredict_back(false,t2,i,j); /* 1/z0 */
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i][i1] += t2[i1]/2; 
#else
			t[i][i1] = sf_cadd(t[i][i1],sf_crmul(t2[i1],0.5)); 
#endif    
                        /* s = e + U[d] */
		    }
		} else {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];  
		    }
		    ocpredict_forw(true,t1,i,j); 
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i][i1] += t1[i1];
#else
			t[i][i1] = sf_cadd(t[i][i1],t1[i1]); 
#endif 
			/* s = e + P'[d] */
		    }
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t2[i1] = -t[i][i1]/2; 
#else
			t2[i1] = sf_crmul(t[i][i1],-0.5); 
#endif 
		    }
		    ocpredict_back(true,t2,i,j); 
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] += t2[i1];
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],t2[i1]); 
#endif 
			/* d = o - U'[s] */
		    }
		}
	    }
	}
    } else {
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
		    t2[i1] = t[i+j][i1];
		}
		ocpredict_back(false,t2,i,j); 
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i][i1] -= t2[i1]/2; 
#else
		    t[i][i1] = sf_csub(t[i][i1],sf_crmul(t2[i1],0.5)); 
#endif 
		    /* e = s - U[d] */
		}
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i][i1];
		}
		ocpredict_forw(false,t1,i,j); 
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i+j][i1] += t1[i1];
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],t1[i1]); 
#endif      
		    /* o = d + P[e] */
		}
	    }
	}
    } 
}

static void oclinear(bool adj)
/* Lifting Linear transform in place */
{
    int i, j, i1;
    
    if (adj) {
	for (j=1; j <= h/2; j *= 2) {
	    if (inv) {
		for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
			t2[i1] = t[i+2*j][i1];
		    }
		    ocpredict_forw(false,t1,i,j);
		    ocpredict_back(false,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] -= (t1[i1]+t2[i1])/2; 
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),-0.5)); 
#endif
                        /* d = o - P[e] */
		    }
		}		
		if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(false,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] -= t1[i1];
#else
			t[i+j][i1] = sf_csub(t[i+j][i1],t1[i1]); 
#endif
			/* d = o - P[e] */
		    }		    
		}
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(false,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[0][i1] += t1[i1]/2;
#else
		    t[0][i1] = sf_cadd(t[0][i1],
					 sf_crmul(t1[i1],0.5)); 
#endif
		}
		for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(false,t1,i,j);
		    ocpredict_forw(false,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i][i1] += (t1[i1]+t2[i1])/4;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),0.25)); 
#endif
			/* s = e + U d */
		    }
		}
	    } else {
		for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
			t2[i1] = t[i+2*j][i1];
		    }
		    ocpredict_forw(true,t1,i,j);
		    ocpredict_back(true,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] -= (t1[i1]+t2[i1])/2; 
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),-0.5)); 
#endif
                        /* d = o - P[e] */
		    }
		}		
		if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(true,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i+j][i1] -= t1[i1];
#else
			t[i+j][i1] = sf_csub(t[i+j][i1],t1[i1]); 
#endif
			/* d = o - P[e] */
		    }		    
		}
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(true,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[0][i1] += t1[i1]/2;
#else
		    t[0][i1] = sf_cadd(t[0][i1],
					 sf_crmul(t1[i1],0.5)); 
#endif
		}
		for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(true,t1,i,j);
		    ocpredict_forw(true,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
			t[i][i1] += (t1[i1]+t2[i1])/4;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),0.25)); 
#endif
			/* s = e + U d */
		    }
		}
	    }
	}
    } else {
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=2*j; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i+j][i1];
		    t2[i1] = t[i-j][i1];
		}
		ocpredict_back(false,t1,i,j);
		ocpredict_forw(false,t2,i-j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i][i1] -= (t1[i1]+t2[i1])/4;
#else
		    t[i][i1] = sf_cadd(t[i][i1],
				       sf_crmul(sf_cadd(t1[i1],t2[i1]),-0.25)); 
#endif
		    /* e = s - U d */
		}
	    }
	    for (i1=0; i1 < nx; i1++) {
		t1[i1] = t[j][i1];
	    }
	    ocpredict_back(false,t1,0,j);
	    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		t[0][i1] -= t1[i1]/2;
#else
		t[0][i1] = sf_cadd(t[0][i1],
				   sf_crmul(t1[i1],-0.5)); 
#endif
	    }
	    for (i=0; i < h-2*j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i][i1];
		    t2[i1] = t[i+2*j][i1];
		}
		ocpredict_forw(false,t1,i,j);
		ocpredict_back(false,t2,i+j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i+j][i1] += (t1[i1]+t2[i1])/2; 
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],
					 sf_crmul(sf_cadd(t1[i1],t2[i1]),0.5)); 
#endif
		    /* o = d + P[e] */
		}
	    }	 
	    if (i+j < h) {
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i][i1];
		}
		ocpredict_forw(false,t1,i,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i+j][i1] += t1[i1];
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],t1[i1]); 
#endif
		    /* o = d + P[e] */
		}		    
	    }
	}
    }
}

static void ocbiorthogonal(bool adj)
/* Lifting Biorthogonal transform in place */
{
    int i, j, i1;
    float a;    
    if (adj) {
	for (j=1; j <= h/2; j *= 2) {
	    if (inv) {
		a = -1.586134342f;
	    	for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
                        t1[i1] = t[i][i1];
                        t2[i1] = t[i+2*j][i1];
                    }
		    ocpredict_forw(false,t1,i,j);
		    ocpredict_back(false,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Predict 1 */
                    }
	        }	 
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(false,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
                        t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/ 
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif 
		    }		    
		}
                a= -0.05298011854f;
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(false,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
	            t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		    t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif 
                }

	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(false,t1,i,j);
		    ocpredict_forw(false,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					   sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Update 1 */
                    }
	        }
                /* Step 1 */
		a = 0.8829110762f;
	    	for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
                        t1[i1] = t[i][i1];
                        t2[i1] = t[i+2*j][i1];
                    }
		    ocpredict_forw(false,t1,i,j);
		    ocpredict_back(false,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Predict 2 */
                    }
	        }	 
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(false,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
                        t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/  
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif 
		    }		    
		}
                a= 0.4435068522f;
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(false,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
	            t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		    t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif 
                }

	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(false,t1,i,j);
		    ocpredict_forw(false,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					   sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Update 2 */
                    }
	        }
                /* Step 2 */
                a= 1/(1.230174105f);
	        for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i+j][i1] *= a;
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif
                    }
	        }
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
                        t[i+j][i1] *= a;  /*right boundary*/ 
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif 
		    }		    
		}
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
	            t[0][i1] /= a;        /*left boundary*/
#else
		    t[0][i1] = sf_crmul(t[0][i1],1./a); 
#endif 

                }
	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i+j][i1] /= a;
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],1./a); 
#endif 
                    }
	        }
		       /* Scale */
	    } else {

		a = -1.586134342f;
	    	for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
                        t1[i1] = t[i][i1];
                        t2[i1] = t[i+2*j][i1];
                    }
		    ocpredict_forw(true,t1,i,j);
		    ocpredict_back(true,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H	
		        t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Predict 1 */
                    }
	        }	 
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(true,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/  
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif
		    }		    
		}
                a= -0.05298011854f;
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(true,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	            t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		    t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif
                }

	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(true,t1,i,j);
		    ocpredict_forw(true,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        t[i][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					   sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Update 1 */
                    }
	        }
                /* Step 1 */
		a = 0.8829110762f;
	    	for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
                        t1[i1] = t[i][i1];
                        t2[i1] = t[i+2*j][i1];
                    }
		    ocpredict_forw(true,t1,i,j);
		    ocpredict_back(true,t2,i+j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],
					     sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Predict 2 */
                    }
	        }	 
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i][i1];
		    }
		    ocpredict_forw(true,t1,i,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/  
#else
			t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif
		    }		    
		}
                a= 0.4435068522f;
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[j][i1];
		}
		ocpredict_back(true,t1,0,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	            t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		    t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif
                }

	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
			t1[i1] = t[i+j][i1];
			t2[i1] = t[i-j][i1];
		    }
		    ocpredict_back(true,t1,i,j);
		    ocpredict_forw(true,t2,i-j,j);
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        t[i][i1] += (t1[i1]+t2[i1])*a;
#else
			t[i][i1] = sf_cadd(t[i][i1],
					   sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Update 2 */
                    }
	        }
                     /* Step 2 */
                a= 1/(1.230174105f);
	        for (i=0; i < h-2*j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        t[i+j][i1] *= a;
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif
                    }
	        }
	        if (i+j < h) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                        t[i+j][i1] *= a;  /*right boundary*/ 
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif 
		    }		    
		}
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	            t[0][i1] /= a;        /*left boundary*/
#else
		    t[0][i1] = sf_crmul(t[0][i1],1./a); 
#endif 
                }
	        for (i=2*j; i < h-j; i += 2*j) {
		    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		        t[i+j][i1] /= a;
#else
			t[i+j][i1] = sf_crmul(t[i+j][i1],1./a); 
#endif 
                    }
	        }
  		       /* Scale */
	    }

	}
    } else {
	for (j=h/2; j >= 1; j /= 2) {

            a= 1.230174105f;
	    for (i=2*j; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j][i1] /= a;
#else
		    t[i+j][i1] = sf_crmul(t[i+j][i1],1./a); 
#endif 
                }
	    }
            for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	        t[0][i1] /= a;        /*left boundary*/
#else
		t[0][i1] = sf_crmul(t[0][i1],1./a); 
#endif 
            }
	    for (i=0; i < h-2*j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j][i1] *= a;
#else
		    t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif 
                }
	    }
	    if (i+j < h) {
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                    t[i+j][i1] *= a;  /*right boundary*/ 
#else
		    t[i+j][i1] = sf_crmul(t[i+j][i1],a); 
#endif  
		}		    
            }
		    /* Undo Scale */

            a= -0.4435068522f;
	    for (i=2*j; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
	            t1[i1] = t[i+j][i1];
		    t2[i1] = t[i-j][i1];
		}
		ocpredict_back(false,t1,i,j);
		ocpredict_forw(false,t2,i-j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i][i1] += (t1[i1]+t2[i1])*a;
#else
		    t[i][i1] = sf_cadd(t[i][i1],
				       sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Undo Update 2 */
                }
	    }
            for (i1=0; i1 < nx; i1++) {
		t1[i1] = t[j][i1];
	    }
	    ocpredict_back(false,t1,0,j);
	    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	        t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif
            }

	    a = -0.8829110762f;
	    for (i=0; i < h-2*j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
                    t1[i1] = t[i][i1];
                    t2[i1] = t[i+2*j][i1];
                }
		ocpredict_forw(false,t1,i,j);
		ocpredict_back(false,t2,i+j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],
					 sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Undo Predict 2 */
                }
	    }	 
	    if (i+j < h) {
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i][i1];
		}
		ocpredict_forw(false,t1,i,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                    t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/  
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif  
		}		    
	    }
                   /* Undo Step 2 */

            a= 0.05298011854f;
	    for (i=2*j; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
	            t1[i1] = t[i+j][i1];
		    t2[i1] = t[i-j][i1];
		}
		ocpredict_back(false,t1,i,j);
		ocpredict_forw(false,t2,i-j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i][i1] += (t1[i1]+t2[i1])*a;
#else
		    t[i][i1] = sf_cadd(t[i][i1],
				       sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		        /* Undo Update 1 */
                }
	    }
            for (i1=0; i1 < nx; i1++) {
		t1[i1] = t[j][i1];
	    }
	    ocpredict_back(false,t1,0,j);
	    for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
	        t[0][i1] += 2*a*t1[i1];      /*left boundary*/
#else
		t[0][i1] = sf_cadd(t[0][i1],sf_crmul(t1[i1],2*a)); 
#endif
            }
	    a = 1.586134342f;
	    for (i=0; i < h-2*j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
                    t1[i1] = t[i][i1];
                    t2[i1] = t[i+2*j][i1];
                }
		ocpredict_forw(false,t1,i,j);
		ocpredict_back(false,t2,i+j,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j][i1] += (t1[i1]+t2[i1])*a;
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],
					 sf_crmul(sf_cadd(t1[i1],t2[i1]),a)); 
#endif
		         /* Undo Predict 1 */
                }
	    }	 
	    if (i+j < h) {
		for (i1=0; i1 < nx; i1++) {
		    t1[i1] = t[i][i1];
		}
		ocpredict_forw(false,t1,i,j);
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
                    t[i+j][i1] += 2*a*t1[i1];  /*right boundary*/  
#else
		    t[i+j][i1] = sf_cadd(t[i+j][i1],sf_crmul(t1[i1],2*a)); 
#endif
		}		    
	    }
                   /* Undo Step 1 */
	}
    }
}


void oclet_init(int nx_in      /* midpoint number   */, 
		int nh_in      /* offset number     */, 
		float dx_in    /* midpoint interval */,
		float dh_in    /* offset interval   */,
		float dw_in    /* freqency interval */,
		float h0_in    /* minumum offset    */,
		bool inv1      /* inversion flag    */, 
		bool unit1     /* weighting flag    */,
		char type      /* transform type    */) 
/*< allocate space >*/
{
    int i,j;
    float wi;

    inv = inv1;
    unit = unit1;
    nx = nx_in;
    nh = nh_in;
    dh = dh_in;
    h0 = h0_in;
    dw = dw_in;

    for (h=1; h < nh; h *= 2) ;
    t = sf_complexalloc2(nx,h);

    t1 = sf_complexalloc(nx);
    t2 = sf_complexalloc(nx);

    switch(type) {
	case 'h': 
	    transform = ochaar;
	    break;
	case 'l':
	    transform = oclinear;
	    break;
	case 'b':
	    transform = ocbiorthogonal;
	    break;
	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }

    if (unit) {
	wei = sf_complexalloc(h);

	wei[0] = sf_cmplx(sqrtf((float) h),0.);
	wi = 0.5;	
	for (j=1; j <= h/2; j *= 2, wi *= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		wei[i+j] = sf_cmplx(sqrtf(wi),0.);
	    }
	}
    }
}

void oclet_close(void) 
/*< deallocate space >*/
{
    free (*t);
    free (t);
    free (t1);
    free (t2);
    if (unit) free(wei);
}

void oclet_lop(bool adj, bool add, int nx1, int ny1, sf_complex *x, sf_complex *y, float w_in)
/*< linear operator >*/
{
    int it, i, j, i1;
    w = w_in;

    sf_cadjnull (adj,add,nx1,ny1,x,y);

    if (adj) {
	for (it=0; it < nx1; it++) {
	    t[0][it] = y[it];
	}
	for (it=nx1; it < nx*h; it++) {
	    t[0][it] = 0.;
	}
    } else {
	for (i1=0; i1 < nx; i1++) {
	    t[0][i1] = x[i1];
	}
	it = nx;
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
		    if (it < ny1) {
			t[i+j][i1]=x[it];
			it++;
		    } else {
			t[i+j][i1]=0.;
		    }
		}
	    }	    	    
	}

	if (unit) {
	    for (it=0; it < h; it++) {
		for (i1=0; i1 < nx; i1++) {
		    if (inv) {
#ifdef SF_HAS_COMPLEX_H
			t[it][i1] /= wei[it];
#else
			t[it][i1] = sf_cdiv(t[it][i1],wei[it]); 
#endif
		    } else {
#ifdef SF_HAS_COMPLEX_H
			t[it][i1] *= wei[it];
#else
			t[it][i1] = sf_cmul(t[it][i1],wei[it]); 
#endif
		    }
		}
	    }
	}
    }

    transform(adj);

    if (adj) {
	if (unit) {
	    for (it=0; it < h; it++) {
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    t[it][i1] *= wei[it];
#else
		    t[it][i1] = sf_cmul(t[it][i1],wei[it]); 
#endif
		}
	    }
	}

	for (i1=0; i1 < nx; i1++) {
	    x[i1] += t[0][i1];
	}
	it = nx;
	for (j=h/2; j >= 1; j /= 2) {
	    for (i=0; i < h-j; i += 2*j) {
		for (i1=0; i1 < nx; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    x[it] += t[i+j][i1];
#else
		    x[it] = sf_cadd(x[it],t[i+j][it]); 
#endif
		    it++;
		    if (it >= ny1) return;
		}
	    }	    	    
	}
    } else {
	for (it=0; it < nx1; it++) {
#ifdef SF_HAS_COMPLEX_H
	    y[it] += t[0][it];
#else
	    y[it] = sf_cadd(y[it],t[0][it]); 
#endif
	}
    }
}
