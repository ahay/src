/* Generate internal multiples by using virtual seismic data. */
/*
  Copyright (C) 2017 Jilin University
   
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

int main(int argc, char *argv[])
{
    int nw, n1, n2, n4, iw, i1, i2, i4, s1=0, x1=0, s2=0, x2=0, m, fold;
    int jumpo, jumps, newn1, newn2, tn;
    float d1, d2, newd1, newd2;
    bool verb, stack, both;

    sf_complex *dd, *ref=NULL, *mm=NULL, *mtemp=NULL;
    sf_file in, out, dif=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    
    if (d1!=d2) sf_error("Need d1==d2");

    n4 = sf_leftsize(in,3);
    fold = 0;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getbool("stack",&stack)) stack=false;
    /* stack flag, if y, no common multiple gather */

    if (!sf_getbool("both",&both)) both=false;
    /* receiver flag, if y, receiver with both sides */

    if (!stack) {
	if (!sf_getint("jumpo",&jumpo)) jumpo=1;
	/* jump in offset dimension, only for stack=n */
	
	if (!sf_getint("jumps",&jumps)) jumps=1;
	/* jump in shot dimension, only for stack=n  */
    }
    newn1 = n1;
    newn2 = n2;
    if (!stack) {
	if (!both) {
	    sf_putint(out,"n1",(2*n1-1));
	    sf_putfloat(out,"d1",d1);
	    sf_putfloat(out,"o1",(1-n1)*d1);
	} else {
	    sf_putint(out,"n1",n1);
	    sf_putfloat(out,"d1",d1);
	    sf_putfloat(out,"o1",(-1*n1/2)*d1);
	}
	(void) sf_shiftdim(in, out, 1);
	if (n1%jumpo == 0) {
	    newn1 = n1 / jumpo;
	} else {
	    newn1 = n1 / jumpo +1;
	}
	if (n2%jumps == 0) {
	    newn2 = n2 / jumps;
	} else {
	    newn2 = n2 / jumps +1;
	}
	newd1 = (float) (d1 * jumpo);
	newd2 = (float) (d2 * jumps);
	sf_putint(out,"n2",newn1);
	sf_putfloat(out,"d2",newd1);
	sf_putint(out,"n3",newn2);
	sf_putfloat(out,"d3",newd2);
    }

    if (NULL != sf_getstring ("dif")) {
	dif = sf_input("dif");
	ref = sf_complexalloc(n1*n2);
    } else {
	dif = NULL;
    }	
    
    dd = sf_complexalloc(n1*n2);

    if (stack) {
	mm = sf_complexalloc(n1*n2);
    } else {
	if (!both) {
	    mm = sf_complexalloc((2*n1-1)*n1*n2);
	    mtemp = sf_complexalloc((2*n1-1)*newn1*newn2);
	} else {
	    mm = sf_complexalloc(n1*n1*n2);
	    mtemp = sf_complexalloc(n1*newn1*newn2);
	}
    }

    /* loop over n4 */
    for (i4=0; i4 < n4; i4++) {
	if (verb) sf_warning("slice %d of %d",i4+1,n4);
	for (iw=0; iw < nw; iw++) { /* loop over frequency */
	    if (verb) sf_warning("frequency %d of %d;",iw+1,nw);
	    sf_complexread(dd,n1*n2,in);
	    if (NULL != dif) {
		sf_complexread(ref,n1*n2,dif);
	    }
	    if (!both) {
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			mm[i2*n1+i1] = sf_cmplx(0.,0.);
			fold = 0;
			for (m=(-1*n1+1); m < n1; m++) {
			    if (m >= 0 && (i1 -m ) >= 0) {
				x1 = m;
				s1 = i2;
				x2 = i1 - m;
				s2 = m + i2;
			    } else {
				mm[i2*n1+i1] = sf_cmplx(0.,0.);
			    }
			    if (stack) {
				if (s1 >= 0 && s1 < n2 && x1 >= 0 && x1 < n1 && 
				    s2 >= 0 && s2 < n2 && x2 >= 0 && x2 < n1 ) {
#ifdef SF_HAS_COMPLEX_H
				    if(NULL != dif) {
					mm[i2*n1+i1] += dd[s1*n1+x1]*ref[s2*n1+x2];
				    } else {
					mm[i2*n1+i1] += dd[s1*n1+x1]*dd[s2*n1+x2];
				    }
#else
				    if(NULL != dif) {
					mm[i2*n1+i1] = sf_cadd(mm[i2*n1+i1],sf_cmul(dd[s1*n1+x1],ref[s2*n1+x2]));
				    } else {
					mm[i2*n1+i1] = sf_cadd(mm[i2*n1+i1],sf_cmul(dd[s1*n1+x1],dd[s2*n1+x2]));
				    }
#endif
				} else {
				    mm[i2*n1+i1] = sf_cmplx(0.,0.);
				}
				if (0.0 != cabsf(mm[i2*n1+i1])) fold++;	
			    } else {
				if (s1 >= 0 && s1 < n2 && x1 >= 0 && x1 < n1 && 
				    s2 >= 0 && s2 < n2 && x2 >= 0 && x2 < n1 ) {
#ifdef SF_HAS_COMPLEX_H
				    if(NULL != dif) {
					mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = dd[s1*n1+x1]*ref[s2*n1+x2];
				    } else {
					mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = dd[s1*n1+x1]*dd[s2*n1+x2];
				    }
#else
				    if(NULL != dif) {
					mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = sf_cmul(dd[s1*n1+x1],ref[s2*n1+x2]);
				    } else {
					mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = sf_cmul(dd[s1*n1+x1],dd[s2*n1+x2]);
				    }
#endif
				} else {
				    mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = sf_cmplx(0.,0.);
				}
			    }
			}

			if (stack) {
#ifdef SF_HAS_COMPLEX_H
			    mm[i2*n1+i1] = mm[i2*n1+i1]/(fold+SF_EPS);
#else
			    mm[i2*n1+i1] = sf_crmul(mm[i2*n1+i1],1.0/(fold+SF_EPS));
#endif
			}

		    }
		}
		
		if (!stack) {
		    tn = 0;
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    if ((i2 % jumps == 0) && (i1 % jumpo == 0)) {
				for (m=(-1*n1+1); m < n1; m++) {
				    mtemp[tn] = mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1];
				    tn ++;
				}
			    }
			}
		    }
		    if (tn!=(2*n1-1)*newn1*newn2) sf_error("jump error!");		    
		    sf_complexwrite(mtemp,tn,out);
		} else {
		    sf_complexwrite(mm,n1*n2,out);
		}
	    } else {
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			mm[i2*n1+i1] = sf_cmplx(0.,0.);
			fold = 0;
			for (m=0; m < n1; m++) {
			    x1 = m;
			    s1 = i2;
			    x2 = i1 - m + n1/2;
			    s2 = m + i2 - n1/2;
			    if (stack) {
				if (s1 >= 0 && s1 < n2 && x1 >= 0 && x1 < n1 && 
				    s2 >= 0 && s2 < n2 && x2 >= 0 && x2 < n1 ) {
#ifdef SF_HAS_COMPLEX_H
				    if(NULL != dif) {			
					mm[i2*n1+i1] += dd[s1*n1+x1]*ref[s2*n1+x2];
				    } else {
					mm[i2*n1+i1] += dd[s1*n1+x1]*dd[s2*n1+x2];
				    }
#else
				    if(NULL != dif) {	
					mm[i2*n1+i1] = sf_cadd(mm[i2*n1+i1],sf_cmul(dd[s1*n1+x1],ref[s2*n1+x2]));
				    } else {
					mm[i2*n1+i1] = sf_cadd(mm[i2*n1+i1],sf_cmul(dd[s1*n1+x1],dd[s2*n1+x2]));
				    }
#endif
				} else {
				    mm[i2*n1+i1] = sf_cmplx(0.,0.);
				}
				if (0.0 != cabsf(mm[i2*n1+i1])) fold++;				
			    } else {
				if (s1 >= 0 && s1 < n2 && x1 >= 0 && x1 < n1 && 
				    s2 >= 0 && s2 < n2 && x2 >= 0 && x2 < n1 ) {
#ifdef SF_HAS_COMPLEX_H
				    if(NULL != dif) {			
					mm[i2*n1*n1+i1*n1+m] = dd[s1*n1+x1]*ref[s2*n1+x2];
				    } else {
					mm[i2*n1*n1+i1*n1+m] = dd[s1*n1+x1]*dd[s2*n1+x2];
				    }
#else
				    if(NULL != dif) {	
					mm[i2*n1*n1+i1*n1+m] = sf_cmul(dd[s1*n1+x1],ref[s2*n1+x2]);
				    } else {
					mm[i2*n1*n1+i1*n1+m] = sf_cmul(dd[s1*n1+x1],dd[s2*n1+x2]);
				    }
#endif
				} else {
				    mm[i2*n1*n1+i1*n1+m] = sf_cmplx(0.,0.);
				}
			    }
			}
			if (stack) {
#ifdef SF_HAS_COMPLEX_H
			    mm[i2*n1+i1] = mm[i2*n1+i1]/(fold+SF_EPS);
#else
			    mm[i2*n1+i1] = sf_crmul(mm[i2*n1+i1],1.0/(fold+SF_EPS));
#endif
			}

		    }
		}
		if (!stack) {
		    tn = 0;
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    if (i2 % jumps == 0 && i1 % jumpo == 0) {
				for (m=0; m < n1; m++) {
				    mtemp[tn] = mm[i2*n1*n1+i1*n1+m];
				    tn ++;
				}
			    }
			}
		    }	
		    if (tn!=n1*newn1*newn2) sf_error("jump error!");
		    sf_complexwrite(mtemp,tn,out);
		} else {
		    sf_complexwrite(mm,n1*n2,out);
		}
	    }		
	}   
	if (verb) sf_warning(".");
    
    }

    exit(0);
}
