/* Half thresholding using exact value.*/

/*
  Copyright (C) 2013 University of Texas at Austin

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
  
  Reference:
  Pengliang Yang et al, 2013, An iterative half thresholding method for seismic data interpolation,
  83rd Annual International Meeting, SEG, Expanded Abstracts, 3579-3582.
  
*/
#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, n, n1;		/*index and length of the total array and temporary integer*/ 
    int ifverb;			/* if print the threshold value */
    int ifperc;			/* if use percentile thresholding */
    float *dat=NULL;		/* float array*/
    float *adat=NULL; 		/* float array (for absolute value )*/
    float *diff;		/* float array (for difference )*/
    float t;			/* threshold value */
    float tau;			/* regularization term : t= 3/2tau^(2/3) */ 
    float thr;			/* input threshold value (redundent) */
    float d;			/* absolute value (temporary) */
    float multiplier=0;    	/* thresholding multiplier */	
    char *thre;			/*name of the difference section file*/
    sf_complex *cdat=NULL;	/*complex array*/
    sf_complex *cdiff=NULL;	/*complex array (for difference) */
    sf_file in=NULL;		/*input file*/
    sf_file out=NULL; 		/*output thresholded file*/
    sf_file other=NULL;		/*output difference section*/


/***************************************************/
/*	Initialization				   */
/***************************************************/
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

/***************************************************/
/*		Allocate memory 		   */
/***************************************************/

    n = sf_filesize(in);
    adat = sf_floatalloc(n);
    
    
/***************************************************/
/*	Getting parameters 		   	*/
/***************************************************/
    if(!sf_getint("ifverb",&ifverb)) ifverb=0;
    /* 0, not print threshold value; 1, print threshold value. */

    if(!sf_getint("ifperc",&ifperc)) ifperc=1;
    /* 0, exact-value thresholding; 1, percentile thresholding. */

    if (!sf_getfloat("thr",&thr)) sf_error("Need thr=");
    /* thresholding level */ 
    
    if (NULL != (thre=sf_getstring("other"))){other = sf_output("other");}
    /* If output the difference between the thresholded part and the original one */
    
/***************************************************/
/*Deciding thresholding level, data type and reading data.*/
/***************************************************/
    if(ifperc==1)
    {
    	n1 = 0.5+n*(1.-0.01*thr);
    	if (n1 < 0) n1=0;
    	if (n1 >= n) n1=n-1;

    	if (SF_FLOAT == sf_gettype(in)) {
		dat = sf_floatalloc(n);
		sf_floatread(dat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = fabsf(dat[i]);
		}
    	} else if (SF_COMPLEX == sf_gettype(in)) {
		cdat = sf_complexalloc(n);
		sf_complexread(cdat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = cabsf(cdat[i]);
		}
    	} else {
		sf_error("Need float or complex input");
    	}
    	t = sf_quantile(n1,n,adat); /* lambda = 3/2 * tau ^ 2/3 */
    }
    else 
    {
    	if (SF_FLOAT == sf_gettype(in)) {
		dat = sf_floatalloc(n);
		sf_floatread(dat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = fabsf(dat[i]);
		}
    	} else if (SF_COMPLEX == sf_gettype(in)) {
		cdat = sf_complexalloc(n);
		sf_complexread(cdat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = cabsf(cdat[i]);
		}
    	} else {
		sf_error("Need float or complex input");
    	}
    	t=thr;			/* lambda = 3/2 * tau ^ 2/3 */
    }

    tau=powf(2.0/3*t,3/2); 	/* tau = (lambda/1.5) ^ 3/2 */

    if(ifverb==1) sf_warning("Threshold=%g",t);   


/***************************************************/
/* Implementing thresholding and writing output	   */
/***************************************************/
    if (NULL != dat) {
	if(thre !=NULL ) {diff=sf_floatalloc(n);}
	for (i=0; i < n; i++) {
	    d = dat[i];
	    if (d >= -t && d<=t ) {
		dat[i] = 0.;
		if(thre!=NULL) diff[i]=d;
	    }
	    else {
	    multiplier=2.0/3*(1+cosf(2.0*SF_PI/3-2.0/3*acosf(tau/8*powf(adat[i]/3,-1.5))));
	    dat[i] *= multiplier;
	    }   
	}
	if(thre != NULL) 
		sf_floatwrite(diff,n,other); /* write the difference */		
	sf_floatwrite(dat,n,out);
    } else {
	if(thre !=NULL ){ cdiff=sf_complexalloc(n);}
	for (i=0; i < n; i++) {
	    d = cabsf(cdat[i]);	
	    if (d >= -t && d<=t) {
		if(thre!=NULL)   cdiff[i] = cdat[i];
		cdat[i] = sf_cmplx(0.,0.);
	    } else {
#ifdef SF_HAS_COMPLEX_H
		multiplier=(2.0/3.0)*(1.0+cosf(2.0*SF_PI/3-2.0/3*acosf(tau/8.0*powf(fabsf(d)/3,-1.5))));
		cdat[i] *= multiplier;
#else
		multiplier=(2.0/3.0)*(1.0+cosf(2.0*SF_PI/3-2.0/3*acosf(tau/8.0*powf(fabsf(d)/3,-1.5))));
		cdat[i] = sf_crmul(cdat[i],multiplier);
#endif
	    }

	}
	
	if(thre!=NULL) sf_complexwrite(cdiff,n,other);
	sf_complexwrite(cdat,n,out);
    }

    exit(0);
}
