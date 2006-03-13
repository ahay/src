/* Display dataset attributes.
   
Sample output from "sfspike n1=100 | sfbandpass fhi=60 | sfattr"
******************************************* 
rms = 0.992354 
mean value = 0.987576 
norm value = 9.92354 
variance = 0.00955481 
standard deviation = 0.0977487 
maximum value = 1.12735 at 97 
minimum value = 0.151392 at 100 
number of nonzero samples = 100 
total number of samples = 100 
******************************************* 

rms                = sqrt[ sum(data^2) / n ]
mean               = sum(data) / n
norm               = sqrt[ sum(data^2) ]
variance           = [ sum(data^2) - n*mean^2 ] / [ n-1 ]
standard deviation = sqrt [ variance ]
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include <rsf.h>

static void location(size_t loc, size_t dim, const int *n);

int main(int argc, char* argv[])
{
    sf_file in;
    char *want, buf[BUFSIZ];
    int n[SF_MAX_DIM], nzero;
    size_t i, nsiz, nbuf, nleft, dim, minloc=0, maxloc=0;
    size_t bufsiz=BUFSIZ, minloc1=0, minloc2=0, maxloc1=0, maxloc2=0;
    float f, fmin, fmax;
    double fsum, fsqr, frms, fmean, fnorm, fvar, fstd;
    float complex c=0., cmin1=0., cmin2=0., cmax1=0., cmax2=0.;
    sf_datatype type;
    
    sf_init (argc,argv);
    want = sf_getstring("want");
    /* 'all' (default), 'rms', 'mean', 'norm', 'var', 'std', 'max', 'min', 'short' 
        want=   'rms' displays the root mean square
        want=   'var' displays the variance
        want=   'std' displays the standard deviation
        want=   'short' displays a short one-line version
     */ 
    if (NULL != want && 0==strcmp(want,"all")) want=NULL;
    
    in = sf_input("in");

    dim = (size_t) sf_filedims (in,n);
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }
    bufsiz /= sf_esize(in);

    type = sf_gettype (in);
    if (SF_COMPLEX==type) {
	cmin1 = +FLT_MAX; cmin2 = +I*FLT_MAX;
	cmax1 = -FLT_MAX; cmax2 = -I*FLT_MAX;
    }
    
    fsum = fsqr = 0.0;
    fmin = +FLT_MAX;
    fmax = -FLT_MAX;
    nzero = 0;

    for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
	nbuf = (bufsiz < nleft)? bufsiz: nleft;
	switch (type) {
	    case SF_FLOAT: 
		sf_floatread((float*) buf,nbuf,in);
		break;
	    case SF_INT:
		sf_intread((int*) buf,nbuf,in);
		break;
	    case SF_COMPLEX:
		sf_complexread((float complex*) buf,nbuf,in);
		break;
	    case SF_UCHAR:
		sf_ucharread((unsigned char*) buf,nbuf,in);
		break;
	    case SF_CHAR:
	    default:
		sf_charread(buf,nbuf,in);
		break;
	}
	for (i=0; i < nbuf; i++) {
	    switch (type) {
		case SF_FLOAT: 
		    f=((float*)buf)[i]; 
		    break;
		case SF_INT: 
		    f=(float) ((int*)buf)[i]; 
		    break;
		case SF_COMPLEX: 
		    c=((float complex*)buf)[i];  
		    f=cabsf(c); 
		    break;
		case SF_UCHAR:
		    f=(float) ((unsigned char*) buf)[i]; 
		    break;
		case SF_CHAR:
		default: 
		    f=(float) buf[i]; 
		    break;
	    }
	    fsum += f;
	    fsqr += f*f;
	    if (0. == f) nzero++;
	    if (SF_COMPLEX==type) {
		if (crealf(c) > crealf(cmax1)) {
		    cmax1 = c;
		    maxloc1 = nsiz-nleft+i;
		}
		if (cimagf(c) > cimagf(cmax2)) {
		    cmax2 = c;
		    maxloc2 = nsiz-nleft+i;
		}
		if (crealf(c) < crealf(cmin1)) {
		    cmin1 = c;
		    minloc1 = nsiz-nleft+i;
		}
		if (cimagf(c) < cimagf(cmin2)) {
		    cmin2 = c;
		    minloc2 = nsiz-nleft+i;
		}
	    } else {
		if (f > fmax) {
		    fmax = f;
		    maxloc = nsiz-nleft+i;
		}
		if (f < fmin) {
		    fmin = f;
		    minloc = nsiz-nleft+i;
		}
	    }
	}
    }
    fmean = fsum/nsiz;
    fnorm = sqrt(fsqr);
    frms  = sqrt(fsqr/nsiz);
    if (nsiz > 1) fvar = (fsqr-nsiz*fmean*fmean)/(nsiz-1);
    else          fvar = 0.0;
    fstd = sqrt(fvar);

    if( NULL==want){
	printf("******************************************* \n");
	if(SF_COMPLEX==type) 
	    printf("   rms, mean, norm, var, and std refer to amplitude\n\n");
    }
    if(NULL==want || 0==strcmp(want,"rms"))
	printf("rms = %g \n",(float) frms);
    if(NULL==want || 0==strcmp(want,"mean"))
	printf("mean value = %g \n",(float) fmean);
    if(NULL==want || 0==strcmp(want,"norm"))
	printf("norm value = %g \n",(float) fnorm);
    if(NULL==want || 0==strcmp(want,"var"))
	printf("variance = %g \n",(float) fvar);
    if(NULL==want || 0==strcmp(want,"std"))
	printf("standard deviation = %g \n",(float) fstd);
    if(SF_COMPLEX==type) {
	if (NULL==want || 0==strcmp(want,"max")) {
	    printf("maximum real value      = (%g,%g) at ",
		   crealf(cmax1),cimagf(cmax1));
	    location(maxloc1,dim,n);
	    printf("maximum imaginary value = (%g,%g) at ",
		   crealf(cmax2),cimagf(cmax2));
	    location(maxloc2,dim,n);
	}
	if (NULL==want || 0==strcmp(want,"min")) {
	    printf("minimum real value      = (%g,%g) at ",
		   crealf(cmin1),cimagf(cmin1));
	    location(minloc1,dim,n);
	    printf("minimum imaginary value = (%g,%g) at ",
		   crealf(cmin2),cimagf(cmin2));
	    location(minloc2,dim,n);
	}
    } else {
	if(NULL==want || 0==strcmp(want,"max")) {
	    printf("maximum value = %g at ",fmax);
	    location(maxloc,dim,n);
	}
	if(NULL==want || 0==strcmp(want,"min")) {
	    printf("minimum value = %g at ",fmin);
	    location(minloc,dim,n);
	}
    }    
    if(NULL==want) {
	printf("number of nonzero samples = %d \n",(int) nsiz-nzero);
	printf("total number of samples = %d \n",(int) nsiz);
	printf("******************************************* \n");
    }
    if(NULL != want && 0==strcmp( want,"short")) {
	printf("%6.2f%% zeros min: %g max: %g \n",
	       100.*((double)nzero)/((double)nsiz),fmin,fmax);  
    }
    
    exit(0);
}

static void location(size_t loc, size_t dim, const int *n)
{
    size_t i;
    int ni;

    for (ni=1, i=0; i < dim; ni *= n[i], i++) {
      printf("%d ",(int) (1+(loc/ni)%n[i]));
    }
    printf("\n");
}

/* 	$Id$	 */

