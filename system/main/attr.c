/* Display dataset attributes.

Sample output from "sfspike n1=100 | sfbandpass fhi=60 | sfattr"
*******************************************
     rms =      0.992354
    mean =      0.987576
  2-norm =       9.92354
variance =    0.00955481
 std dev =     0.0977487
     max =       1.12735 at 97
     min =      0.151392 at 100
nonzero samples = 100
  total samples = 100
*******************************************

rms                = sqrt[ sum(data^2) / n ]
mean               = sum(data) / n
norm               = sum(abs(data)^lval)^(1/lval)
variance           = [ sum(data^2) - n*mean^2 ] / [ n-1 ]
standard deviation = sqrt [ variance ]
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2009 Colorado School of Mines

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

static void location(size_t loc, size_t dim, const off_t *n);

int main(int argc, char* argv[])
{
    sf_file in=NULL;
    char *want=NULL, *buf;
    off_t n[SF_MAX_DIM], nsiz, nzero;
    int lval;
    size_t i, nbuf, nleft, dim, minloc=0, maxloc=0;
    size_t bufsiz, minloc1=0, minloc2=0, maxloc1=0, maxloc2=0;
    float f, fmin, fmax;
    double fsum, fsqr, flval, frms, fmean, fnorm, fvar, fstd;
    sf_complex c, cmin1, cmin2, cmax1, cmax2;
    sf_datatype type;

    sf_init (argc,argv);
    want = sf_getstring("want");
    /* 'all'(default), 'rms', 'mean', 'norm', 'var', 
       'std', 'max', 'min', 'nonzero', 'samples', 'short' 
        want=   'rms' displays the root mean square
        want=   'norm' displays the square norm, otherwise specified by lval.
        want=   'var' displays the variance
        want=   'std' displays the standard deviation
        want=   'nonzero' displays number of nonzero samples
        want=   'samples' displays total number of samples
        want=   'short' displays a short one-line version
     */
    if (NULL != want && 0==strcmp(want,"all")) want=NULL;

    if (!sf_getint("lval",&lval)) lval=2;
    /* norm option, lval is a non-negative integer, computes the vector lval-norm */

    in = sf_input("in");

    dim = (size_t) sf_largefiledims (in,n);
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }
    bufsiz = sf_bufsiz(in);
    buf = sf_charalloc(bufsiz);
    bufsiz /= sf_esize(in);

    type = sf_gettype (in);

    cmin1 = sf_cmplx(+FLT_MAX,0.); cmin2 = sf_cmplx(0.,+FLT_MAX);
    cmax1 = sf_cmplx(-FLT_MAX,0.); cmax2 = sf_cmplx(0.,-FLT_MAX);

    fsum = fsqr = flval = 0.0;
    fmin = +FLT_MAX;
    fmax = -FLT_MAX;
    nzero = 0;
    c = sf_cmplx(0.,0.);

    for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
	nbuf = (bufsiz < nleft)? bufsiz: nleft;
	switch (type) {
	    case SF_FLOAT: 
		sf_floatread((float*) buf,nbuf,in);
		break;
	    case SF_INT:
		sf_intread((int*) buf,nbuf,in);
		break;
            case SF_SHORT:
                sf_shortread((short*) buf,nbuf,in);
		break;
	    case SF_COMPLEX:
		sf_complexread((sf_complex*) buf,nbuf,in);
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
                case SF_SHORT:
                    f=(float) ((short*)buf)[i];
                    break;
		case SF_COMPLEX:
		    c=((sf_complex*)buf)[i];
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
	    fsqr += (double) f*f;
	    if (lval != 2 || lval != 0) flval += pow(fabs(f),lval);
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
    if (lval==2)      fnorm = sqrt(fsqr);
    else if (lval==0) fnorm = nsiz-nzero;
    else              fnorm = pow(flval,1./lval);
    frms = sqrt(fsqr/nsiz);
    if (nsiz > 1) fvar = fabs(fsqr-nsiz*fmean*fmean)/(nsiz-1);
    else          fvar = 0.0;
    fstd = sqrt(fvar);

    if(NULL==want){
	printf("******************************************* \n");
	if(SF_COMPLEX==type) 
	    printf("   rms, mean, norm, var, and std refer to amplitude\n\n");
    }
    if(NULL==want || 0==strcmp(want,"rms"))
	printf("     rms = %13.6g \n",(float) frms);
    if(NULL==want || 0==strcmp(want,"mean"))
	printf("    mean = %13.6g \n",(float) fmean);
    if(NULL==want || 0==strcmp(want,"norm"))
	printf("  %d-norm = %13.6g \n",lval,(float) fnorm);
    if(NULL==want || 0==strcmp(want,"var"))
	printf("variance = %13.6g \n",(float) fvar);
    if(NULL==want || 0==strcmp(want,"std"))
	printf(" std dev = %13.6g \n",(float) fstd);
    if(SF_COMPLEX==type) {
	if (NULL==want || 0==strcmp(want,"max")) {
	    printf("max real = (%g,%g) at ",
		   crealf(cmax1),cimagf(cmax1));
	    location(maxloc1,dim,n);
	    printf("max imag = (%g,%g) at ",
		   crealf(cmax2),cimagf(cmax2));
	    location(maxloc2,dim,n);
	}
	if (NULL==want || 0==strcmp(want,"min")) {
	    printf("min real = (%g,%g) at ",
		   crealf(cmin1),cimagf(cmin1));
	    location(minloc1,dim,n);
	    printf("min imag = (%g,%g) at ",
		   crealf(cmin2),cimagf(cmin2));
	    location(minloc2,dim,n);
	}
    } else {
	if(NULL==want || 0==strcmp(want,"max")) {
	    printf("     max = %13.6g at ",fmax);
	    location(maxloc,dim,n);
	}
	if(NULL==want || 0==strcmp(want,"min")) {
	    printf("     min = %13.6g at ",fmin);
	    location(minloc,dim,n);
	}
    }
#if defined(__cplusplus) || defined(c_plusplus)
    if(NULL==want || 0==strcmp(want,"nonzero"))
	printf("nonzero samples = %ld \n",(long int) (nsiz-nzero));
    if(NULL==want || 0==strcmp(want,"samples"))
	printf("  total samples = %ld \n",(long int) nsiz);
#else
    if(NULL==want || 0==strcmp(want,"nonzero"))
	printf("nonzero samples = %lld \n",(long long int) (nsiz-nzero));
    if(NULL==want || 0==strcmp(want,"samples"))
	printf("  total samples = %lld \n",(long long int) nsiz);
#endif
    if(NULL==want) {
	printf("******************************************* \n");
    }
    if(NULL != want && 0==strcmp(want,"short")) {
	printf("%6.2f%% zeros; min: %g; max: %g \n",
	       100.*((double)nzero)/((double)nsiz),fmin,fmax);  
    }

    exit (0);
}

static void location(size_t loc     /* liner location index */, 
		     size_t dim     /* number of dimensions */, 
		     const off_t *n /* hypercube dimensions [dim] */)
/* print out an event location */
{
    size_t i;
    off_t ni;

    for (ni=1, i=0; i < dim; ni *= n[i], i++) {
#if defined(__cplusplus) || defined(c_plusplus)
	printf("%ld ",(long int) (1+(loc/ni)%n[i]));
#else
	printf("%lld ",(long long int) (1+(loc/ni)%n[i]));
#endif
    }
    printf("\n");
}
