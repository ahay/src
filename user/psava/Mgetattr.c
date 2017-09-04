/* Output dataset attributes. */
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_file Fi, Fo;
    char *want=NULL;
    
    char *buf;
    off_t n[SF_MAX_DIM], nsiz;
    size_t i, nbuf, nleft, dim;
    size_t bufsiz;

    float f, fmin, fmax, amax;
    double fsum, fsqr, frms, fmean, fvar, fstd;

    float dou[1];
    
    /*------------------------------------------------------------*/
    sf_init (argc,argv);

    want = sf_getstring("want");
    if( want==NULL ) want="amax";
    
    Fi = sf_input ("in");
    
    Fo = sf_output("out");
    sf_putint  (Fo,"n1",1);
    sf_putfloat(Fo,"d1",1.0);
    sf_putfloat(Fo,"o1",0.0);
    
    sf_putint  (Fo,"n2",1);
    sf_putfloat(Fo,"d2",1.0);
    sf_putfloat(Fo,"o2",0.0);

    sf_putint  (Fo,"n3",1);
    sf_putfloat(Fo,"d3",1.0);
    sf_putfloat(Fo,"o3",0.0);

    sf_putint  (Fo,"n4",1);
    sf_putfloat(Fo,"d4",1.0);
    sf_putfloat(Fo,"o4",0.0);
    
    /*------------------------------------------------------------*/
    
    /* get dimensions */
    dim = (size_t) sf_largefiledims (Fi,n);
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }
    
    bufsiz = sf_bufsiz(Fi);
    buf = sf_charalloc(bufsiz);
    bufsiz /= sf_esize(Fi);

    /*------------------------------------------------------------*/
    
    fsum = fsqr = 0.0;
    fmin = +FLT_MAX;
    fmax = -FLT_MAX;
    amax = 0.;

    for (nleft = nsiz; nleft > 0; nleft -= nbuf) {
	nbuf = (bufsiz < nleft)? bufsiz: nleft;

	sf_floatread((float*) buf,nbuf,Fi);
		
	for (i=0; i < nbuf; i++) {
	    f = ((float*)buf)[i];

	    if( 0==strcmp(want, "amax") )
		if (fabs(f) > fabs(amax)) amax = f;
		    
	    if( 0==strcmp(want, "max") )
		if (f > fmax) fmax = f;

	    if( 0==strcmp(want, "min") )
		if (f < fmin) fmin = f;

	    if( 0==strcmp(want,"mean") || 0==strcmp(want, "var") || 0==strcmp(want, "std") )
		fsum += f;

	    if( 0==strcmp(want,"rms") || 0==strcmp(want, "var") || 0==strcmp(want, "std") )
		fsqr += (double) f*f;

	}
    }

    /*------------------------------------------------------------*/
    if(0==strcmp(want,"amax")) dou[0] = amax;
    if(0==strcmp(want, "max")) dou[0] = fmax;
    if(0==strcmp(want, "min")) dou[0] = fmin;

    if( 0==strcmp(want,"mean") || 0==strcmp(want, "var") || 0==strcmp(want, "std") )
	fmean = fsum/nsiz;
    if(0==strcmp(want,"mean")) dou[0] = fmean;
    
    if( 0==strcmp(want, "rms") || 0==strcmp(want, "var") || 0==strcmp(want, "std") )
	frms = sqrt(fsqr/nsiz);
    if(0==strcmp(want, "rms")) dou[0] = frms;
    
    if( 0==strcmp(want, "var") || 0==strcmp(want, "std")) {
	if (nsiz > 1) fvar = fabs(fsqr-nsiz*fmean*fmean)/(nsiz-1);
	else          fvar = 0.0;
    }
    if(0==strcmp(want, "var")) dou[0] = fvar;

    if( 0==strcmp(want, "std") )
	fstd = sqrt(fvar);
    if(0==strcmp(want, "std")) dou[0] = fstd;

    /*------------------------------------------------------------*/
    
    sf_floatwrite(dou,1,Fo);
	
    exit (0);
}

