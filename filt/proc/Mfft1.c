/* Fast Fourier Transform along the first axis.

Takes: < input.rsf > output.rsf
*/

#include <rsf.h>

int main (int argc, char *argv[])
{
    bool cos, inv;
    int n1, nt, nw, i1, i2, n2;
    float dw, *p, *cc=NULL, d1, o1;
    float complex *pp;
    sf_file in, out;

    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("cos",&cos)) cos=false;
    /* if y, perform cosine transform */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse transform */
    
    if (cos) {  
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    } else if (inv) {
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype (out,SF_FLOAT);
    } else {
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    	sf_settype (out,SF_COMPLEX);
    }

    n2 = sf_leftsize(in,1);

    if (!inv) {
	if (!sf_histint(in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	/* determine wavenumber sampling (for real to complex FFT) */
	nt = sf_npfar(cos? 2*(n1-1): n1);
	nw = nt/2+1;
	dw = 1./(nt*d1);

	sf_putint(out,"n1",nw);
	sf_putfloat(out,"o1",0.);
	sf_putfloat(out,"d1",dw);

	sf_putfloat(out,"t0",o1);
    } else {
	if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"t0",&o1)) o1=0.; 

	nt = 2*(nw-1);
	d1 = 1./(nt*dw);
	n1 = cos? 1+nt/2:nt;

	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
    }	
    
    p = sf_floatalloc(nt);
    pp = sf_complexalloc(nw);
    if (cos) cc = sf_floatalloc(nw);

    for (i2=0; i2 < n2; i2++) {
	if (!inv) {
	    sf_read (p,sizeof(float),n1,in);
	    if (!cos) {
		for (i1=n1; i1 < nt; i1++) {
		    p[i1]=0.0;
		}
	    } else {
		for (i1=n1; i1 <= nt/2; i1++) {
		    p[i1]=0.0;
		}
		for (i1=nt/2+1; i1 < nt; i1++) {
		    p[i1] = p[nt-i1];
		}
	    }
	    
	    sf_pfarc(1,nt,p,pp);
	    
	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    pp[i1] *= cexpf(I*2.0*SF_PI*i1*dw*o1);
		}
	    }
	    
	    if (cos) {
		for (i1=0; i1 < nw; i1++) {
		    cc[i1] = crealf(pp[i1]);
		}
		sf_write(cc,sizeof(float),nw,out);
	    } else {
		sf_write(pp,sizeof(float complex),nw,out);
	    }
	} else {
	    if (cos) {
		sf_read(cc,sizeof(float),nw,in);
		for (i1=0; i1 < nw; i1++) {
		    pp[i1] = cc[i1];
		}
	    } else {	    
		sf_read(pp,sizeof(float complex),nw,in);
	    }

	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    pp[i1] *= cexpf(-I*2.0*SF_PI*i1*dw*o1);
		}
	    }

	    sf_pfacr(-1,nt,pp,p);

	    for (i1=0; i1 < n1; i1++) {
		p[i1] /= nt;
	    }

	    sf_write (p,sizeof(float),n1,out);
	}
    }
    
    exit (0);
}

/* 	$Id: Mfft1.c,v 1.8 2004/03/13 06:00:33 fomels Exp $	 */
