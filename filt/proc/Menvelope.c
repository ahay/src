/* Compute data envelope.

Takes: < data.rsf > envelope.rsf
*/

#include <rsf.h>

#include "triangle.h"

int main (int argc, char* argv[])
{
    bool freq;
    int n1,n2,n3, i1,i2,i3, tc1, tc2, nw;
    float *data, **bot=NULL, **top=NULL, den;
    float complex *cdat, *ctop=NULL;
    triangle tr1=NULL, tr2=NULL;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1; 
    n3 = sf_leftsize(in,2);

    /* determine frequency sampling (for real to complex FFT) */
    nw = sf_npfa(n1);
     
    if (!sf_getbool("freq",&freq)) freq=false;
    /* if y, compute instantenous frequency */

    data = sf_floatalloc(n1);
    cdat = sf_complexalloc (nw);

    if (freq) {
	ctop = sf_complexalloc (nw);
	top = sf_floatalloc2(n1,n2);
	bot = sf_floatalloc2(n1,n2);

	if (!sf_getint("tc1",&tc1)) tc1=1;
	if (!sf_getint("tc2",&tc2)) tc2=1;
	/* smoothing triangle size for instanteneous frequency (if freq=y) */
    
	tr1 = triangle_init (tc1,n1);
	tr2 = triangle_init (tc2,n2);
    }

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_read(data,sizeof(float),n1,in);
	    for (i1=0; i1 < n1; i1++) {
		cdat[i1]=data[i1];
	    }
	    for (i1=n1; i1 < nw; i1++) {
		cdat[i1]=0.;
	    }
	    sf_pfacc (1,nw,cdat);
	    cdat[0] *= 0.5;
	    cdat[nw/2] *= 0.5;
	    for (i1=nw/2+1; i1 < nw; i1++) {
		cdat[i1]=0.;
	    }
	    if (freq) {
		for (i1=0; i1 < nw; i1++) {
		    ctop[i1] = (2.*SF_PI*I*i1/n1) * cdat[i1];
		}
		sf_pfacc (-1,nw,ctop);
		sf_pfacc (-1,nw,cdat);
		for (i1=0; i1 < n1; i1++) {
		    bot[i2][i1] = crealf(conjf(cdat[i1])*cdat[i1]);
		    top[i2][i1] = crealf(-conjf(cdat[i1])*ctop[i1]*I);
		}
		smooth(tr1,i2*n1,1,false,top[0]);
		smooth(tr1,i2*n1,1,false,bot[0]);
	    } else {
		sf_pfacc (-1,nw,cdat);
		for (i1=0; i1 < n1; i1++) {
		    data[i1] = 2.*cabsf(cdat[i1])/nw;
		}
		sf_write(data,sizeof(float),n1,out);
	    }
	}
	if (freq) {
	    smooth(tr2,0,n1,false,bot[0]);
	    smooth(tr2,0,n1,false,top[0]);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    den = bot[i2][i1];
		    if (den != 0.) top[i2][i1] /= den;
		}
	    }
	    sf_write(top[0],sizeof(float),n1*n2,out);
	}
    }

    exit(0);
}

/* 	$Id: Menvelope.c,v 1.4 2003/11/17 19:42:01 fomels Exp $	 */
