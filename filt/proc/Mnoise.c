/* Add random noise to the data.

Takes: < data.rsf > noisy.rsf
*/

#include <math.h>
#include <time.h>

#include <rsf.h>

#include "randn.h"

int main (int argc, char* argv[])
{
    float mean, var, range, a, b, *dat;
    int nbuf, nsiz, seed, i;
    bool normal, rep;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    srand(seed);

    if (!sf_getbool("type",&normal)) normal=true;
    /* noise distribution, y: normal, n: uniform */
    if (!sf_getfloat("var",&var)) {
	/* noise variance */
	if (!sf_getfloat("range",&range)) {
	    /* noise range (default=1) */
	    a = 1.;
	} else {
	    a = normal? 2.*range/9. : 2.*range;
	}
    } else {
	a = normal? sqrt(var): sqrtf(12*var);
    }

    if (!sf_getfloat("mean",&mean)) mean=0;
    /* noise mean */
    b = normal? mean: mean - 0.5*a;

    if (!sf_getbool("rep",&rep)) rep=false;
    /* if y, replace data with noise */

    nbuf = BUFSIZ/sizeof(float);
    dat = sf_floatalloc (nbuf);

    for (nsiz = sf_filesize(in); nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;

	if (rep) {
	    if (normal) {
		for (i=0; i < nbuf; i++) {
		    dat[i] = a*randn_one() + b;
		}
	    } else {
		for (i=0; i < nbuf; i++) {
		    dat[i] = a*random_one() + b;
		}
	    }
	} else {
	    sf_read(dat,sizeof(float),nbuf,in);
	    
	    if (normal) {
		for (i=0; i < nbuf; i++) {
		    dat[i] += a*randn_one() + b;
		}
	    } else {
		for (i=0; i < nbuf; i++) {
		    dat[i] += a*random_one() + b;
		}
	    }
	}

	sf_write(dat,sizeof(float),nbuf,out);  
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mnoise.c,v 1.4 2004/03/22 05:43:25 fomels Exp $	 */
