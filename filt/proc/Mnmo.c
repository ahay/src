#include <math.h>

#include <rsf.h>

#include "stretch.h"

int main (int argc, char* argv[])
{
    map nmo;
    bool inv;
    int it, ix, nt,nx;
    float dt, t0, eps, x, h, h0, f;
    float *str, *out, *inp, *vel, *off;
    sf_file cmp, nmod, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    offset = sf_input("offset");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    nx = sf_leftsize(cmp,1);

    if (!sf_getbool("inv",&inv)) inv=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getfloat ("h0",&h0)) h0=0.;

    inp = sf_floatalloc(nt);
    out = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);
    off = sf_floatalloc(nx);

    nmo = stretch_init (nt, t0, dt, nt, eps);

    sf_read (off,sizeof(float),nx,offset);
    for (ix = 0; ix < nx; ix++) {
	sf_read (vel,sizeof(float),nt,velocity);	
	x = off[ix]; 
	h = x*x - h0*h0;

	for (it=0; it < nt; it++) {
	    f = t0 + it*dt;
	    f = f*f - h*vel[it]*vel[it];
	    str[it] = (f > 0.)? sqrtf (f): t0-2.*dt;
	}

	stretch_define (nmo, str);

	sf_read (inp,sizeof(float),nt,cmp);

	if (inv) {
	    stretch_invert (nmo, out, inp);
	} else {
	    stretch_apply (nmo, inp, out);
	}

	sf_write (out,sizeof(float),nt,nmod);
    }

    exit (0);
}

