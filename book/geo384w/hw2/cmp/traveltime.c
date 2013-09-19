/* Compute traveltime in a V(z) model. */
#include <rsf.h>

int main(int argc, char* argv[])
{
    char *type;
    int ih, nh, it, nt, ir, nr, *r, iter, niter;
    float h, dh, h0, dt, t0, t2, h2, v2, s, p, hp, tp;
    float *v, *t;
    sf_file vel, tim;

    /* initialize */
    sf_init(argc,argv);

    /* input and output */
    vel = sf_input("in");
    tim = sf_output("out");
    
    /* time axis from input */
    if (!sf_histint(vel,"n1",&nt)) sf_error("No n1=");
    if (!sf_histfloat(vel,"d1",&dt)) sf_error("No d1=");

    /* offset axis from command line */
    if (!sf_getint("nh",&nh)) nh=1;
    /* number of offsets */
    if (!sf_getfloat("dh",&dh)) dh=0.01;
    /* offset sampling */
    if (!sf_getfloat("h0",&h0)) h0=0.0;
    /* first offset */

    /* get reflectors */
    if (!sf_getint("nr",&nr)) nr=1;
    /* number of reflectors */
    r = sf_intalloc(nr);
    if (!sf_getints("r",r,nr)) sf_error("Need r=");

    if (NULL == (type = sf_getstring("type")))
	type = "hyperbolic";
    /* traveltime computation type */

    if (!sf_getint("niter",&niter)) niter=10;
    /* maximum number of shooting iterations */

    /* put dimensions in output */
    sf_putint(tim,"n1",nh);
    sf_putfloat(tim,"d1",dh);
    sf_putfloat(tim,"o1",h0);
    sf_putint(tim,"n2",nr);

    /* read velocity */
    v = sf_floatalloc(nt);
    sf_floatread(v,nt,vel);
    /* convert to velocity squared */
    for (it=0; it < nt; it++) {
	v[it] *= v[it];
    }

    t = sf_floatalloc(nh);

    for (ir=0; ir<nr; ir++) {
	nt = r[ir];
	t0 = nt*dt; /* zero-offset time */
	t2 = t0*t0;

	p = 0.0;

	for (ih=0; ih<nh; ih++) {
	    h = h0+ih*dh; /* offset */
	    h2 = h*h; 

	    switch(type[0]) {
		case 'h': /* hyperbolic approximation */
		    v2 = 0.0;
		    for (it=0; it < nt; it++) {
			v2 += v[it];
		    }
		    v2 /= nt;

		    t[ih] = sqrtf(t2+h2/v2);
		    break;
		case 's': /* shifted hyperbola */

		    /* !!! MODIFY BELOW !!! */

		    s = 0.0;

		    v2 = 0.0;
		    for (it=0; it < nt; it++) {
			v2 += v[it];
		    }
		    v2 /= nt;

		    t[ih] = sqrtf(t2+h2/v2);
		    break;
		case 'e': /* exact */
		    
		    /* !!! MODIFY BELOW !!! */

		    for (iter=0; iter < niter; iter++) {
			hp = 0.0;
			for (it=0; it < nt; it++) {
			    v2 = v[it];
			    hp += v2/sqrtf(1.0-p*p*v2);
			}
			hp *= p*dt;

			/* !!! SOLVE h(p)=h !!! */
		    }

		    tp = 0.0;
		    for (it=0; it < nt; it++) {
			v2 = v[it];
			tp += dt/sqrtf(1.0-p*p*v2);
		    }
		    
		    t[ih] = tp;
		    break;
		default:
		    sf_error("Unknown type");
		    break;
	    }
	}

	sf_floatwrite(t,nh,tim);
    }

    exit(0);
}
