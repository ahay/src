/* Compute traveltime in a V(z) model. */

#include <rsf.h>

int main(int argc, char* argv[])
{
    char *type;
    int ih, nh, it, nt, ir, nr, *r, iter, niter;
    float h, dh, h0, dt, t0, t2, h2, v2, va, p, hp, ghp, tp;
    float *v, *t;
    sf_file vel, tim;

    /**************************************/
    /* ADD variables for auxiliary files! */
    /**************************************/

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

    /***************************************/
    /* CHANGE to read reflector from file! */
    /***************************************/
    /* reflector axis from command line */
    if (!sf_getint("nr",&nr)) nr=1;
    /* number of reflectors */
    r = sf_intalloc(nr);
    if (!sf_getints("r",r,nr)) sf_error("Need r=");

    if (NULL == (type = sf_getstring("type")))
	type = "accelerate";
    /* traveltime computation type */

    if (!sf_getint("niter",&niter)) niter=10;
    /* maximum number of Newton iterations */

    /* put dimensions in output */
    sf_putint(tim,"n1",nh);
    sf_putfloat(tim,"o1",h0);
    sf_putfloat(tim,"d1",dh);
    sf_putstring(tim,"label1","Offset");
    sf_putstring(tim,"unit1","km");
    sf_putint(tim,"n2",nr);
    sf_putfloat(tim,"o2",1.);
    sf_putfloat(tim,"d2",1.);
    sf_putstring(tim,"label2","Reflector");
    sf_putstring(tim,"unit2","number");

    /***********************************/
    /* ADD lines for parameter output! */
    /***********************************/

    /* read velocity */
    v = sf_floatalloc(nt);
    sf_floatread(v,nt,vel);

    /* convert to velocity squared */
    for (it=0; it < nt; it++) {
	v[it] *= v[it];
    }

    /* allocate temporary memory */
    t = sf_floatalloc(nh);

    for (ir=0; ir<nr; ir++) { /* loop over reflectors */
	nt = r[ir];

	t0 = nt*dt;
	t2 = t0*t0;

	p = 0.0;

	for (ih=0; ih<nh; ih++) { /* loop over offsets */
	    h = h0+ih*dh;
	    h2 = h*h; 

	    switch(type[0]) {		
		case 'a': /* accelerated RMS approximation */
		    v2 = 0.0;
		    va = 0.0;
		    for (it=0; it < nt; it++) {
			v2 += v[it];
			va += v[it]*v[it];
		    }
		    v2 /= nt;
		    va /= nt;

		    p = -1/(4*t2)+va/(4*t2*v2*v2);

		    t[ih] = sqrtf(t2+h2/(v2+p*h2));
		    break;

		case 'e': /* exact */
		    for (iter=0; iter < niter; iter++) {
			hp = 0.0;
			ghp = 0.0;
			for (it=0; it < nt; it++) {
			    v2 = v[it];
			    hp += v2/sqrtf(1.0-p*p*v2);
			    ghp += v2/(sqrtf(1.0-p*p*v2)*(1.0-p*p*v2));
			}
			hp *= p*dt;
			ghp *= dt;

			p = p-(hp-h)/ghp;
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

	/* write output */
	sf_floatwrite(t,nh,tim);

	/********************************/
	/* ADD line to write parameter! */
	/********************************/
    }

    exit(0);
}
