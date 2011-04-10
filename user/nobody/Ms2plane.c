#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int n1, n2, i1, i2;
    float a, z0, y, d2, d1, t, z, s, g, s0, w;
    float sina, cosa, sin2a, sina2, sina4, cosa2;
    float *trace;
    sf_file out;

    /* plane reflector: z = z0 + (x-x0)*tan(a) */
    /* slowness: s = sqrt(s0^2 + 2*g*z) */

    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    if (!sf_getint("n2",&n2)) sf_error("Need n2=");

    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    if (!sf_getfloat("d2",&d2)) sf_error("Need d2=");

    sf_putint(out,"n1",n1); sf_putfloat(out,"d1",d1); sf_putfloat(out,"o1",0.);
    sf_putint(out,"n2",n2); sf_putfloat(out,"d2",d2); sf_putfloat(out,"o2",0.);
    sf_setformat(out,"native_float");

    if (!sf_getfloat("s0",&s0)) sf_error("Need s0=");
    s0 *= s0;
    if (!sf_getfloat("g",&g)) g=0;

    if (!sf_getfloat("z0",&z0)) sf_error("Need z0=");
    if (!sf_getfloat("a",&a)) a=0.;
    a *= SF_PI/180.; /* degrees to radians */

    sina = sinf(a);
    cosa = cosf(a);
    sin2a = 2.*sina*cosa;
    sina2 = sina*sina;
    sina4 = sin2a*sin2a;
    cosa2 = cosa*cosa;
    
    trace = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = 0.;
	}

	if (sina==0.) {
	    z = z0;
	} else {
	    y = z0*cosa/sina + i2*d2;

	    z = s0*sin2a - 4.*g*y*sina4;

	    if (s0*s0+g*y*z < 0.) {
		sf_floatwrite(trace,n1,out);
		continue;
	    }
	    
	    z = y*sina/(1.+3.*sina2)*
		(cosa*(1.+2.*sina2) + sina*z/(s0+sqrtf(s0*s0+g*y*z)));
	}

	s = sqrtf(s0+2.*g*z);

	sf_warning("depth=%g slowness=%g",z,s);

	if (z < z0) {
	    sf_floatwrite(trace,n1,out);
	    continue;
	}

	if (s*s*cosa2 - 2.*g*z < 0.) {
	    sf_floatwrite(trace,n1,out);
	    continue;
	}

	t = sqrtf(s*s*cosa2 - 2.*g*z);

	t = 2.*z*(2.*(s*s - g*z) - (s*cosa - t)*(s*cosa - t)/3.)/(s*cosa + t);

	w = t/d1;
	i1 = w;
	w -= i1;

	if (i1>=0 && i1 < n1-1) {
	    trace[i1] = (1.-w)/t;
	    trace[i1+1] = w/t;
	}

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}

/* 	$Id$	 */
