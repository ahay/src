/* Stretch of the time axis.

Takes: < input.rsf > output.rsf
*/

#include <math.h>
#include <string.h>
#include <float.h>

#include <rsf.h>

#include "fint1.h"

static float t0, h=0.;

static float t2(float t) { return t*t; }
static float log_frw(float t) { return logf(t/t0); }
static float log_inv(float t) { return t0*expf(t); }
static float nmo(float t) 
{ 
    t = t*t + h;
    return (t > 0.)? sqrtf(t): 0.; 
} 
static float lmo(float t) { return t + h; } 

int main(int argc, char* argv[])
{
    fint1 str;
    bool inv, half;
    int i2, n1,n2, i3, n3, n, dens, nw;
    float d1, o1, d2, o2, *trace, *stretched, h0, dh, v0;
    float (*forward)(float), (*inverse)(float);
    char *rule, *prog;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse stretching */
    if (!sf_getint("dens",&dens)) dens=1;
    /* axis stretching factor */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    prog = sf_getprog();
    if (NULL != strstr (prog, "log")) {
	rule="Log";
    } else if (NULL != strstr (prog, "t2")) {
	rule="2";
    } else if (NULL != strstr (prog, "lmo")) {
	rule="lmo";
    } else if (NULL != strstr (prog,"nmo") || 
	       NULL == (rule = sf_getstring("rule"))) {
	/* Stretch rule:
	   n - normal moveout
	   l - linear moveout
	   L - logstretch
	   2 - t^2 stretch
	*/
	rule="nmo";
    }

    if ('n'==rule[0] || 'l'==rule[0]) {
	if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input"); 
	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input"); 
	if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
	/* moveout velocity */
	sf_putfloat(out,"v0",v0);
	
	if (!sf_getbool("half",&half)) half=true;
	/* if y, the second axis is half-offset instead of full offset */
	if (half) v0 *= 0.5;

	h0 /= v0; 
	dh /= v0;
    }

    switch (rule[0]) {
	case 'n':
	    forward = inverse = nmo;
	    break;
	case 'l':
	    forward = inverse = lmo;
	    if (!sf_getfloat("delay",&t0)) sf_error("Need delay=");
	    /* time delay for rule=lmo */
	    h = -t0;

	    break;
	case 'L':
	    forward = log_frw;
	    inverse = log_inv;

	    if (o1 < FLT_EPSILON) o1=FLT_EPSILON;
	    if (!sf_getfloat ("t1", &t0) && 
		!sf_histfloat(in,"t1",&t0)) t0=o1;
	    sf_putfloat(out,"t1",t0);
	    break;
	case '2':
	    forward = t2;
	    inverse = sqrtf;
	    
	    if (o1 < FLT_EPSILON) o1=FLT_EPSILON;
	    break;
	default:
	    sf_error("rule=%s is not implemented",rule);
	    break;
    }

    if (inv) {
	if (!sf_histint(in,"nin",&n)) n=n1/dens;

	o2 = inverse(o1);
	d2 = (inverse(o1+(n1-1)*d1) - o2)/(n-1);
    } else {
	if (!sf_getint("nout",&n)) n=dens*n1;
	/* output axis length (if inv=n) */
	sf_putint(out,"nin",n1);

	o2 = forward(o1);
	d2 = o1+(n1-1)*d1;
	d2 = (forward(d2) - o2)/(n-1);
    }

    sf_putint(out,"n1",n);
    sf_putfloat(out,"o1",o2);
    sf_putfloat(out,"d1",d2);

    trace = sf_floatalloc(n1);
    stretched = sf_floatalloc(n);

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */
    str = fint1_init(nw,n1);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if ('l' == rule[0] || 'n' == rule[0]) {
		h = h0+i2*dh;
		if ('n' == rule[0]) h *= h;
		if (inv) h = -h;
	    } 

	    sf_read (trace,sizeof(float),n1,in);
	    fint1_set(str,trace);
	    
	    stretch(str,inv? forward: inverse, 
		    n1, d1, o1, n, d2, o2, stretched);

	    sf_write (stretched,sizeof(float),n,out);
	}
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mstretch.c,v 1.1 2004/04/01 02:17:14 fomels Exp $	 */

