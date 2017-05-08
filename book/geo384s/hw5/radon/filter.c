/* Mute data in the Fourier domain based on velocity and 
 * frequency. */

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nw, nx, iw, ix, i3, n3;
    float x0,dx,x, w0,dw,w, vel;
    float v1, v2, f1, f2, factor, *ctrace;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) 
      sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nw)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");

    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o1",&w0)) sf_error("No o1=");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2=");

    n3 = sf_leftsize(in,2);

    ctrace = sf_floatalloc(2*nw);

    if (!sf_getfloat("v1",&v1)) v1=0.;
    if (!sf_getfloat("v2",&v2)) v2=0.1;
    /* Velocity gate */

    if (v1>=v2) sf_error("Need v1 < v2"); 

    /* input frequencies */
    if (!sf_getfloat("f1",&f1)) f1=999;
    if (!sf_getfloat("f2",&f2)) f2=999.1;

    if (f1>=f2) sf_error("Need f1 < f2"); 


    /* Loop over shots */
    for (i3=0; i3 < n3; i3++) { 

	/* Loop over wavenumber */
	for (ix=0; ix < nx; ix++) {
	    x = fabsf(x0 + ix*dx)+dx*FLT_EPSILON;

	    sf_floatread(ctrace,2*nw,in);

	    /* Loop over frequency */
	    for (iw=0; iw < nw; iw++) {
		w = w0+iw*dw;	    
		vel = w/x; 

		if ((vel>=-v2) && (vel<=-v1)) {
		    factor=1.0f-
		      sinf(0.5*SF_PI*(vel+v2)/(v2-v1));
		} else if ((vel>=-v1) && (vel<=v1)) {
		    factor=0.0f; /* reject */
		} else if ((vel>=v1) && (vel<=v2)) {
		    factor=sinf(0.5*SF_PI*(vel-v1)/(v2-v1));
		} else {
		    factor=1.0f; /* pass */
		}

        /* low pass filter */
        if ((w >= f1) && (w <= f2)){
            factor = 1.0f - (1/(f2-f1) * (w-f1));
        } else if ((w <= -f1) && (w >= -f2)){
            factor = 1/(f2-f1) * (w-f1);
        } else if (w > f2){
            factor=0.0f;
        }
        
        
		/* real and imaginary parts */
		ctrace[2*iw]   *= factor;
		ctrace[2*iw+1] *= factor;
	    } /* iw */

	    sf_floatwrite(ctrace,2*nw,out);
	} /* ix */
    } /* i3 */

    exit(0);
}

