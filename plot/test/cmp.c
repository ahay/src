#include <math.h>

#include <rsfplot.h>

static int iseed=1996;

static float rand01(void);

int main (void)
{
    float dx, xmin=-4., xmax=7., tmax=9.0, x,t,z, x0,t0, v=1., alfa, ca, sa;
    int ix, iz, nz=100, nx=100, degrees=85;

    vp_init();

    t0 = 0;	
    x0= xmin;
    dx= (xmax-xmin)/(nx-1);
    v = 1.;
    vp_uorig( -1.+xmin, 0.);

    vp_utext(  xmin/2, tmax+.45, 18, 0, "Common Midpoint");
    vp_utext( xmax+.1, tmax-.45, 14, 0, "h");
    vp_utext( .25     ,  .1    , 14, 0, "t");

    vp_uclip (xmin, .6, xmax, tmax);
    vp_umove( xmin, tmax); 	
    vp_udraw( xmax, tmax);
    vp_umove(   0., tmax); 	
    vp_udraw( 0., tmax-tmax);

    vp_color(6);
    for (iz=0; iz < nz; iz++) {
	alfa = degrees *  2 * 3.14 / 360 * rand01(); 
	ca = cosf( alfa);
	sa = sinf( alfa);
	z =  tmax * (.1 + .9*rand01());
	z =  tmax * rand01();
	vp_penup();
	for (ix=0; ix < nx; ix++) {
	    x = x0 + ix*dx;
	    t = hypotf( z,ca*x/v);
	    if( x > -v*t &&
		x <  v*t)
		vp_upendn(x, tmax-t);
	}
    }

    return 0;
}

float rand01(void)
{
    const int ia = 727, im = 524287;
    iseed = (iseed*ia)%im;
    return ((float) iseed - 0.5)/((float) im - 1.);
}
