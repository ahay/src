#include <math.h>

#include <rsfplot.h>

static int iseed=1996;

static float rand01(void);

int main(void)
{
    int ix,nx=100, iz, nz=100, degrees=85;
    float dx, xmin=-4., xmax=7., tmax=9., x,t, x0,t0=0.;
    float alfa, c2a, orig, arg;

    t0 = 0;	
    x0= xmin;
    dx= (xmax-xmin)/(nx-1);
    vp_uorig( -1.+xmin, 0.);

    vp_utext(  xmin/2, tmax+.45, 18, 0, "Common Shot");
    vp_utext( xmax+.1, tmax-.45, 14, 0, "g");
    vp_utext( .25     ,  .1    , 14, 0, "t");

    vp_uclip (xmin, .6, xmax, tmax);
    vp_umove( xmin, tmax); 	
    vp_udraw( xmax, tmax);
    vp_umove(   0., tmax); 	
    vp_udraw( 0., tmax-tmax);
    vp_color(6);
    for (iz=0; iz < nz; iz++) {
	orig = 3 * ((xmax-xmin) * rand01() +xmin);
	alfa = degrees *  2 * 3.14 / 360 * rand01(); 
	c2a = cos( 2*alfa);
	vp_penup();
	for (ix=0; ix < nx; ix++) {			/* x=g; s=0 */
	    x = x0 + ix*dx;
	    arg = orig*orig +(x-orig)*(x-orig) + 2*orig*(x-orig) * c2a;
	    t = sqrtf( arg); 
	    if( t < fabsf(x)) t = fabsf(x);
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
