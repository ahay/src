#include <math.h>

#include <rsfplot.h>

int main(void)
{
    float v=1.,t,z,x,x0, theta;
    float top=3.2, c1=.9, c2=6.8;

    vp_init();

    vp_uorig (-c1,-.5);
    vp_uclip (0.,top-3.,4.,top);
    vp_umove (0.,top-0.);  
    vp_udraw (0.,top-4.);
    vp_umove (0.,top-0.);  
    vp_udraw (4.,top-0.);

    for (z=.4; z < 4.; z += .4) {
	vp_penup ();
	x0 = z * tanf( 6.283*45./360.);
	for (x=0.; x < 4.; x += .01) {
	    t = hypotf(z,x-x0)/v;
	    vp_upendn (x,top-t);
	}
    }

    vp_uorig (-c2,-.5);
    vp_uclip (0.,top-3.,4.,top);
    vp_umove (0.,top-0.);  
    vp_udraw (0.,top-4.);
    vp_umove (0.,top-0.);   
    vp_udraw (4.,top-0.);

    for(t=.4; t<6.; t += .4) {
	vp_penup ();
	x0 = t / sinf( 6.283*45./360.);
	for (theta=-89.5; theta<89.5; theta += 1.) {
	    z =      t * cosf (6.28*theta/360.);
	    x = x0 + t * sinf (6.28*theta/360.);
	    vp_upendn (x,top-z);
	}
    }

    return 0;
}
