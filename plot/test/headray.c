#include <stdlib.h>
#include <math.h>

#include <rsfplot.h>

static void arrow( float x1, float y1, float x2, float y2);

int main (void)
{
    float top=10.;
    float xmin=0.;
    float ymin=0.;
    float xmax = 13.;
    float ymax = 3.;
    float d = 2.;

    vp_uorig (-.1,-.1);
    vp_color (4);

/* axes */
    vp_fat(7);
    vp_umove( xmin, top-  0.);  
    vp_udraw( xmax, top-  0.);
    vp_umove( 0.,   top-ymin);  
    vp_udraw( 0.,   top-ymax);
    vp_fat(3);

/* rays */
    vp_color(7);
    vp_umove( xmin,           top-ymin);
    vp_udraw( xmin+d,         top-( ymin+d ));	/* diagonal down */
    vp_udraw( xmin+d+d+d+d+d, top-( ymin+d ));	/* along bed */

    arrow(xmin+3*d, top-(ymin+d),
	  xmin+4*d, top-(ymin  ));
    arrow(xmin+4*d, top-(ymin+d),
	  xmin+5*d, top-(ymin  ));
    arrow(xmin+5*d, top-(ymin+d),
	  xmin+6*d, top-(ymin  ));

    exit(0);
}

static void arrow( float x1, float y1, float x2, float y2)
{
    float dx, dy, r, backx, backy, perpx, perpy, tipx, tipy, noise;
    int i, n=200;

    dx = x2 - x1;
    dy = y2 - y1;
    r = hypotf(dx,dy);
    if( r < .5)  r = .5;

    backx = -.40 * dx / r;
    backy = -.40 * dy / r;
    perpx =  .14 * dy / r;
    perpy = -.14 * dx / r;

    vp_umove(x1,y1);
    noise = 1/600.;
    for (i=0; i < n; i++) {
	vp_udraw(x1 + ((i+1)*(x2-x1))/n + noise, 
		 y1 + ((i+1)*(y2-y1))/n + noise);
	noise = -noise;
    }

    tipx = x2 + backx + perpx;
    tipy = y2 + backy + perpy;
    vp_umove(x2, y2);		
    vp_udraw(tipx,tipy);

    tipx = x2 + backx - perpx;
    tipy = y2 + backy - perpy;
    vp_umove(x2, y2);		
    vp_udraw(tipx, tipy);
}
