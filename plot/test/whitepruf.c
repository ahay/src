/* */
/*
Copyright (C) 2004 University of Texas at Austin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include <rsfplot.h>

static void filter(float x0, float y0, char *letter);
static void arrow(float x1, float y1, float x2, float y2);

static const float dx=.5, dy=.5;

int main (void)
{
    vp_init();

    vp_fat(3);
    
    vp_color(7);    
    filter(8., 4., "A");

    vp_color(6);
    filter(6., 6., "B");

    return 0;
}

static void filter(float x0, float y0, char *letter)
{
    float theta, rad=2.0, x,y,xp,yp, e=0.4;

    /* vertical half line */
    vp_move( x0+.5*dx,  y0-.5*dy);
    vp_draw( x0+.5*dx,  y0+(2.5+e)*dy);

    /* vertical lines */
    vp_move( x0- .5*dx, y0-(2.5+e)*dy);
    vp_draw( x0- .5*dx, y0+(2.5+e)*dy);

    vp_move( x0-1.5*dx, y0-(2.5+e)*dy);
    vp_draw( x0-1.5*dx, y0+(2.5+e)*dy);

    vp_move( x0-2.5*dx, y0-(2.5+e)*dy);
    vp_draw( x0-2.5*dx, y0+(2.5+e)*dy);

    /* short horiz lines */
    vp_move( x0-(2.5+e)*dx, y0-2.5*dy);
    vp_draw( x0-     .5*dx, y0-2.5*dy);

    vp_move( x0-(2.5+e)*dx, y0-1.5*dy);
    vp_draw( x0-     .5*dx, y0-1.5*dy);

    /* longer horiz lines */
    vp_move( x0-(2.5+e)*dx, y0- .5*dy);
    vp_draw( x0+     .5*dx, y0- .5*dy);

    vp_move( x0-(2.5+e)*dx, y0+ .5*dy);
    vp_draw( x0+     .5*dx, y0+ .5*dy);

    vp_move( x0-(2.5+e)*dx, y0+1.5*dy);
    vp_draw( x0+     .5*dx, y0+1.5*dy);

    vp_move( x0-(2.5+e)*dx, y0+2.5*dy);
    vp_draw( x0+     .5*dx, y0+2.5*dy);

    vp_penup();

    /* circles */
    for( theta = 90.; theta <=270; theta += 10) {
        x = x0 + rad * cosf( theta * SF_PI/180.);
	y = y0 + rad * sinf( theta * SF_PI/180.);
        vp_pendn( x, y);
    }

    /* arrows */
    for( theta = 90.+180/10. ; theta <=270; theta += 180./5.) {
        x  = x0 + rad      * cosf( theta * SF_PI/180.);
        y  = y0 + rad      * sinf( theta * SF_PI/180.);
        xp = x0 + (rad+.5) * cosf( theta * SF_PI/180.);
        yp = y0 + (rad+.5) * sinf( theta * SF_PI/180.);
        arrow( x, y, xp, yp );
    }

    vp_text(x0-.3*dx, y0-.3*dy, 10, 0, letter);
}

static void arrow(float x1, float y1, float x2, float y2)
{
    float dx, dy, r, backx, backy, perpx, perpy, tipx, tipy;
    dx = x2 - x1;
    dy = y2 - y1;
    r = hypotf(dx,dy);
    if ( r < .5)  r = .5;

    backx = -.15 * dx / r;
    backy = -.15 * dy / r;
    perpx =  .05 * dy / r;
    perpy = -.05 * dx / r;
    vp_umove( x1, y1 );         
    vp_udraw( x2, y2 );

    tipx = x2 + backx + perpx;
    tipy = y2 + backy + perpy;
    vp_umove( x2, y2 );         
    vp_udraw( tipx, tipy);

    tipx = x2 + backx - perpx;
    tipy = y2 + backy - perpy;
    vp_umove( x2, y2 );         
    vp_udraw( tipx, tipy);
}

/* 	$Id: whitepruf.c 7107 2011-04-10 02:04:14Z ivlad $	 */
