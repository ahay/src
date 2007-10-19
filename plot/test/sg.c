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

#include <stdlib.h>
#include <stdio.h>

#include <rsfplot.h>

int main (int argc, char* argv[])
{
    char number[100];
    int	size, is, ns, ih,nh;
    float h, s, g, ds,dg,s0,g0, dgods;

    vp_init();

    s0=.5; 
    g0=.5; 

    dgods = strtod(argv[1],NULL);

    ns=13; size=10;
    nh=10;
	
    ds= (10.24-1.)/ns;
    dg = ds * dgods ;
    vp_clip( 0.,0., 10.24/.75, 10.24);

    vp_move( g0, 1.);
    vp_draw( g0, 8.5);
    vp_text( g0, 9.,  size+4, 0, "s");

    vp_move( 8., s0);
    vp_draw( 12., s0);
    vp_text( 12.5, s0, size+4, 0, "g");

    for( is=0; is<ns; is++ )  {  
	s= s0 + ds*is;
	for( ih=0; ih<nh; ih++ )  {  
	    h=      dg*ih;
	    g = s + h;
	    if( h < .1+ 9*ds) {
		sprintf(number, "%d", ih );
		vp_text( g, s, size, 0, number);
	    }
	}
    }
	
    exit(0);
}

/* 	$Id$	 */
