#include <stdlib.h>
#include <stdio.h>

#include <rsfplot.h>

int main (int argc, char* argv[])
{
    char number[100];
    int	size, is, ns,ng, ih,nh;
    float h, s, g, ds,dg,s0,g0, dgods;

    vp_init();

    s0=.5; 
    g0=.5; 

    dgods = strtod(argv[1],NULL);

    ns=20; ng=20; size=5;
    ns=10; ng=10; size=10;
    ns=13; ng=18; size=10;
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
