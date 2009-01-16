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

#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int job;
    float x,y, rads, theta;
    sf_complex cs, cz;
    const int nfat=1, ntheta= 271;
    const float eps=0.1, k2=0.05;

    vp_init();

    for (job=0; job < 5; job++) {
	vp_uorig ( -1.75-2.20*job, -1.5 );
	vp_uclip (-1.,-1.,1.3,1.1);
	vp_fat (nfat);

	vp_umove (0.,-.9);  
	vp_udraw (0.,4.9);
	vp_umove (-.5,0.);  
	vp_udraw (1.2,0.);

	vp_utext(  .1,  .9,  4,0, "Im");
	vp_utext( 1.05, -.2, 4,0, "Re");

	switch (job) {
	    case 0:
		vp_utext( .1, -.7 , 4,0, "\\F10 w=+p");
		vp_utext( .1,  .6 , 4,0, "\\F10 w=-p");
		vp_utext( .1,  .05, 4,0, "\\F10 w=0");
		break;
	    case 1:
		vp_utext( .4, -.65, 4,0, "\\F10 w=p/2");
		vp_utext( .4,  .55, 4,0, "\\F10 w=-p/2");
		vp_utext( .1,  .05, 4,0, "\\F10 w=0");
		break;
	    default:
		break;
	}

	vp_fat(0);
	vp_penup();

	for( theta = -180.; theta<180.1; theta += 360./ntheta )  {
	    rads = 2. * SF_PI * theta / 360.;
	    cz = cexpf( sf_cmplx(0.,rads) );
#ifdef SF_HAS_COMPLEX_H	 
	    cs = (1.+eps - cz )/2.;
#else
	    cs = sf_crmul(sf_cadd(sf_cmplx(1.+eps,0.),sf_cneg(cz)),0.5);
#endif
	    
	    switch (job) {
		case 0: cs = sf_cmplx( .05, .8*theta/180.); break;
		case 1: break;
#ifdef SF_HAS_COMPLEX_H	 
		case 2: cs *= cs; break;
		case 3: cs = cs*cs + k2; break;
		case 4: cs = csqrtf( cs*cs + k2 ); break;
#else
		case 2: cs = sf_cmul(cs,cs); 
		    break;
		case 3: cs = sf_cadd(sf_cmul(cs,cs),sf_cmplx(k2,0.)); 
		    break;
		case 4: cs = csqrtf(sf_cadd(sf_cmul(cs,cs),sf_cmplx(k2,0.))); 
		    break;
#endif 
		default: break;
	    }

	    x = crealf( cs ); 
	    y = cimagf( cs );

	    vp_upendn ( x, -y);
	}
    }

    return 0;
}

