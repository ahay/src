/* Linear move-out in the frequency domain
*/
/*
  Copyright (C) 2010 Ioan Vlad

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <rsf.h>
#include <math.h>

int main(int argc, char* argv[])
{
    /* Frequency axis = axis 1 of input, output */ 
    int   nw, iw;
    float dw;

    /* Axis 2 of input */
    int   nx, ix;
    float ox, dx, x;
    float p;    
    
    float *wp = NULL; 
    
    sf_complex *tr = NULL;
    
    sf_file in, out;

    /**********************************************************************/

    /* Read parameters of input file */

    sf_init (argc,argv);

    in  = sf_input ("in" );   
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if(!sf_histint  (in,"n1",&nw)) sf_error ("No n1= in input"); 
    if(!sf_histfloat(in,"d1",&dw)) sf_error ("No d1= in input");
    if(!sf_histint  (in,"n2",&nx)) sf_error ("No n2= in input");
    if(!sf_histfloat(in,"o2",&ox)) sf_error ("No o2= in input");
    if(!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");

    /* Read parameters from the command line */

    if (!sf_getfloat("p", &p)) sf_error("Need p=");
    /* Slope of LMO */

    /* The frequency axis */

    wp = sf_floatalloc( nw );
    for (iw=0; iw < nw; iw++) {
        wp[iw] = p * iw * dw * 2*SF_PI;
    }

    tr = sf_complexalloc( nw );

    for (ix=0; ix < nx; ix++) {

        x = ox + ix * dx;
        
        sf_complexread ( tr, nw, in ); /* Read the gather */

        for (iw=1; iw < nw; iw++) {
#ifdef SF_HAS_COMPLEX_H
            tr[iw] *= cexpf(sf_cmplx(0.,-wp[iw] * x ));
#else
	    tr[iw] = sf_cmul(tr[iw],cexpf(sf_cmplx(0.,-wp[iw] * x )));
#endif
        }
    
        sf_complexwrite( tr, nw , out );

    }
}
