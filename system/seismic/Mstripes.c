/* Model the positions and dips of the constant offset, source, midpoint, and receiver strikes in a source vs. offset space. */
/*
  Copyright (C) 1994 The Board of Trustees of Stanford University
  
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

#include <assert.h>
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])

{
    float *tabout                            ;
    float  oh, os, oy, or, dh, ds, dy, dr    ;
    float  sloc, hloc, yloc, rloc            ;
    int    nh, ns, ny, nr                    ;
    int    mo, ms, my, mr                    ;
    int    i, j, k, l                        ;
    sf_file inp, out;

    sf_init( argc, argv ) ;
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp, "n1"   , &nh)) sf_error("No n1= in input");
    if (!sf_histint(inp, "n2"   , &ns)) sf_error("No n2= in input");
    if (!sf_histfloat(inp, "o1" , &oh)) sf_error("No o1= in input");
    if (!sf_histfloat(inp, "o2" , &os)) sf_error("No o2= in input");
    if (!sf_histfloat(inp, "d1" , &dh)) sf_error("No d1= in input");
    if (!sf_histfloat(inp, "d2" , &ds)) sf_error("No d2= in input");
		    
    if (!sf_getint("mo", &mo)) sf_error("Need mo=");
    /* offset parameter, a constant offset line will appear in the output every o offset */
    if (!sf_getint("ms", &ms)) sf_error("Need ms=");
    /* source parameter, a constant source line will appear in the output every s source */
    if (!sf_getint("my", &my)) sf_error("Need my=");
    /* midpoint parameter, a constant midpoint line will appear in the output every y midpoint */
    if (!sf_getint("mr", &mr)) sf_error("Need mr=");
    /* receiver parameter, a constant receiver line will appear in the output every r receiver */

    ny = 2 * (ns - 1) + nh ;
    oy = oh/2              ;
    dy = dh/2              ;
    
    nr = ns - 1 + nh ;
    or = os + oh     ;
    dr = dh          ;

    tabout = sf_floatalloc( nh * ns);

    for ( k = 0 ; k < nh * ns ; k++) 
	tabout[k] = 1.0f;

    for ( i = 0 ; i < nh ; i ++ ) {
	for ( j = 0 ; j < ns ; j ++ ) {
	    hloc = oh + i * dh   ;
	    sloc = os + j * ds   ;
	    yloc = sloc + hloc/2 ;
	    rloc = sloc + hloc   ;

	    k = (int) floorf( (yloc - oy)/dy + 0.5 ) ;
	    l = (int) floorf( (rloc - or)/dr + 0.5 ) ;

	    assert( k >= 0 && k < ny ) ;
	    assert( l >= 0 && l < nr ) ;
	    
	    if ( !(i % mo) || !(j % ms) || !(k % my) || !(l % mr) ) {
		tabout[nh * j + i] = 0.0 ;
	    }
	}
    }

    sf_floatwrite(tabout, nh * ns, out);
    
    exit(0);
}

