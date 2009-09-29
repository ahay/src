/* Apply dip correction for angle-gathers computed with absolute offset */
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

/* 

tan t_0 = k_|h|/k_z

tan t = tan t_0 * 1/sqrt{ 1 + (dz/dx)^2 }

*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    bool inv, verb;
    int nz, na, iz, ia, nw, n3, i3;
    float **gather, *trace, *modl, *coord, *dzdx; 
    float da, a0, t,d;
    sf_file in, out, dip;

    sf_init (argc,argv);

    if (!sf_getbool("verb",&verb)) verb=false;

    in   = sf_input (  "in");
    out  = sf_output( "out");
    dip  = sf_input ( "dip");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint  (in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint  (in,"n2",&na)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&da)) sf_error("No d2= in input"); 
    if (!sf_histfloat(in,"o2",&a0)) sf_error("No o2= in input"); 

    if (!sf_getint("nw",&nw)) nw=4; /* accuracy level */

    sf_prefilter_init (nw,na,2*na);

    n3 = sf_leftsize(in,2);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    gather = sf_floatalloc2(nz,na);
    trace  = sf_floatalloc (   na);
    coord  = sf_floatalloc (   na);
    modl   = sf_floatalloc (   na);
    dzdx   = sf_floatalloc (nz   );

    for (i3=0; i3 < n3; i3++) { /* loop over CIG */
	if(verb) sf_warning("%d of %d",i3+1,n3);

	sf_floatread( dzdx    ,nz   ,dip );
	sf_floatread(gather[0],nz*na,in  );

	for (iz=0; iz < nz; iz++) { /* loop over depth */
	    d =  dzdx[iz];
	    
	    for (ia=0; ia < na; ia++) { /* loop over tan */
		t = a0+ia*da;

		/* 
		   loop over output t;
		   find map to input t0, i.e. t0=f(t);
		*/
		coord[ia] = t * sqrtf(1.+d*d);

		trace[ia] = gather[ia][iz];
	    }

	    sf_prefilter_apply (na,trace);
	    sf_int1_init (coord, a0, da, na, sf_spline_int, nw, na);
	    sf_int1_lop (false,false,na,na,trace,modl);
	    
	    for (ia=0; ia < na; ia++) {
		gather[ia][iz] = modl[ia];
	    }
	}
	sf_floatwrite(gather[0],nz*na,out);
    }

    exit(0);
}

