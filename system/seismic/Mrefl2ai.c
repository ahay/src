/* Convert reflectivity to acoustic impedance. */
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
#include <float.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int nt, it, n2, i2;
    float r, a, *imp, *sig;
    sf_file ref, ai, a0;

    sf_init (argc,argv);
    ref  = sf_input("in");
    ai = sf_output("out");
    a0 = sf_input("a0"); /* impedance on the surface */

    if (!sf_histint(ref,"n1",&nt)) sf_error("No n1= in input");
    n2 = sf_leftsize(ref,1);

    imp = sf_floatalloc (nt);
    sig = sf_floatalloc (nt);
    
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(&a,1,a0);
	sf_floatread(sig,nt,ref);

	for (it=0; it < nt; it++) {
	    imp[it] = a;
	    r = sig[it];
	    a *= (1.0f+r)/(1.0f-r);
	}
	
	sf_floatwrite(imp,nt,ai);
    }
    
    exit (0);
}


