/* Convert acoustic impedance to reflectivity. 

August 2013 program of the month:
http://www.ahay.org/rsflog/index.php?/archives/350-Program-of-the-month-sfai2refl.html
*/
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
    float imp1, imp2, *imp, *sig;
    sf_file ai, mod;

    sf_init (argc,argv);
    ai  = sf_input("in");
    mod = sf_output("out");

    if (!sf_histint(ai,"n1",&nt)) sf_error("No n1= in input");
    n2 = sf_leftsize(ai,1);

    imp = sf_floatalloc (nt);
    sig = sf_floatalloc (nt);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(imp,nt,ai);

	imp1=imp[0];
	for (it=0; it < nt-1; it++) {
	    imp2 = imp[it+1];
	    sig[it] = (imp2-imp1)/(imp2+imp1+FLT_EPSILON);
	    imp1 = imp2;
	}
	sig[nt-1] = 0.;

	sf_floatwrite(sig,nt,mod);
    }


    exit (0);
}

/* 	$Id: Mai2refl.c 12812 2014-06-09 19:51:49Z sfomel $	 */
