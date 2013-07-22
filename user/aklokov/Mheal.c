/* Heal data */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

int main (int argc, char* argv[]) {

    int n1, n2;
    int i1, i2;

// Initialize RSF 
    sf_init (argc,argv);

// Input files
    sf_file inFile = sf_input ("in");
// Output file
    sf_file outFile = sf_output ("out");

    // data parameters
    if ( !sf_histint   (inFile, "n1", &n1) )  sf_error ("Need n1= in input");
    if ( !sf_histint   (inFile, "n2", &n2) )  sf_error ("Need n2= in input");

    float* trace = sf_floatalloc (n1);	

    sf_floatread  (trace, n1, inFile);		
    sf_floatwrite (trace, n1, outFile);

    for (i2 = 1; i2 < n2 - 1; ++i2) {
	    sf_floatread (trace, n1, inFile);		
		for (i1 = 1; i1 < n1 - 1; ++i1) {
			if ( fabs(trace[i1]) > 1e-6)
				continue; // good sample
			if ( fabs(trace[i1-1]) < 1e-6 || fabs(trace[i1+1]) < 1e-6) 
				continue; // too big gap
			trace[i1] = 0.5 * (trace[i1-1] + trace[i1+1]);
		}
	    sf_floatwrite (trace, n1, outFile);
	}

    sf_floatread  (trace, n1, inFile);		
    sf_floatwrite (trace, n1, outFile);

    sf_fileclose (inFile);
    sf_fileclose (outFile);

    free (trace);

    return 0;
}
