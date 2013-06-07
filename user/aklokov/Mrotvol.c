/* Rotate data */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

// Initialize RSF 
    sf_init (argc, argv);
// Input files
    sf_file inFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (inFile) ) sf_error ("Need float input");
    /* input volume */
// Output file
    sf_file outFile = sf_output("out");
    /* output volume */

// Depth/time axis 
	int n1, n2, n3;
	float d1, d2, d3, o1, o2, o3;
    if ( !sf_histint   (inFile, "n1", &n1) ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (inFile, "d1", &d1) ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (inFile, "o1", &o1) ) sf_error ("Need o1= in input");
    if ( !sf_histint   (inFile, "n2", &n2) ) sf_error ("Need n2= in input");
    if ( !sf_histfloat (inFile, "d2", &d2) ) sf_error ("Need d2= in input");
    if ( !sf_histfloat (inFile, "o2", &o2) ) sf_error ("Need o2= in input");
    if ( !sf_histint   (inFile, "n3", &n3) ) sf_error ("Need n3= in input");
    if ( !sf_histfloat (inFile, "d3", &d3) ) sf_error ("Need d3= in input");
    if ( !sf_histfloat (inFile, "o3", &o3) ) sf_error ("Need o3= in input");


	float xf, zf, theta;			
    /* x-coord of the fixed point */
    if ( !sf_getfloat ("xf", &xf) ) xf = 0.f;
    /* z-coord of the fixed point */
    if ( !sf_getfloat ("zf", &zf) ) zf = 0.f;
    /* rotation angle */	
	if ( !sf_getfloat ("theta", &theta) ) theta = 0.f;	

	const float thetaRad = theta * M_PI / 180.f;
	const float sint = sin (-thetaRad);
	const float cost = cos (-thetaRad);

	const float zLim = n1 - 1;
	const float xLim = n2 - 1;

	const int sliceSize = n1 * n2;
    float* sliceIn  = sf_floatalloc (sliceSize);	
    float* sliceOut = sf_floatalloc (sliceSize);	

	for (int i3 = 0; i3 < n3; ++i3) {
		// read data
	    sf_floatread (sliceIn, sliceSize, inFile);					
	
		for (int i2 = 0; i2 < n2; ++i2) {
			float x = o2 + i2 * d2 - xf;
			for (int i1 = 0; i1 < n1; ++i1) {
				float z = o1 + i1 * d1 - zf;

				float zr = cost * z - sint * x + zf;
				const int i1r = (zr - o1) / d1;
				if (i1r < 0 || i1r > zLim) continue;

				float xr = sint * z + cost * x + xf;
				const int i2r = (xr - o2) / d2;
				if (i1r < 0 || i1r > xLim) continue;

				const int indIn  = i1r + i2r * n1;			
				const int indOut = i1 + i2 * n1;			
				
				sliceOut [indOut] = sliceIn [indIn];
			}
		}
	    
		sf_floatwrite (sliceOut, sliceSize, outFile);					
	}

	sf_fileclose (inFile);
    sf_fileclose (outFile);

    free (sliceIn);
    free (sliceOut);

    return 0;
}
