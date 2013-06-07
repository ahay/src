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

int n1, n2, n3;
float d1, d2, d3, o1, o2, o3;
float* trace;
sf_file inFile, outFile;
/*void getTriangle (float x1, float y1, 	
				  float* px1, float* py1,
				  float* px2, float* py2,
				  float* px3, float* py3) {
	

}*/

void getTrace (float xr, float yr) {

	const int i2r = (xr - o2) / d2;
	const int i3r = (yr - o3) / d3;

	const size_t posr = (i3r * n2 + i2r) * n1 * sizeof (float);
	sf_seek (inFile, posr, SEEK_SET);		
    sf_floatread (trace, n1, inFile);	
	
	return;
}


int main (int argc, char* argv[]) {

// Initialize RSF 
    sf_init (argc, argv);
// Input files
    inFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (inFile) ) sf_error ("Need float input");
    /* input volume */
// Output file
    outFile = sf_output("out");
    /* output volume */

// Depth/time axis 

    if ( !sf_histint   (inFile, "n1", &n1) ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (inFile, "d1", &d1) ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (inFile, "o1", &o1) ) sf_error ("Need o1= in input");
    if ( !sf_histint   (inFile, "n2", &n2) ) sf_error ("Need n2= in input");
    if ( !sf_histfloat (inFile, "d2", &d2) ) sf_error ("Need d2= in input");
    if ( !sf_histfloat (inFile, "o2", &o2) ) sf_error ("Need o2= in input");
    if ( !sf_histint   (inFile, "n3", &n3) ) sf_error ("Need n3= in input");
    if ( !sf_histfloat (inFile, "d3", &d3) ) sf_error ("Need d3= in input");
    if ( !sf_histfloat (inFile, "o3", &o3) ) sf_error ("Need o3= in input");


	float xf, yf, theta;			
    /* x-coord of the fixed point */
    if ( !sf_getfloat ("xf", &xf) ) xf = 0.f;
    /* z-coord of the fixed point */
    if ( !sf_getfloat ("yf", &yf) ) yf = 0.f;
    /* rotation angle */	
	if ( !sf_getfloat ("theta", &theta) ) theta = 0.f;	

	const float thetaRad = theta * M_PI / 180.f;
	const float sint = sin (-thetaRad);
	const float cost = cos (-thetaRad);

	// limits
	const float yMin = o3;
	const float yMax = o3 + (n3-1)*d3;
	const float xMin = o2;
	const float xMax = o2 + (n2-1)*d2;

    trace  = sf_floatalloc (n1);	

	for (int i3 = 0; i3 < n3; ++i3) {
		float y = o3 + i3 * d3 - yf;
		for (int i2 = 0; i2 < n2; ++i2) {
			float x = o2 + i2 * d2 - xf;

			for (int it = 0; it < n1; ++it)
				trace[it] = 0.f;

			float xr = cost * x - sint * y + xf;
			float yr = sint * x + cost * y + yf;

			if (yr > yMin && yr < yMax && xr > xMin && xr < xMax) 
				getTrace (xr, yr);
			
			sf_floatwrite (trace, n1, outFile);					
		}
	}

	sf_fileclose (inFile);
    sf_fileclose (outFile);

    return 0;
}
