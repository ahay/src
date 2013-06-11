/* 3D volume rotation about the vertical axes */
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

int n1, n2, n3;
float d1, d2, d3, o1, o2, o3;
float* trace;
sf_file inFile, outFile;

bool isPointInsideTriangle (float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {

	double x01 = x0 - x1;
	double y01 = y0 - y1;
	if (fabs (x01) < 1e-6 && fabs (y01) < 1e-6) return true;

	double x02 = x0 - x2;
	double y02 = y0 - y2;
	if (fabs (x02) < 1e-6 && fabs (y02) < 1e-6) return true;

	double x03 = x0 - x3;
	double y03 = y0 - y3;
	if (fabs (x03) < 1e-6 && fabs (y03) < 1e-6) return true;

	//
		
	double x31 = x3 - x1; 
	double y31 = y3 - y1;

	double x21 = x2 - x1;
	double y21 = y2 - y1;

	// dot products

	double dot00 = x31*x31 + y31*y31;
	double dot01 = x31*x21 + y31*y21;
	double dot02 = x31*x01 + y31*y01;
	double dot11 = x21*x21 + y21*y21;
	double dot12 = x21*x01 + y21*y01;

	// Compute barycentric coordinates
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1);
}

void getTrace (float xr, float yr, float* restrace) {

	int ind1x = xr / d2;
	int ind1y = yr / d3;
	int ind2x = ind1x + 1;
	int ind2y = ind1y;
	int ind3x = ind1x;
	int ind3y = ind1y + 1;

	float x1 = ind1x * d2;
	float y1 = ind1y * d3;
	float x2 = (ind1x + 1) * d2;
	float y2 = ind1y * d3;
	float x3 = ind1x * d2;
	float y3 = (ind1y + 1) * d3;

	// check triangle
	bool res = isPointInsideTriangle (xr, yr, x1, y1, x2, y2, x3, y3);
	if (!res) {
		ind1x += 1;
		ind1y += 1;
		x1 = ind1x * d2;
		y1 = ind1y * d3;
	}

	const float x21 = x2 - x1;
	const float x31 = x3 - x1;
	const float y21 = y2 - y1;
	const float y31 = y3 - y1;

	const float denom3 = y31*x21 - x31*y21;
	if (fabs (denom3) < 1e-9) return;
	const float denom2 = y21*x31 - x21*y31;
	if (fabs (denom2) < 1e-9) return;
	
	const float m3 = ( (yr - y1)*x21 - (xr - x1)*y21 ) / denom3;
	const float m2 = ( (yr - y1)*x31 - (xr - x1)*y31 ) / denom2;
	const float m1 = 1.f - m2 - m3;
	
	// interpolation
	float* trace1 = sf_floatalloc (n1);	
	float* trace2 = sf_floatalloc (n1);	
	float* trace3 = sf_floatalloc (n1);	

	int xshift = o2 / d2;
	int yshift = o3 / d3;

	ind1x -= xshift;
	ind2x -= xshift;
	ind3x -= xshift;
	ind1y -= yshift;
	ind2y -= yshift;
	ind3y -= yshift;

	// trace1 
	size_t posr = (ind1y * n2 + ind1x) * n1 * sizeof (float);
	sf_seek (inFile, posr, SEEK_SET);		
	sf_floatread (trace1, n1, inFile);	
	// trace2
	posr = (ind2y * n2 + ind2x) * n1 * sizeof (float);
	sf_seek (inFile, posr, SEEK_SET);		
	sf_floatread (trace2, n1, inFile);	
	// trace3 
	posr = (ind3y * n2 + ind3x) * n1 * sizeof (float);
	sf_seek (inFile, posr, SEEK_SET);		
	sf_floatread (trace3, n1, inFile);	

	float m = m1 + m2 + m3;

	for (int i = 0; i < n1; ++i) {
		restrace[i] = m1 * trace1[i] + m2 * trace2[i] + m3 * trace3[i];
		
	}

	free (trace1);
	free (trace2);
	free (trace3);

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
    if ( !sf_getfloat ("xf", &xf) ) xf = 0.f;
    /* inline-coord of the vertical axes */
    if ( !sf_getfloat ("yf", &yf) ) yf = 0.f;
    /* xline-coord of the vertical axes */
	if ( !sf_getfloat ("theta", &theta) ) theta = 0.f;	
    /* rotation angle */	

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
				getTrace (xr, yr, trace);
			
			sf_floatwrite (trace, n1, outFile);					
		}
	}

	free (trace);

	sf_fileclose (inFile);
    sf_fileclose (outFile);

    return 0;
}
