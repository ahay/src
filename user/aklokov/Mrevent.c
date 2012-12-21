/* Compute reflection event
*/

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
    sf_init (argc,argv);
// Input files
    sf_file reflFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (reflFile) ) sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

	sf_file derivFile;
    if ( NULL != sf_getstring ("deriv") ) {
	/* first derivative estimated along the reflection boundary */
		derivFile  = sf_input ("deriv");
	} else {
		sf_error ("Need float input: derivatives");
	}
// Output file
    sf_file dataFile = sf_output("out");

	int n1;
	float d1, o1;
// Depth/time axis 
    if ( !sf_histint   (reflFile, "n1", &n1) ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (reflFile, "d1", &d1) ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (reflFile, "o1", &o1) ) sf_error ("Need o1= in input");

//	
	int tNum, hNum, sNum;
	float tStart, hStart, sStart, tStep, hStep, sStep;

    if ( !sf_getint ("tn", &tNum) ) tNum = 1001;
    if ( !sf_getint ("hn", &hNum) ) hNum = 51;
    if ( !sf_getint ("sn", &sNum) ) sNum = 1;

    if ( !sf_getfloat ("to", &tStart) ) tStart = 0.f;
    if ( !sf_getfloat ("ho", &hStart) ) hStart = 0.f;
    if ( !sf_getfloat ("so", &sStart) ) sStart = 0.f;

    if ( !sf_getfloat ("td", &tStep) ) tStep = 0.004f;
    if ( !sf_getfloat ("hd", &hStep) ) hStep = 0.05f;
    if ( !sf_getfloat ("sd", &sStep) ) sStep= 0.025f;

	float eps;
    if ( !sf_getfloat ("eps", &eps) ) eps = 0.01f;

	float vel;
    if ( !sf_getfloat ("vel", &vel) ) vel = 2.f;


	const int dataSize = hNum * tNum;

   	sf_putint   (dataFile, "n1", tNum);
   	sf_putfloat (dataFile, "d1", tStep);
   	sf_putfloat (dataFile, "o1", tStart);

   	sf_putint   (dataFile, "n2", hNum);
   	sf_putfloat (dataFile, "d2", hStep);
   	sf_putfloat (dataFile, "o2", hStart);

   	sf_putint   (dataFile, "n3", sNum);
   	sf_putfloat (dataFile, "d3", sStep);
   	sf_putfloat (dataFile, "o3", sStart);


// read data
	float* depth = sf_floatalloc (n1);
	sf_seek (reflFile, 0, SEEK_SET);		
	sf_floatread (depth, n1, reflFile);

	float* deriv = sf_floatalloc (n1);
	sf_seek (derivFile, 0, SEEK_SET);		
	sf_floatread (deriv, n1, derivFile);

	for (int is = 0; is < sNum; ++is) {
		const float sPos = sStart + is * sStep;
		float* data = sf_floatalloc (dataSize);
		memset (data, 0, dataSize * sizeof (float));
		for (int ih = 0; ih < hNum; ++ih) {
			const float offset = hStart + ih * hStep;
			for (int ir = 0; ir < n1; ++ir) {
				const float z  = depth [ir];
				const float z1 = deriv [ir];
				const float x = o1 + ir * d1 - sPos;

				const float a1 = z - x * z1;
				const float a2 = z * (1 - z1*z1) - 2*x*z1;

				const float l = 2 * (x + z*z1)*a1/a2;
		
				if (fabs (l - offset) > eps ) continue;
		
				const float t = 2 * sqrt (x*x + z*z) * a1 / (vel * a2);

				// put impulse
				const int tInd = (t - tStart) / tStep;
				if (tInd < 0 || tInd > tNum - 1)
					continue;

				const float bef = (t - tInd * tStep) / tStep;
				const float aft = 1.f - bef;

				data [ih*tNum + tInd]     += aft;
				data [ih*tNum + tInd + 1] += bef;
			}
		}
		sf_floatwrite (data, dataSize, dataFile);
		free (data);
	}

	sf_fileclose (reflFile);
	sf_fileclose (derivFile);
	sf_fileclose (dataFile);

	free (depth);
	free (deriv);

	return 0;
}
