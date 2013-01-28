/* Compute reflection event */

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

int main (int argc, char* argv[]) 
{
    sf_file reflFile, derivFile=NULL, dataFile;
    int n1;
    float d1, o1;
    //	
    int tNum, hNum, sNum;
    float tStart, hStart, sStart, tStep, hStep, sStep;
    float eps;
    float vel;
    int dataSize;
    float *depth, *deriv, *data;
    int is, ih, ir, tInd;
    float sPos, offset, z, z1, rx, x, a1, a2, l, dx, way, dt, t, fullWay, bef, aft;

// Initialize RSF 
    sf_init (argc,argv);
// Input files
    reflFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (reflFile) ) sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

    if ( NULL != sf_getstring ("deriv") ) {
	/* first derivative estimated along the reflection boundary */
	derivFile  = sf_input ("deriv");
    } else {
	sf_error ("Need float input: derivatives");
    }
// Output file
    dataFile = sf_output("out");

// Depth/time axis 
    if ( !sf_histint   (reflFile, "n1", &n1) ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (reflFile, "d1", &d1) ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (reflFile, "o1", &o1) ) sf_error ("Need o1= in input");

    if ( !sf_getint ("tn", &tNum) ) tNum = 1001;
    /* number of time samples */
    if ( !sf_getint ("hn", &hNum) ) hNum = 51;
    /* number of offsets */
    if ( !sf_getint ("sn", &sNum) ) sNum = 1;
    /* number of sources */
    if ( !sf_getfloat ("to", &tStart) ) tStart = 0.f;
    /* start time (in s) */
    if ( !sf_getfloat ("ho", &hStart) ) hStart = 0.f;
    /* start offset (in s) */
    if ( !sf_getfloat ("so", &sStart) ) sStart = 0.f;
    /* start source position (in s) */
    if ( !sf_getfloat ("td", &tStep) ) tStep = 0.004f;
    /* step in time (in s) */
    if ( !sf_getfloat ("hd", &hStep) ) hStep = 0.05f;
    /* step in offset (in km) */
    if ( !sf_getfloat ("sd", &sStep) ) sStep= 0.025f;
    /* step in source position (in km) */

    if ( !sf_getfloat ("eps", &eps) ) eps = 0.5 * hStep;
    /* receiver position accuracy (in km) */

    if ( !sf_getfloat ("vel", &vel) ) vel = 2.f;
    /* constant velocity value (in km/s) */

    dataSize = hNum * tNum;

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
    depth = sf_floatalloc (n1);
    sf_seek (reflFile, 0, SEEK_SET);		
    sf_floatread (depth, n1, reflFile);

    deriv = sf_floatalloc (n1);
    sf_seek (derivFile, 0, SEEK_SET);		
    sf_floatread (deriv, n1, derivFile);

    for (is = 0; is < sNum; ++is) {
	sPos = sStart + is * sStep;
	data = sf_floatalloc (dataSize);
	memset (data, 0, dataSize * sizeof (float));
	for (ih = 0; ih < hNum; ++ih) {
	    offset = hStart + ih * hStep;
	    for (ir = 0; ir < n1; ++ir) {
		z  = depth [ir];
		z1 = deriv [ir];
		rx = o1 + ir * d1;
		x = rx - sPos;

		a1 = z - x * z1;
		a2 = z * (1 - z1*z1) - 2*x*z1;

		l = 2 * (x + z*z1)*a1/a2;
		
		if (fabs (l - offset) > eps ) continue;
		
		// time correction
		dx  = sPos + l - rx;
		way = sqrt (dx*dx + z*z);
		dt  = dx * (l - offset) / (vel * way);
		t = 2 * sqrt (x*x + z*z) * a1 / (vel * a2) - dt;
		fullWay = vel * t;

		// put impulse
		tInd = (t - tStart) / tStep;
		if (tInd < 0 || tInd > tNum - 1)
		    continue;

		bef = (t - tInd * tStep) / tStep;
		aft = 1.f - bef;

		data [ih*tNum + tInd]     += aft / fullWay;
		data [ih*tNum + tInd + 1] += bef / fullWay;
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
