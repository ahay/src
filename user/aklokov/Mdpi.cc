/* */

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

#include <rsf.hh>
#include <float.h>
#include "support.hh"

int pind;

bool isPointInsideTriangle (float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {

	float x01 = x0 - x1;
	float y01 = y0 - y1;
	if (fabs (x01) < 1e-6 && fabs (y01) < 1e-6) return true;

	float x02 = x0 - x2;
	float y02 = y0 - y2;
	if (fabs (x02) < 1e-6 && fabs (y02) < 1e-6) return true;

	float x03 = x0 - x3;
	float y03 = y0 - y3;
	if (fabs (x03) < 1e-6 && fabs (y03) < 1e-6) return true;

	//
		
	float x31 = x3 - x1; 
	float y31 = y3 - y1;

	float x21 = x2 - x1;
	float y21 = y2 - y1;

	// dot products

	float dot00 = x31*x31 + y31*y31;
	float dot01 = x31*x21 + y31*y21;
	float dot02 = x31*x01 + y31*y01;
	float dot11 = x21*x21 + y21*y21;
	float dot12 = x21*x01 + y21*y01;

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1);
}

int main (int argc, char* argv[]) {
   
// INIT RSF 
    sf_init (argc, argv);
// INPUT FILES
    sf_file xEscFile = sf_input ("in");
	// check that the input is float 
    if ( SF_FLOAT != sf_gettype (xEscFile) ) sf_error ("Need float input: escape-positions file");
    /* escape-positions file */
	sf_file tEscFile;
    if ( NULL != sf_getstring("esct") ) {
		/* escape-time file */
		tEscFile  = sf_input ("esct");
	}
// OUTPUT FILES
    sf_file xResFile = sf_output ("out");
	/*x-values*/
	sf_file zResFile;
    if ( NULL != sf_getstring ("zres") ) {
		zResFile  = sf_output ("zres");
    }

	// INPUT PARAMETERS

// escape volumes dimensions
	int zNum; float zStart; float zStep;
	int pNum; float pStart; float pStep;
	int xNum; float xStart; float xStep;
// depth axis 
    if ( !sf_histint   (xEscFile, "n1", &zNum) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (xEscFile, "d1", &zStep) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (xEscFile, "o1", &zStart) ) sf_error ("Need o1= in input");
// dip axis
    if ( !sf_histint   (xEscFile, "n2", &pNum) )     sf_error ("Need n2= in input");
    if ( !sf_histfloat (xEscFile, "d2", &pStep) )    sf_error ("Need d2= in input");
    if ( !sf_histfloat (xEscFile, "o2", &pStart) )   sf_error ("Need o2= in input");
// x-axis 
    if ( !sf_histint   (xEscFile, "n3", &xNum) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (xEscFile, "d3", &xStep) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (xEscFile, "o3", &xStart) ) sf_error ("Need o3= in input");

	float x0, z0, p0;
    if (!sf_getfloat ("x0", &x0)) x0 = 0.f;
	/* x-coordinate of the diffraction point */
    if (!sf_getfloat ("z0", &z0)) z0 = 0.f;
	/* z-coordinate of the diffraction point */
    if (!sf_getfloat ("p0", &p0)) p0 = 0.f;
	/* migration angle */

	const int zInd = (z0 - zStart) / zStep;
	const int pInd = (p0 - pStart) / pStep;
	const int xInd = (x0 - xStart) / xStep;

	// OUTPUT PARAMETERS

  	sf_putint (xResFile, "n1", pNum); 
  	sf_putint (xResFile, "n2", 1); 
  	sf_putint (xResFile, "n3", 1); 
  	sf_putint (xResFile, "n4", 1); 

  	sf_putint (zResFile, "n1", pNum); 
  	sf_putint (zResFile, "n2", 1); 
  	sf_putint (zResFile, "n3", 1); 
  	sf_putint (zResFile, "n4", 1); 

	// 

	// read escape volumes
	const int volSize = zNum * pNum * xNum;
	float* xVol = sf_floatalloc (volSize);
	sf_seek (xEscFile, 0, SEEK_SET);		
	sf_floatread (xVol, volSize, xEscFile);
	float* tVol = sf_floatalloc (volSize);
	sf_seek (tEscFile, 0, SEEK_SET);		
	sf_floatread (tVol, volSize, tEscFile);

	// escape line for the diffraction point
	float* et = sf_floatalloc (pNum);
	memset ( et, 0, pNum * sizeof (float) );
	float* ex = sf_floatalloc (pNum);
	memset ( ex, 0, pNum * sizeof (float) );

	for (int ip = 0; ip < pNum; ++ip) {
		pind = xInd * pNum * zNum + ip * zNum + zInd;
		ex[ip] = xVol [pind];
		et[ip] = tVol [pind];
	}

	// constant-dip panel extraction
	const int panelSize = zNum * xNum;	
	float* xPanel = sf_floatalloc (panelSize);
	float* tPanel = sf_floatalloc (panelSize);
	memset ( xPanel, 0, panelSize * sizeof (float) );
	memset ( tPanel, 0, panelSize * sizeof (float) );

	float* pXPanel = xPanel;
	float* pTPanel = tPanel;
	for (int ix = 0; ix < xNum; ++ix) {
		for (int iz = 0; iz < zNum; ++iz, ++pXPanel, ++pTPanel) {
			pind = ix * pNum * zNum + pInd * zNum + iz;
			*pXPanel = xVol [pind];
			*pTPanel = tVol [pind];
		}
	}

	// line-in-depth construction
	float* xRes = sf_floatalloc (pNum);
	float* zRes = sf_floatalloc (pNum);
	memset ( xRes, 0, pNum * sizeof (float) );
	memset ( zRes, 0, pNum * sizeof (float) );

	// loop over image points

	for (int ix = 0; ix < xNum-1; ++ix) {
		for (int iz = 0; iz < zNum-1; ++iz) {

			// upper triangle

			Point2D dataPoint [3];	
			int ind0 = (ix + 1) * zNum + iz;
			dataPoint[0].setX (xPanel [ind0]);
			dataPoint[0].setY (tPanel [ind0]);
			int ind1 = ix * zNum + iz + 1;
			dataPoint[1].setX (xPanel [ind1]);
			dataPoint[1].setY (tPanel [ind1]);
			int ind2 = ix * zNum + iz;
			dataPoint[2].setX (xPanel [ind2]);
			dataPoint[2].setY (tPanel [ind2]);

			// go along data line

			for (int ip = 0; ip < pNum; ++ip) {

				float curT = et[ip];
				float curX = ex[ip];

				Point2D p0 (curX, curT);

				if ( isPointInsideTriangle (curX, curT,
											dataPoint[0].getX (), dataPoint[0].getY (), 
											dataPoint[1].getX (), dataPoint[1].getY (),
											dataPoint[2].getX (), dataPoint[2].getY ()) ) {
		
					const float x13 = dataPoint[0].getX () - dataPoint[2].getX ();
					const float x32 = dataPoint[2].getX () - dataPoint[1].getX ();
	
					const float y23 = dataPoint[1].getY () - dataPoint[2].getY ();
					const float y13 = dataPoint[0].getY () - dataPoint[2].getY ();

					const float denom = y23 * x13 + x32 * y13; 

					float w0 = ( y23 * (curX - dataPoint[2].getX())  + x32 * (curT - dataPoint[2].getY()) ) / denom;
					float w1 = ( -y13 * (curX - dataPoint[2].getX()) + x13 * (curT - dataPoint[2].getY()) ) / denom;
					float w2 = 1 - w0 - w1;

					Point2D imagePoint [3];	
					imagePoint[0].setX (xStart + (ix + 1) * xStep);
					imagePoint[0].setY (zStart + iz * zStep);
					imagePoint[1].setX (xStart + ix * xStep);
					imagePoint[1].setY (zStart + (iz + 1) * zStep);
					imagePoint[2].setX (xStart + ix * xStep);
					imagePoint[2].setY (zStart + iz * zStep);


					float x = w0 * imagePoint[0].getX () + w1 * imagePoint[1].getX () + w2 * imagePoint[2].getX ();
					float z = w0 * imagePoint[0].getY () + w1 * imagePoint[1].getY () + w2 * imagePoint[2].getY ();						

					xRes [ip] = x;
					zRes [ip] = z;
				}
			}
	
			// lower triangle

			ind0 = (ix + 1) * zNum + iz;
			dataPoint[0].setX (xPanel [ind0]);
			dataPoint[0].setY (tPanel [ind0]);
			ind1 = ix * zNum + iz + 1;
			dataPoint[1].setX (xPanel [ind1]);
			dataPoint[1].setY (tPanel [ind1]);
			ind2 = (ix + 1) * zNum + iz + 1;
			dataPoint[2].setX (xPanel [ind2]);
			dataPoint[2].setY (tPanel [ind2]);

			// go along data line

			for (int ip = 0; ip < pNum; ++ip) {
				float curT = et[ip];
				float curX = ex[ip];

				Point2D p0 (curX, curT);

				if ( isPointInsideTriangle (curX, curT,
											dataPoint[0].getX (), dataPoint[0].getY (), 
											dataPoint[1].getX (), dataPoint[1].getY (),
											dataPoint[2].getX (), dataPoint[2].getY ()) ) {
					
					const float x13 = dataPoint[0].getX () - dataPoint[2].getX ();
					const float x32 = dataPoint[2].getX () - dataPoint[1].getX ();
	
					const float y23 = dataPoint[1].getY () - dataPoint[2].getY ();
					const float y13 = dataPoint[0].getY () - dataPoint[2].getY ();

					const float denom = y23 * x13 + x32 * y13; 

					float w0 = ( y23 * (curX - dataPoint[2].getX())  + x32 * (curT - dataPoint[2].getY()) ) / denom;
					float w1 = ( -y13 * (curX - dataPoint[2].getX()) + x13 * (curT - dataPoint[2].getY()) ) / denom;
					float w2 = 1 - w0 - w1;

					Point2D imagePoint [3];	
					imagePoint[2].setX (xStart + ix * xStep);
					imagePoint[2].setY (zStart + iz * zStep);
					imagePoint[0].setX (xStart + (ix + 1) * xStep);
					imagePoint[0].setY (zStart + iz * zStep);
					imagePoint[1].setX (xStart + ix * xStep);
					imagePoint[1].setY (zStart + (iz + 1) * zStep);

					float x = w0 * imagePoint[0].getX () + w1 * imagePoint[1].getX() + w2 * imagePoint[2].getX();
					float z = w0 * imagePoint[0].getY () + w1 * imagePoint[1].getY() + w2 * imagePoint[2].getY();						

					xRes [ip] = x;
					zRes [ip] = z;
				}
			}
		}
	}

	// finish 

    sf_floatwrite (xRes, pNum, xResFile);
    sf_floatwrite (zRes, pNum, zResFile);

	free (ex);
	free (et);
	free (xPanel);
	free (tPanel);	
	free (xRes);
	free (zRes);

    sf_fileclose (xResFile);
	sf_fileclose (zResFile);

	return 0;
}
