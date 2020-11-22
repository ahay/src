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
#include "iTracer2D.hh"
#include <list>

int main (int argc, char* argv[]) {
   
// INIT RSF 
    sf_init (argc, argv);
// INPUT FILES
    sf_file xEscFile = sf_input ("in");
	// check that the input is float 
    if ( SF_FLOAT != sf_gettype (xEscFile) ) sf_error ("Need float input: escape-positions file");
    /* escape-positions file */
	sf_file tEscFile=NULL;
    if ( NULL != sf_getstring("esct") ) {
		/* escape-time file */
		tEscFile  = sf_input ("esct");
	}
// OUTPUT FILES
    sf_file xResFile = sf_output ("out");
	/*x-values*/
	sf_file zResFile=NULL;
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

	float x0, z0, p0, dx, dt, sa0;
    if (!sf_getfloat ("x0", &x0)) x0 = 0.f;
	/* x-coordinate of the diffraction point */
    if (!sf_getfloat ("z0", &z0)) z0 = 0.f;
	/* z-coordinate of the diffraction point */
    if (!sf_getfloat ("p0", &p0)) p0 = 0.f;
	/* migration angle */
    if (!sf_getfloat ("sa0", &sa0)) sa0 = 0.f;
	/* scattering-angle */
    if (!sf_getfloat ("dx", &dx)) dx = 5*xStep;
	/* x-range for point detection */
    if (!sf_getfloat ("dt", &dt)) dt = 0.02f;
	/* time-range for point detection */

	// OUTPUT PARAMETERS

  	sf_putint (xResFile, "n1", pNum); 
  	sf_putint (xResFile, "n2", 1); 
  	sf_putint (xResFile, "n3", 1); 
  	sf_putint (xResFile, "n4", 1); 

  	sf_putint (zResFile, "n1", pNum); 
  	sf_putint (zResFile, "n2", 1); 
  	sf_putint (zResFile, "n3", 1); 
  	sf_putint (zResFile, "n4", 1); 

	// READ ESCAPE VOLUMES
	const int volSize = zNum * pNum * xNum;
	float* xVol = sf_floatalloc (volSize);
	sf_seek (xEscFile, 0, SEEK_SET);		
	sf_floatread (xVol, volSize, xEscFile);
	float* tVol = sf_floatalloc (volSize);
	sf_seek (tEscFile, 0, SEEK_SET);		
	sf_floatread (tVol, volSize, tEscFile);

	// IMAGE TRACING
	float* xRes = sf_floatalloc (pNum);
	float* zRes = sf_floatalloc (pNum);
	memset ( xRes, 0, pNum * sizeof (float) );
	memset ( zRes, 0, pNum * sizeof (float) );

	ITracer2D iTracer;
	iTracer.init (zNum, zStart, zStep, 
  			      pNum, pStart, pStep,
			      xNum, xStart, xStep,
				  dx, dt);
	// trace image
	list<float> xpnts;
	list<float> zpnts;

	iTracer.traceImage (xVol, tVol, x0, z0, p0, sa0, &xpnts, &zpnts);
	const int lsize = xpnts.size ();

	std::list<float>::iterator iterx = xpnts.begin ();
	std::list<float>::iterator iterz = zpnts.begin ();

	const int pLim = lsize < pNum ? lsize : pNum;
	for (int ip = 0; ip < pLim; ++ip, ++iterx, ++iterz) {
		xRes[ip] = *iterx;
		zRes[ip] = *iterz;
	}
	
	// WRITE RESULT
    sf_floatwrite (xRes, pNum, xResFile);
    sf_floatwrite (zRes, pNum, zResFile);

	// FINISH
	free (xRes);
	free (zRes);
	free (xVol);
	free (tVol);

    sf_fileclose (xResFile);
	sf_fileclose (zResFile);

	return 0;
}
