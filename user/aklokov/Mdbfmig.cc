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
#include "depthBackfitMigrator2D.hh"

int main (int argc, char* argv[]) {
   
// INIT RSF 

    sf_init (argc, argv);

// INPUT FILES

    sf_file piFile = sf_input ("in");
	// check that the input is float 
    if ( SF_FLOAT != sf_gettype (piFile) ) sf_error ("Need float input: partial-images file");
    /* partial-images file */

	sf_file xEscFile;
    if ( NULL != sf_getstring("escx") ) {
		/* escape-positions file */
		xEscFile  = sf_input ("escx");
	}

	sf_file tEscFile;
    if ( NULL != sf_getstring("esct") ) {
		/* escape-time file */
		tEscFile  = sf_input ("esct");
	}

// OUTPUT FILES
    sf_file resFile = sf_output ("out");
	/* backfit-migrated images */

	// INPUT PARAMETERS

// escape volumes dimensions
	int zNum; float zStart; float zStep;
	int pNum; float pStart; float pStep;
	int xNum; float xStart; float xStep;
	int rNum; float rStart; float rStep;
// depth axis 
    if ( !sf_histint   (piFile, "n1", &zNum) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (piFile, "d1", &zStep) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (piFile, "o1", &zStart) ) sf_error ("Need o1= in input");
// x-axis 
    if ( !sf_histint   (piFile, "n2", &xNum) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (piFile, "d2", &xStep) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (piFile, "o2", &xStart) ) sf_error ("Need o2= in input");
// dip axis
    if ( !sf_histint   (piFile, "n3", &pNum) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (piFile, "d3", &pStep) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (piFile, "o3", &pStart) )   sf_error ("Need o3= in input");
// r-axis
    if ( !sf_histint   (xEscFile, "n2", &rNum) )     sf_error ("Need n2= in input");
    if ( !sf_histfloat (xEscFile, "d2", &rStep) )    sf_error ("Need d2= in input");
    if ( !sf_histfloat (xEscFile, "o2", &rStart) )   sf_error ("Need o2= in input");

	// OUTPUT PARAMETERS

	// MEMORY ALLOCATION

	const int piSize  = zNum * pNum;
	const int volSize = zNum * rNum * xNum;

	float* tVol = sf_floatalloc (volSize);
	float* xVol = sf_floatalloc (volSize);

	float* piData  = sf_floatalloc (piSize);
	float* piImage = sf_floatalloc (piSize);

	// READ ESCAPE VOLUMES

	sf_seek (xEscFile, 0, SEEK_SET);		
	sf_floatread (xVol, volSize, xEscFile);
	sf_seek (tEscFile, 0, SEEK_SET);		
	sf_floatread (tVol, volSize, tEscFile);

	// MAIN LOOP

	DepthBackfitMigrator2D dbfmig;
	dbfmig.init (zNum, zStart, zStep, 
  	 	         pNum, pStart, pStep,
			     xNum, xStart, xStep,
			     rNum, rStart, rStep);

	for (int ip = 0; ip < pNum; ++ip) {
		const float curP = pStart + ip * pStep;			

		memset ( piData,  0, piSize * sizeof (int) );
		memset ( piImage, 0, piSize * sizeof (int) );

		// read partial image
		const size_t startPos = ip * piSize * sizeof(float);
	    sf_seek (piFile, startPos, SEEK_SET);
	    sf_floatread (piData, piSize, piFile);

		// get image
		dbfmig.processPartialImage (piData, curP, xVol, tVol, piImage);

		// write result
	 	sf_seek (resFile, startPos, SEEK_SET);
	    sf_floatwrite (piImage, piSize, resFile);
	}

	// FINISH

	free (xVol);
	free (tVol);
	
	free (piData);
	free (piImage);

    sf_fileclose (resFile);
	sf_fileclose (piFile);
    sf_fileclose (tEscFile);
	sf_fileclose (xEscFile);

	return 0;
}
