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
#include <time.h>
int main (int argc, char* argv[]) {
   
// INIT RSF 

    sf_init (argc, argv);

// INPUT FILES

    sf_file piFile = sf_input ("in");
    // check that the input is float 
    if ( SF_FLOAT != sf_gettype (piFile) ) sf_error ("Need float input: partial-images file");
    /* partial-images file */

    sf_file xEscFile = NULL;
    if ( NULL != sf_getstring("escx") ) {
	/* escape-positions file */
	xEscFile  = sf_input ("escx");
    }

    sf_file tEscFile = NULL;
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

    bool isAA;
    float dx, dt, xlim, xapert;
    int ppn; float ppo, ppd;
    int izn, ixn; float izo, ixo, izd, ixd;
    int pj;
    int sNum; float sStart; float sStep;

    if (!sf_getint ("ppn", &ppn)) ppn = pNum;
    /* number of processed partial images */
    if (!sf_getfloat ("ppo", &ppo)) ppo = pStart;
    /* first processed partial image */
    if (!sf_getfloat ("ppd", &ppd)) ppd = pStep;
    /* step in processed partial images */

    // IMAGE PARAMS
    if (!sf_getint ("izn", &izn))        izn = zNum;	
    /* number of imaged depth samples */
    if (!sf_getint ("ixn", &ixn))        ixn = xNum;	
    /* number of imaged positions */
    if (!sf_getfloat ("izo", &izo))      izo = zStart;
    /* first imaged depth (in meters) */
    if (!sf_getfloat ("ixo", &ixo))      ixo = xStart;
    /* first imaged position (in meters) */
    if (!sf_getfloat ("izd", &izd))      izd = zStep;
    /* step in depth (in meters) */
    if (!sf_getfloat ("ixd", &ixd))      ixd = xStep;
    /* step in positions (in meters) */

    if (!sf_getint ("sn", &sNum)) sNum = 1;
    /* number of scattering-angles */
    if (!sf_getfloat ("so", &sStart)) sStart = 0.f;
    /* first scattering-angle */
    if (!sf_getfloat ("sd", &sStep)) sStep = 1.f;
    /* step in scattering-angles */

    if ( !sf_getbool ("isAA", &isAA) ) isAA = true;
    /* if y, apply anti-aliasing */
    if (!sf_getfloat ("dx", &dx)) dx = xStep;
    /* x-range for point detection */
    if (!sf_getfloat ("dt", &dt)) dt = 0.008f;
    /* time-range for point detection */
    if (!sf_getfloat ("xlim", &xlim)) xlim = 2 * xStep;
    /* maximum distance between depth-line points */
    if (!sf_getfloat ("xapert", &xapert)) xapert = xNum * xStep;
    /* migration aperture size */
    if (!sf_getint ("pj", &pj)) pj = 1;
    /* jump in points */


    // OUTPUT PARAMETERS
    sf_putint (resFile, "n1", izn); 
    sf_putint (resFile, "n2", ixn); 
    sf_putint (resFile, "n3", ppn); 
    sf_putint (resFile, "n4", 1); 

    sf_putfloat (resFile, "d1", izd); 
    sf_putfloat (resFile, "d2", ixd); 
    sf_putfloat (resFile, "d3", ppd); 
    sf_putfloat (resFile, "d4", 1); 

    sf_putfloat (resFile, "o1", izo); 
    sf_putfloat (resFile, "o2", ixo); 
    sf_putfloat (resFile, "o3", ppo); 
    sf_putfloat (resFile, "o4", 1); 

    // MEMORY ALLOCATION

    const int piSize    = zNum * xNum;
    const int volSize   = zNum * rNum * xNum;
    const int piResSize = izn * ixn;

    float* tVol = sf_floatalloc (volSize);
    float* xVol = sf_floatalloc (volSize);

    float* piData  = sf_floatalloc (piSize);
    float* piImage = sf_floatalloc (piResSize);

    // READ ESCAPE VOLUMES

    sf_seek (xEscFile, 0, SEEK_SET);		
    sf_floatread (xVol, volSize, xEscFile);
    sf_fileclose (xEscFile);

    sf_seek (tEscFile, 0, SEEK_SET);		
    sf_floatread (tVol, volSize, tEscFile);
    sf_fileclose (tEscFile);

    // MAIN LOOP

    DepthBackfitMigrator2D dbfmig;
    dbfmig.init (zNum, zStart, zStep, 
		 pNum, pStart, pStep,
		 xNum, xStart, xStep,
		 rNum, rStart, rStep,
		 sNum, sStart, sStep,
		 izn, izo, izd,
		 ixn, ixo, ixd,
		 dx, dt, xlim, xapert, pj,
		 xVol, tVol, isAA);

    for (int ip = 0; ip < ppn; ++ip) {
	clock_t begin=clock();
	const float curP = ppo + ip * ppd;			

	memset ( piData,  0, piSize * sizeof (int) );
	memset ( piImage, 0, piResSize * sizeof (int) );

	// read partial image
	const int pind = (curP - pStart) / pStep;
	size_t startPos = pind * piSize * sizeof(float);
	sf_seek (piFile, startPos, SEEK_SET);
	sf_floatread (piData, piSize, piFile);

	// get image
	dbfmig.processPartialImage (piData, curP, piImage);

	// write result
	startPos = ip * piResSize * sizeof(float);
	sf_seek (resFile, startPos, SEEK_SET);
	sf_floatwrite (piImage, piResSize, resFile);

	clock_t end=clock();
	float time = (end - begin) / (11 * CLOCKS_PER_SEC);
	sf_warning ("pimage %d of %d (%g s);", ip + 1, ppn, time);		}
    sf_warning(".");

    // FINISH

    free (xVol);
    free (tVol);
	
    free (piData);
    free (piImage);


    exit(0);
}
