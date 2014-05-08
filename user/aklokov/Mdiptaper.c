/* Aperture optimization for migrated gathers in the dip-angle domain.

   Estimates a constructive imaging part of a reflection event in the dip-angle domain.
   Basing on the estimation defines a stacking weight for every migrated sample.

   Input:
   dipFile.rsf - dips esitimated in constant-dip subimages. The dips are in degree (!).
   A positive dip corresponds to an ascending boundary, a negative dip - to a descending boundary.
   A constant-dip subimage consists of migrated traces correspondig to the same dip-angle.

   Output:
   taperFile.rsf - optimal weights for the migrated samples
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

int main (int argc, char* argv[]) 
{
/* dip-angle gathers dimensions */
    int zNum_;   float zStart_;   float zStep_;
    int dipNum_; float dipStart_; float dipStep_;
    int xNum_;	 float xStart_;   float xStep_;   
    sf_file dipFile, taperFile;
    float dz = 0.f;
    float greyarea = 0.f;
    int panelSize, id, startPos, iz, ix, ind;
    float *dipPanel, *taperPanel;
    float CONVRATIO, curDip, curDipRad, curDipSin, curDipCos;
    float curZ, ksi, ksicos2, a, b, c, sqrtD, sin1, sin2, dip1, dip2, dipLeft, dipRight;
    float secondLeft, secondRight, temp, dip, taper;

/* Initialize RSF  */
    sf_init (argc,argv);
/* Input files */
    dipFile = sf_input ("in");
/* check that the input is float  */
    if ( SF_FLOAT != sf_gettype (dipFile) ) sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

/* Output file */
    taperFile = sf_output("out");

/* Depth/time axis  */
    if ( !sf_histint   (dipFile, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dipFile, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dipFile, "o1", &zStart_) ) sf_error ("Need o1= in input");
/* Dip angle axis  */
    if ( !sf_histint   (dipFile, "n2", &xNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dipFile, "d2", &xStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dipFile, "o2", &xStart_) ) sf_error ("Need o2= in input");
/* x axis  */
    if ( !sf_histint   (dipFile, "n3", &dipNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (dipFile, "d3", &dipStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (dipFile, "o3", &dipStart_) )   sf_error ("Need o3= in input");

/* tapering paremeters */
    if ( !sf_getfloat ("dz", &dz) ) dz = 20.f;
    /* half of a migrated wave length */
    if (dz < 0) {sf_warning ("dz value is changed to 20"); dz = 20.f;}

    if ( !sf_getfloat ("greyarea", &greyarea) ) greyarea = 10.f;
    /* width of event tail taper (in degree) */
    if ( greyarea < 0 ) {sf_warning ("greyarea value is changed to 10"); greyarea = 10.f;}

    panelSize = xNum_ * zNum_;
    dipPanel   = sf_floatalloc (panelSize);
    taperPanel = sf_floatalloc (panelSize);

    CONVRATIO = SF_PI / 180.f;

    for (id = 0; id < dipNum_; ++id) {

	sf_warning ("dip %d of %d;", id + 1, dipNum_);	

	curDip = dipStart_ + id * dipStep_;
	curDipRad = curDip * CONVRATIO; 
	curDipSin = sin (curDipRad);
	curDipCos = cos (curDipRad);		

	memset (dipPanel,   0, panelSize * sizeof (float));
	memset (taperPanel, 0, panelSize * sizeof (float));

	startPos = id * panelSize * sizeof (float);
	sf_seek (dipFile, startPos, SEEK_SET);		
	sf_floatread (dipPanel, panelSize, dipFile);

	/* taper values are constant for every depth level */
	for (iz = 0; iz < zNum_; ++iz) {

	    curZ = zStart_ + iz * zStep_;
	    if (! (curZ - dz > 0) ) continue; /* out from data */

	    ksi = curZ / (curZ - dz); 			
	    ksicos2 = pow (ksi * curDipCos, 2);

	    a = pow (curDipSin, 2) + ksicos2;
	    b = -2 * curDipSin;
	    c = 1 - ksicos2;

	    sqrtD = sqrt (b*b - 4*a*c);
	
	    sin1 = (-b + sqrtD) / (2 * a); 
	    sin2 = (-b - sqrtD) / (2 * a); 
	
	    dip1 = asin (sin1) / CONVRATIO; 
	    dip2 = asin (sin2) / CONVRATIO; 
	
	    dipLeft  = dip1 < dip2 ? dip1 : dip2;
	    dipRight = dip1 < dip2 ? dip2 : dip1;

	    /* define tails tapering intervals */
	    secondLeft  = dipLeft  - greyarea;
	    secondRight = dipRight + greyarea;

	    temp = 0.f;

	    for (ix = 0; ix < xNum_; ++ix) {
		ind = ix * zNum_ + iz;
		dip = dipPanel [ind]; /* local slope in degree */
		taper = 1;
		if (dip > dipLeft && dip < dipRight) /* the slope is in a constructive zone */
		    taper = 1;
		else {	
		    if (dip < secondLeft || dip > secondRight) { /* the slope is out from the tapering zone */
			taper = 0;
		    } else {
			if (dip < dipLeft) temp = -SF_PI + SF_PI * (dip - secondLeft) * 1.f / (dipLeft - secondLeft); 
			else 		       temp = SF_PI * (dip - dipRight) * 1.f / (secondRight - dipRight); 
			taper = 0.5 + 0.5 * cos (temp);
		    }
		}
		taperPanel[ind] = taper;
	    }
	}	
	sf_floatwrite (taperPanel, panelSize, taperFile);
    }	

    sf_warning (".");

    free (dipPanel);
    free (taperPanel);

    exit(0);
}
