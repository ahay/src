/* Reflection event apex protector/removal for dip-angle gathers.

   May be used for migration aperture optimization or for reflected energy
   supression. For the last multiply the output on -1.

   Input:
   dagFile.rsf - input dip-angle gathers;
   dipFile.rsf - dips esitimated in the image domain. The dips are in degree (!)

   Output:
   taperFile.rsf - mask for input dip-angle gathers
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
    sf_file dagFile, dipFile=NULL, taperFile;
    bool isDepthDep = true;
    float pwidth = 10.f;
    float greyarea = 0.f;
    float dz = 0.f;
    int panelSize, taperSize, ix, iz, ind, id, indt;
    float CONVRATIO, temp, dip, curZ, ksi, curDip, dipLeft, dipRight;
    float curDipRad, curDipSin, curDipCos, ksicos2, a, b, c;
    float sqrtD, sin1, sin2, dip1, dip2, secondLeft, secondRight, taper;

    float* dipPanel;
    float* taperPanel;

/* Initialize RSF */
    sf_init (argc,argv);
/* Input files */
    dagFile = sf_input ("in");
/* check that the input is float  */
    if ( SF_FLOAT != sf_gettype (dagFile) ) sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

    if ( NULL != sf_getstring("dips") ) {
	/* dips esitimated in the image domain (in degree) */ 
	dipFile  = sf_input ("dips");
    }

/* Output file */
    taperFile = sf_output("out");


/* Depth/time axis  */
    if ( !sf_histint   (dagFile, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dagFile, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dagFile, "o1", &zStart_) ) sf_error ("Need o1= in input");
/* Dip angle axis  */
    if ( !sf_histint   (dagFile, "n2", &dipNum_) )     sf_error ("Need n2= in input");
    if ( !sf_histfloat (dagFile, "d2", &dipStep_) )    sf_error ("Need d2= in input");
    if ( !sf_histfloat (dagFile, "o2", &dipStart_) )   sf_error ("Need o2= in input");
/* x axis  */
    if ( !sf_histint   (dagFile, "n3", &xNum_) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (dagFile, "d3", &xStep_) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (dagFile, "o3", &xStart_) ) sf_error ("Need o3= in input");

/* tapering paremeters */
    if (!sf_getbool ("ddep", &isDepthDep)) isDepthDep = true;
    /* if y, taper depends on depth; if n, no */

    if ( !sf_getfloat ("pwidth", &pwidth) ) pwidth = 10.f;
    /* protected width (in degree) */

    if ( !sf_getfloat ("greyarea", &greyarea) ) greyarea = 10.f;
    /* width of event tail taper (in degree) */
    if ( greyarea < 0 ) {sf_warning ("greyarea value is changed to 10"); greyarea = 10.f;}
    if ( !sf_getfloat ("dz", &dz) ) dz = 20.f;
    /* half of a migrated wave length */
    if (dz < 0) {sf_warning ("dz value is changed to 20"); dz = 20.f;}

    /* input dips */
    panelSize = xNum_ * zNum_;
    dipPanel   = sf_floatalloc (panelSize);
    sf_seek (dipFile, 0, SEEK_SET);		
    sf_floatread (dipPanel, panelSize, dipFile);
    /* output taper */
    taperSize = dipNum_ * zNum_;
    taperPanel = sf_floatalloc (taperSize);

    CONVRATIO = SF_PI / 180.f;
    temp = 0.f;

    for (ix = 0; ix < xNum_; ++ix) {
	memset (taperPanel, 0, taperSize * sizeof (float));
	for (iz = 0; iz < zNum_; ++iz) {
	    ind = ix * zNum_ + iz;
	    dip = dipPanel [ind]; /* local slope in degree */


	    curZ = zStart_ + iz * zStep_;
	    if (! (curZ - dz > 0) ) continue; /* out from data */

	    ksi = curZ / (curZ - dz); 			

	    for (id = 0; id < dipNum_; ++id) {

		curDip = dipStart_ + id * dipStep_;
		indt = id * zNum_ + iz;
				
		dipLeft  = 0.f;
		dipRight = 0.f;

		if (isDepthDep) { /* depth-dependent taper - based on the equation from Landa et al., 2008 */

		    curDipRad = dip * CONVRATIO; 
		    curDipSin = sin (curDipRad);
		    curDipCos = cos (curDipRad);		
	
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
		} else { /* taper widht is constant along depths */
		    dipLeft  = dip - pwidth;
		    dipRight = dip + pwidth;
		}

		/* define tails tapering intervals */
		secondLeft  = dipLeft  - greyarea;
		secondRight = dipRight + greyarea;

		taper = 1;
		if (curDip > dipLeft && curDip < dipRight) /* the slope is in a constructive zone */
		    taper = 1;
		else {	
		    if (curDip < secondLeft || curDip > secondRight) { /* the slope is out from the tapering zone */
			taper = 0;
		    } else {
			if (curDip < dipLeft) temp = -M_PI + M_PI * (curDip - secondLeft) * 1.f / (dipLeft - secondLeft); 
			else 		       temp = M_PI * (curDip - dipRight) * 1.f / (secondRight - dipRight); 
			taper = 0.5 + 0.5 * cos (temp);
		    }
		}
		taperPanel[indt] = taper;				
	    }
	}
	sf_floatwrite (taperPanel, taperSize, taperFile);
    }	

    sf_warning (".");

    free (dipPanel);
    free (taperPanel);

    exit(0);
}
