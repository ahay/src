/* Edge tapering for dip-angle gathers
 */
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

int main (int argc, char* argv[]) {

    int zNum, dipNum, sdipNum, xNum, yNum;
    float dipStep, dipStart, sdipStep, sdipStart;
    int dataSize, taperSize, idy, idx, iy, ix, it, iz;
    float *data_, *pData, *taper_, *pTaper;
    float length, dipMax, edgeDip, curSDip, sdip2, curDip, dip2;
    float dip, w;
    sf_file dataFile, outFile;

// Initialize RSF 
    sf_init (argc, argv);

// INPUT FILES
    dataFile = sf_input ("in");
    /* common-offset sections */
    outFile = sf_output ("out");
    /*  */

    // data params
    if ( !sf_histint   (dataFile, "n1", &zNum   ) )  sf_error ("Need n1= in input");

    if ( !sf_histint   (dataFile, "n2", &dipNum   ) )  sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile, "d2", &dipStep  ) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile, "o2", &dipStart ) )  sf_error ("Need o2= in input");

    if ( !sf_histint   (dataFile, "n3", &sdipNum   ) )  sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile, "d3", &sdipStep  ) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile, "o3", &sdipStart ) )  sf_error ("Need o3= in input");

    if ( !sf_histint   (dataFile, "n4", &xNum     ) )  sf_error ("Need n4= in input");
    if ( !sf_histint   (dataFile, "n5", &yNum     ) )  sf_error ("Need n5= in input");
    
    if ( !sf_getfloat ("len", &length) ) length = 5.f;
    /* length of the taper function (in degree) */

    // build taper function
    dataSize = sdipNum * dipNum * zNum;
    data_  = sf_floatalloc (dataSize);

    taperSize = sdipNum * dipNum;
    taper_ = sf_floatalloc (taperSize);
    pTaper = taper_;

    dipMax = dipStep * dipNum / 2;
    edgeDip = dipMax - length;

    for (idy = 0; idy < sdipNum; ++idy) {
        curSDip = sdipStart + idy * sdipStep;
	sdip2 = curSDip * curSDip;
	for (idx = 0; idx < dipNum; ++idx, ++pTaper) {
    	    curDip = dipStart + idx * dipStep;
	    dip2 = curDip * curDip;

	    // stacking taper	
	    dip = sqrt (sdip2 + dip2);
	    w = 1.f;		
	    if (dip > edgeDip) {
		if (dip > dipMax) w = 0.f;
		else w = 1.f - (dip - edgeDip) / length;
	    }	
	    *pTaper = w;
	}
    }

    // tapering
    for (iy = 0; iy < yNum; ++iy) {
	for (ix = 0; ix < xNum; ++ix) {
	    // read a gather
	    sf_floatread (data_, dataSize, dataFile);
	    // apply the taper			
	    pTaper = taper_;
	    pData  = data_;	
	    for (it = 0; it < taperSize; ++it, ++pTaper) {
		w = *pTaper;
		for (iz = 0; iz < zNum; ++iz, ++pData)
		    *pData *= w;
	    }
	    // write the gather
	    sf_floatwrite (data_, dataSize, outFile);
	}
    }

    free (data_);
    free (taper_);

    sf_fileclose (dataFile);
    sf_fileclose (outFile);

    return 0;
}
