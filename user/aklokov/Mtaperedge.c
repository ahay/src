/* Taper based on data parameters
Input - "inline/xline" plane
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

    int is, iy, ix, topInd, p, indBot, indLeft, indRight;
    float temp;

// Initialize RSF 
    sf_init (argc, argv);

// INPUT FILES
    sf_file dataFile = sf_input ("in");
    /* common-offset sections */
    sf_file maskFile = sf_output ("out");
    /*  */

    // data params
	int xNum, yNum;                 
    if ( !sf_histint   (dataFile, "n1", &xNum) ) sf_error ("Need n1= in input");
    if ( !sf_histint   (dataFile, "n2", &yNum) ) sf_error ("Need n2= in input");
    
	int length;
	if ( !sf_getint ("len", &length) ) length = 11;
    /* length of the taper function */

    const int surfSize = xNum * yNum;

    float* data = sf_floatalloc (surfSize);
    float* mask = sf_floatalloc (surfSize);

	float* w1 = sf_floatalloc (surfSize);
	float* w2 = sf_floatalloc (surfSize);
    for (is = 0; is < surfSize; ++is) {
		w1[is] = 1.f;
		w2[is] = 1.f;
    }

    sf_floatread (data, surfSize, dataFile);
    for (is = 0; is < surfSize; ++is) {
		mask[is] = data[is] ? 1.f : 0.f;
    }

    // TOP - BOTTOM
	for (iy = 0; iy < yNum; ++iy) {
		// top
		ix = 0;
		while (!mask[iy*xNum + ix] && ix < xNum) ++ix;
		if (ix != xNum) {
			topInd = ix;	
			for (p = 0; p < length && ix < xNum; ++p, ++ix) {
				temp = SF_PI - SF_PI * (ix - topInd) * 1.f / length;
				w1 [iy*xNum + ix] = 0.5 + 0.5 * cos (temp);
	    	}  		
		}	
		// bottom
		ix = xNum - 1;
		while (!mask[iy*xNum + ix] && ix > -1) --ix;
		if (ix > -1) {
		    indBot = ix;	
			for (p = 0; p < length && ix > -1; ++p, --ix) {
				temp = SF_PI - SF_PI * (ix - indBot) * 1.f / length;
				w1 [iy*xNum + ix] = 0.5 + 0.5 * cos (temp);
			}  		
		}
    }

    // LEFT - RIGHT
    for (ix = 0; ix < xNum; ++ix) {
		// left
		iy = 0;
		while (!mask[iy*xNum + ix] && iy < yNum) ++iy;
		if (iy != yNum) {
			indLeft = iy;	
			for (p = 0; p < length && iy < yNum; ++p, ++iy) {
				temp = SF_PI - SF_PI * (iy - indLeft) * 1.f / length;
				w2 [iy*xNum + ix] = 0.5 + 0.5 * cos (temp);
			}  		
		}
		// right
		iy = yNum - 1;
		while (!mask[iy*xNum + ix] && iy > -1) --iy;
		if (iy > -1) {
		    indRight = iy;	
		    for (p = 0; p < length && iy > -1; ++p, --iy) {
				temp = SF_PI - SF_PI * (iy - indRight) * 1.f / length;
				w2 [iy*xNum + ix] = 0.5 + 0.5 * cos (temp);
		    }  		
		}
    }

	// merge two weight functions
    for (is = 0; is < surfSize; ++is) {
		mask[is] *= (w1[is] * w2[is]);
    }

	sf_floatwrite (mask, surfSize, maskFile);

	free (data);
	free (mask);

	free (w1);
	free (w2);

	sf_fileclose (dataFile);
	sf_fileclose (maskFile);

	return 0;
}
