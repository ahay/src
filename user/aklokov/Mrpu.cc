/* Reflection Parametric Unwrapping */
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

#include <rsf.hh>
#include "support.hh"
#include "signalUnwrapper.hh"

VolumeParams     dp; // data params

// files
sf_file dataFile;
sf_file outFile;

int main (int argc, char* argv[]) {

// Initialize RSF 
    sf_init (argc,argv);

// INPUT FILES
    dataFile = sf_input ("in");
    /* common-offset sections */
    outFile  = sf_output("out");
    /*  */

    // data params

    if ( !sf_histint   (dataFile, "n1", &dp.zNum)   ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile, "d1", &dp.zStep)  ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile, "o1", &dp.zStart) ) sf_error ("Need o1= in input");

    if ( !sf_histint   (dataFile, "n2", &dp.hNum)   ) sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile, "d2", &dp.hStep)  ) sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile, "o2", &dp.hStart) ) sf_error ("Need o2= in input");

    if ( !sf_histint   (dataFile, "n3", &dp.xNum)   ) sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile, "d3", &dp.xStep)  ) sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile, "o3", &dp.xStart) ) sf_error ("Need o3= in input");

//    if ( !sf_histint   (dataFile, "n4", &dp.yNum)   ) sf_error ("Need n4= in input");
//    if ( !sf_histfloat (dataFile, "d4", &dp.yStep)  ) sf_error ("Need d4= in input");
//    if ( !sf_histfloat (dataFile, "o4", &dp.yStart) ) sf_error ("Need o4= in input");

    //
    char* corUnit;
    char* unit;
	
    // data
    // time - in ms
    corUnit = (char*) "ms"; unit = sf_histstring (dataFile, "unit1"); if (!unit) sf_error ("unit1 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.zStep *= 1000; dp.zStart *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit2"); if (!unit) sf_error ("unit2 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.xStep *= 1000; dp.xStart *= 1000; }
    // crossline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit3"); if (!unit) sf_error ("unit3 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.yStep *= 1000; dp.yStart *= 1000; }
    // offset - in m
 //   corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit4"); if (!unit) sf_error ("unit4 in data file is not defined");
//    if ( strcmp (corUnit, unit) ) { dp.hStep *= 1000; dp.hStart *= 1000; }

	const int pNum = 201;
	float pStart = 0.f;
	float pStep = 10.0;

	const int tNum = dp.zNum;
	const int hNum = dp.hNum;
	const int xNum = dp.xNum;

	const int zoSize = pNum * tNum;
	const int dataSize = xNum * hNum * tNum;

    sf_putint (outFile, "n1", tNum);
	sf_putint (outFile, "n2", pNum);
    sf_putfloat (outFile, "d1", dp.zStep); 
	sf_putfloat (outFile, "d2", pStep);
    sf_putfloat (outFile, "o1", dp.zStart); 
	sf_putfloat (outFile, "o2", pStart);
    sf_putstring(outFile, "label1", "time"); sf_putstring(outFile, "label2", "inline");
    sf_putstring(outFile, "unit1", "ms"); sf_putstring(outFile, "unit2", "m");

		
    float* zo   = sf_floatalloc (zoSize);
	float* data = sf_floatalloc (dataSize);

    sf_floatread (data, dataSize, dataFile);

	SignalUnwrapper signalUnwrapper;
	signalUnwrapper.setParams (&dp, data);

	signalUnwrapper.unwrap (zo);

	sf_floatwrite (zo, zoSize, outFile);

    free (data);
    free (zo);

    return 0;
}
