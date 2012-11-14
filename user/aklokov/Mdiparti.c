/* diparti
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

// dip-angle gathers dimensions
int zNum_;   float zStart_;   float zStep_;
int dipNum_; float dipStart_; float dipStep_;
int xNum_;	 float xStart_;   float xStep_;

// velocity model dimensions:
int   v_tNum_;                 
float v_tStart_;
float v_tStep_;

int   v_xNum_;
float v_xStart_;
float v_xStep_;

float apert_;

void processPartImage (const float migDip, float* partImage, float* dPartImage, float* velModel) {

	const float curDipRad = migDip * SF_PI / 180.f;

	// loop over x
	for (int ix = 0; ix < xNum_; ++ix) {
		const float curX = xStart_ + ix * xStep_;
		// loop over z
		for (int iz = 0; iz < zNum_; ++iz) {
			const float curT = zStart_ + iz * zStep_;
			// get velocity
			const int velInd = (curT - v_tStart_) / v_tStep_;
			float velMig = *(velModel + velInd);

			const float zd = 0.5 * velMig * curT + 1e-6;		

			float diffStack = 0.f;
			
			// loop over x - the second one 
			for (int ip = 0; ip < xNum_; ++ip) {
				const float curPos = xStart_ + ip * xStep_;	
				const float dx = curX - curPos;

				if (fabs (dx) > apert_) continue;

				const float xi = dx / zd;				

				const float t = curT * (xi * tan (curDipRad) + sqrt ( pow (xi / cos (curDipRad), 2) + 1) );

				const int tInd = (t - zStart_) / zStep_;	
				if (tInd < 0 || tInd >= zNum_)
						continue;

				int ind = ip * zNum_ + tInd;
		
				diffStack += partImage [ind];							
			}

			int sInd = ix * zNum_ + iz;
			dPartImage [sInd] += diffStack;
		}
	}

	return;
}

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc,argv);
// Input files
    sf_file piFile = sf_input ("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (piFile) ) sf_error ("Need float input: constant-dip partial images");
    /* constant-dip partial images */

	sf_file velFile;
	bool  isVelKMS;
    if ( NULL != sf_getstring("vel") ) {
	/* velocity model file (velocity in km/s) */ 
		velFile = sf_input ("vel");
		float firstvel;
		sf_floatread (&firstvel, 1, velFile);
		isVelKMS = true;		
		if (firstvel > 15.f) {
		    sf_warning ("it seems that velocity is in m/s - will be divided by 1000");	
		    isVelKMS = false;					
		}			
    } else { sf_error ("Need input: velocity model"); }

// Output file
    sf_file resFile = sf_output ("out");

    if ( !sf_getfloat ("apert", &apert_) ) apert_ = 1000;
    /* diffraction summation aperture */

// Depth/time axis 
    if ( !sf_histint   (piFile, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (piFile, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (piFile, "o1", &zStart_) ) sf_error ("Need o1= in input");
// Dip angle axis 
    if ( !sf_histint   (piFile, "n2", &xNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (piFile, "d2", &xStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (piFile, "o2", &xStart_) ) sf_error ("Need o2= in input");
// x axis 
    if ( !sf_histint   (piFile, "n3", &dipNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (piFile, "d3", &dipStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (piFile, "o3", &dipStart_) )   sf_error ("Need o3= in input");

    // VELOCITY MODEL PARAMS
    if ( !sf_histint   (velFile, "n1", &v_tNum_) )   sf_error ("Need n1= in velocity file");
    if ( !sf_histfloat (velFile, "d1", &v_tStep_) )  sf_error ("Need d1= in velocity file");
    if ( !sf_histfloat (velFile, "o1", &v_tStart_) ) sf_error ("Need o1= in velocity file");

    if ( !sf_histint   (velFile, "n2", &v_xNum_) )   sf_error ("Need n2= in velocity file");
    if ( !sf_histfloat (velFile, "d2", &v_xStep_) )  sf_error ("Need d2= in velocity file");
    if ( !sf_histfloat (velFile, "o2", &v_xStart_) ) sf_error ("Need o2= in velocity file");

    char* corUnit;
    char* unit;
    corUnit = (char*) "ms"; unit = sf_histstring (velFile, "unit1"); if (!unit) sf_error ("unit1 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_tStep_ *= 1000; v_tStart_ *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit2"); if (!unit) sf_error ("unit2 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_xStep_ *= 1000; v_xStart_ *= 1000; }

// main part

	const int piSize = xNum_ * zNum_;
	float* partImage  = sf_floatalloc (piSize);
	float* dPartImage = sf_floatalloc (piSize);

	const int velSize = v_tNum_ * v_xNum_;
	float* vel = sf_floatalloc (velSize);
	sf_seek (velFile, 0, SEEK_SET);
	sf_floatread (vel, velSize, velFile);
	if (!isVelKMS) {
		float* ptrVel = vel;
		for (int iv = 0; iv < velSize; ++iv, ++ptrVel)
			*ptrVel *= 0.001f;
	}

	for (int id = 0; id < dipNum_; ++id) {

		sf_warning ("dip %d of %d;", id + 1, dipNum_);	

		const float curDip = dipStart_ + id * dipStep_;

		memset ( partImage,  0, piSize * sizeof (float) );
		memset ( dPartImage, 0, piSize * sizeof (float) );

		// read partial image
		const int startPos = id * piSize * sizeof (float);
		sf_seek (piFile, startPos, SEEK_SET);		
		sf_floatread (partImage, piSize, piFile);

		// process partial image
		processPartImage (curDip, partImage, dPartImage, vel);

		// write down the result
		sf_floatwrite (dPartImage, piSize, resFile);
	}	

	sf_warning (".");

	sf_fileclose (piFile);
	sf_fileclose (resFile);

	free (partImage);
	free (dPartImage);

	return 0;
}
