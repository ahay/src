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

#include <rsf.hh>
#include "sembler.hh"
#ifdef _OPENMP
#include <omp.h>
#endif


// dip-angle gathers dimensions
int tNum_;     float tStart_;   float tStep_;
int dipNum_;   float dipStart_; float dipStep_;
int xNum_;	   float xStart_;   float xStep_;
// secondary images
int ppn;       float ppo;       float ppd;
int itn;       float ito;       float itd;
int ixn;       float ixo;       float ixd;
// velocity model dimensions:
int   v_tNum_; float v_tStart_; float v_tStep_;
int   v_xNum_; float v_xStart_; float v_xStep_;
// run parameters
int sembWindow_;
float apert_;
float gamma_;

void processPartImage (const float migDip, float* partImage, float* dPartImage, float* sembMap, float* velModel) {

	const float curDipRad = migDip * SF_PI / 180.f;

	// loop over x
	for (int ix = 0; ix < ixn; ++ix) {
		const float curX = ixo + ix * ixd;
		// loop over z
		float* trace2 = sf_floatalloc (tNum_);
		memset (trace2, 0, tNum_ * sizeof (float) );
#ifdef _OPENMP 
#pragma omp parallel for
#endif		
		for (int iz = 0; iz < itn; ++iz) {
			const float curT = ito + iz * itd;
			// get velocity
			const int velInd = (curT - v_tStart_) / v_tStep_;
			float velMig = *(velModel + velInd);

			const float zd = 0.5 * velMig * curT / gamma_ + 1e-6;		

			float diffStack  = 0.f;
			float diffStack2 = 0.f;
			// loop over x - the second one 
			int count = 0;
			for (int ip = 0; ip < xNum_; ++ip) {
				const float curPos = xStart_ + ip * xStep_;	
				const float dx = curX - curPos;

				if (fabs (dx) > apert_) continue;

				const float xi = dx / zd;				

				const float aux = 1 - pow (gamma_ * sin (curDipRad), 2);
				const float t = curT * cos (curDipRad) * (xi * gamma_ * sin (curDipRad) + sqrt (xi*xi + aux) ) / aux;

				const int tInd = (t - tStart_) / tStep_;	
				if (tInd < 0 || tInd >= tNum_)
						continue;

				const int ind = ip * tNum_ + tInd;
				const float sample = partImage [ind];		

				diffStack  += sample;
				diffStack2 += sample * sample;
				++count;
			}

			const int sInd = ix * itn + iz;
			dPartImage [sInd] += diffStack / count;

			trace2 [iz] += diffStack2 / count;
		}

		float* sembTrace = sf_floatalloc (itn);
		Sembler::getSemblanceForTrace (xNum_, dPartImage + ix*itn, trace2, itn, sembWindow_, sembTrace);		
		float* pMap = sembMap + ix*itn;
		float* pTrace = sembTrace;
		for (int iz = 0; iz < itn; ++iz, ++pMap, ++pTrace)
			*pMap = *pTrace;
		free (sembTrace);
		free (trace2);
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

	sf_file velFile = NULL;
	bool  isVelKMS = true;
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

// Output files
    sf_file resFile = sf_output ("out");

    sf_file sembFile = NULL;
    if ( NULL != sf_getstring("semb") ) {
	/* output file containing semblance */ 
		sembFile  = sf_output ("semb"); 
    } else sf_error ("Need float output: semblance");


    if ( !sf_getfloat ("apert", &apert_) ) apert_ = 1000;
    /* diffraction summation aperture */
    if ( !sf_getint ("sembWindow", &sembWindow_) ) sembWindow_ = 11;
    /* vertical window for semblance calculation (in samples) */
    if ( !sf_getfloat ("gamma", &gamma_) ) gamma_ = 1.f;
    /* velocity-model-accuracy parameter */


// Depth/time axis 
    if ( !sf_histint   (piFile, "n1", &tNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (piFile, "d1", &tStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (piFile, "o1", &tStart_) ) sf_error ("Need o1= in input");
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

    if (!sf_getint ("ppn", &ppn)) ppn = dipNum_;
	/* number of processed partial images */
    if (!sf_getfloat ("ppo", &ppo)) ppo = dipStart_;
	/* first processed partial image */
    if (!sf_getfloat ("ppd", &ppd)) ppd = dipStep_;
	/* step in processed partial images */

    // IMAGE PARAMS
    if (!sf_getint ("itn", &itn))        itn = tNum_;	
    /* number of imaged depth samples */
    if (!sf_getint ("ixn", &ixn))        ixn = xNum_;	
    /* number of imaged positions */
    if (!sf_getfloat ("ito", &ito))      ito = tStart_;
    /* first imaged time (in ms) */
    if (!sf_getfloat ("ixo", &ixo))      ixo = xStart_;
    /* first imaged position (in m) */
    if (!sf_getfloat ("itd", &itd))      itd = tStep_;
    /* step in time (in ms) */
    if (!sf_getfloat ("ixd", &ixd))      ixd = xStep_;
    /* step in positions (in m) */

	// OUTPUT PARAMETERS
  	sf_putint (resFile, "n1", itn); 
  	sf_putint (resFile, "n2", ixn); 
  	sf_putint (resFile, "n3", ppn); 
  	sf_putint (resFile, "n4", 1); 

  	sf_putfloat (resFile, "d1", itd); 
  	sf_putfloat (resFile, "d2", ixd); 
  	sf_putfloat (resFile, "d3", ppd); 
  	sf_putfloat (resFile, "d4", 1); 

	sf_putfloat (resFile, "o1", ito); 
  	sf_putfloat (resFile, "o2", ixo); 
  	sf_putfloat (resFile, "o3", ppo); 
  	sf_putfloat (resFile, "o4", 1); 

  	sf_putint (sembFile, "n1", itn); 
  	sf_putint (sembFile, "n2", ixn); 
  	sf_putint (sembFile, "n3", ppn); 
  	sf_putint (sembFile, "n4", 1); 

  	sf_putfloat (sembFile, "d1", itd); 
  	sf_putfloat (sembFile, "d2", ixd); 
  	sf_putfloat (sembFile, "d3", ppd); 
  	sf_putfloat (sembFile, "d4", 1); 

	sf_putfloat (sembFile, "o1", ito); 
  	sf_putfloat (sembFile, "o2", ixo); 
  	sf_putfloat (sembFile, "o3", ppo); 
  	sf_putfloat (sembFile, "o4", 1); 

// main part

	const int inSize = xNum_ * tNum_;
	float* partImage  = sf_floatalloc (inSize);
	const int outSize = ixn * itn;
	float* dPartImage = sf_floatalloc (outSize);
	float* sembMap  = sf_floatalloc (outSize);

	const int velSize = v_tNum_ * v_xNum_;
	float* vel = sf_floatalloc (velSize);
	sf_seek (velFile, 0, SEEK_SET);
	sf_floatread (vel, velSize, velFile);
	if (!isVelKMS) {
		float* ptrVel = vel;
		for (int iv = 0; iv < velSize; ++iv, ++ptrVel)
			*ptrVel *= 0.001f;
	}

	for (int id = 0; id < ppn; ++id) {

		sf_warning ("dip %d of %d;", id + 1, ppn);	

		const float curDip = ppo + id * ppd;

		memset ( partImage,  0, inSize * sizeof (float) );

		memset ( dPartImage, 0, outSize * sizeof (float) );
		memset ( sembMap,    0, outSize * sizeof (float) );

		// read partial image
		const int idp = (curDip - dipStart_) / dipStep_;

		const int startPos = idp * inSize * sizeof (float);
		sf_seek (piFile, startPos, SEEK_SET);		
		sf_floatread (partImage, inSize, piFile);

		// process partial image
		processPartImage (curDip, partImage, dPartImage, sembMap, vel);

		// write down the result
		sf_floatwrite (dPartImage, outSize, resFile);
		sf_floatwrite (sembMap, outSize, sembFile);
	}	

	sf_warning (".");

	sf_fileclose (piFile);
	sf_fileclose (resFile);

	free (partImage);
	free (dPartImage);
	free (sembMap);

	return 0;
}
