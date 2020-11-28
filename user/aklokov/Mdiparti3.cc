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
int xNum_;	   float xStart_;   float xStep_;
int yNum_;	   float yStart_;   float yStep_;
int dipNum_;   float dipStart_; float dipStep_;
int sdipNum_;  float sdipStart_; float sdipStep_;
// secondary images
int xppn;      float xppo;      float xppd;
int yppn;      float yppo;      float yppd;
int itn;       float ito;       float itd;
int ixn;       float ixo;       float ixd;
int iyn;       float iyo;       float iyd;
// velocity model dimensions:
int   v_tNum_; float v_tStart_; float v_tStep_;
int   v_xNum_; float v_xStart_; float v_xStep_;
int   v_yNum_; float v_yStart_; float v_yStep_;
// run parameters
int sembWindow_;
float apert_;
float gamma_;

void processPartImage (const float migDipY, const float migDipX, float* partImage, 
					   float* dPartImage, float* sembMap, float* velModel) {

	const float curDipRadY = migDipY * SF_PI / 180.f;
	const float curDipRadX = migDipX * SF_PI / 180.f;

	// loop over y
	for (int iy = 0; iy < iyn; ++iy) {
		const float curY = iyo + iy * iyd;
		const int velIndY = (curY - v_yStart_) / v_yStep_ * v_xNum_	* v_tNum_;
			
		// loop over x
		for (int ix = 0; ix < ixn; ++ix) {
			const float curX = ixo + ix * ixd;
			const int velIndX = (curX - v_xStart_) / v_xStep_ * v_tNum_;		
			// loop over z
			float* trace2 = sf_floatalloc (tNum_);
			memset (trace2, 0, tNum_ * sizeof (float) );

			const size_t shiftForTrace = (iy * ixn + ix) * itn;

			int tcountMain = 0;

#ifdef _OPENMP 
#pragma omp parallel for
#endif		
			for (int iz = 0; iz < itn; ++iz) {
				const float curT = ito + iz * itd;
				// get velocity
				const int velIndT = (curT - v_tStart_) / v_tStep_;
				float velMig = *(velModel + velIndY + velIndX + velIndT);

				const float zd = 0.5 * velMig * curT / gamma_ + 1e-6;		

				float diffStack  = 0.f;
				float diffStack2 = 0.f;
				int tcount = 0;				
				// loop over y - the second one 
				for (int ipy = 0; ipy < yNum_; ++ipy) {
					const float curPosY = yStart_ + ipy * yStep_;	
					const float dy = curPosY - curY;
					if (fabs (dy) > apert_) continue;
					// loop over x - the second one 
					for (int ipx = 0; ipx < xNum_; ++ipx) {
						const float curPosX = xStart_ + ipx * xStep_;	
						const float dx = curPosX - curX;
						if (fabs (dx) > apert_) continue;

						// inside aperture

						const float xiX = dx / zd;				
						const float xiY = dy / zd;				
						const float k2 = xiY * tan (curDipRadY) + xiX * tan (curDipRadX);
						const float ke = xiY*xiY + xiX*xiX + 1;

						const float t = curT * (k2 + sqrt (k2 * k2 + ke));

						const int tInd = (t - tStart_) / tStep_;	
						if (tInd < 0 || tInd >= tNum_)
							continue;

						const int ind = (ipy * xNum_ + ipx) * tNum_ + tInd;
						const float sample = partImage [ind];		

						diffStack  += sample;
						diffStack2 += sample * sample;
						++tcount;
					}
				}
					
				if (tcount > tcountMain) tcountMain = tcount;

				const int sInd = shiftForTrace + iz;
				dPartImage [sInd] += diffStack;
				trace2 [iz] += diffStack2;
			}

			float* sembTrace = sf_floatalloc (itn);
			Sembler::getSemblanceForTrace (tcountMain, dPartImage + shiftForTrace, trace2, itn, sembWindow_, sembTrace);		
//			Sembler::getSemblanceForTrace (tcount, dPartImage + shiftForTrace, trace2, itn, sembWindow_, sembMap + shiftForTrace);		
			float* pMap = sembMap + shiftForTrace;
			float* pTrace = sembTrace;
			for (int iz = 0; iz < itn; ++iz, ++pMap, ++pTrace)
				*pMap = *pTrace;

			free (sembTrace);
			free (trace2);
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
// x-axis
    if ( !sf_histint   (piFile, "n2", &xNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (piFile, "d2", &xStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (piFile, "o2", &xStart_) ) sf_error ("Need o2= in input");
// y-axis
    if ( !sf_histint   (piFile, "n3", &yNum_) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (piFile, "d3", &yStep_) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (piFile, "o3", &yStart_) ) sf_error ("Need o3= in input");
// x-dip
    if ( !sf_histint   (piFile, "n4", &dipNum_) )     sf_error ("Need n4= in input");
    if ( !sf_histfloat (piFile, "d4", &dipStep_) )    sf_error ("Need d4= in input");
    if ( !sf_histfloat (piFile, "o4", &dipStart_) )   sf_error ("Need o4= in input");
// y-dip
    if ( !sf_histint   (piFile, "n5", &sdipNum_) )     sf_error ("Need n5= in input");
    if ( !sf_histfloat (piFile, "d5", &sdipStep_) )    sf_error ("Need d5= in input");
    if ( !sf_histfloat (piFile, "o5", &sdipStart_) )   sf_error ("Need o5= in input");

    // VELOCITY MODEL PARAMS
    if ( !sf_histint   (velFile, "n1", &v_tNum_) )   sf_error ("Need n1= in velocity file");
    if ( !sf_histfloat (velFile, "d1", &v_tStep_) )  sf_error ("Need d1= in velocity file");
    if ( !sf_histfloat (velFile, "o1", &v_tStart_) ) sf_error ("Need o1= in velocity file");

    if ( !sf_histint   (velFile, "n2", &v_xNum_) )   sf_error ("Need n2= in velocity file");
    if ( !sf_histfloat (velFile, "d2", &v_xStep_) )  sf_error ("Need d2= in velocity file");
    if ( !sf_histfloat (velFile, "o2", &v_xStart_) ) sf_error ("Need o2= in velocity file");

    if ( !sf_histint   (velFile, "n3", &v_yNum_) )   sf_error ("Need n3= in velocity file");
    if ( !sf_histfloat (velFile, "d3", &v_yStep_) )  sf_error ("Need d3= in velocity file");
    if ( !sf_histfloat (velFile, "o3", &v_yStart_) ) sf_error ("Need o3= in velocity file");


    char* corUnit;
    char* unit;
    corUnit = (char*) "ms"; unit = sf_histstring (velFile, "unit1"); if (!unit) sf_error ("unit1 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_tStep_ *= 1000; v_tStart_ *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit2"); if (!unit) sf_error ("unit2 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_xStep_ *= 1000; v_xStart_ *= 1000; }

	// DIP PARAMS
    if (!sf_getint ("xppn", &xppn)) xppn = dipNum_;
	/* number of processed partial images */
    if (!sf_getfloat ("xppo", &xppo)) xppo = dipStart_;
	/* first processed partial image */
    if (!sf_getfloat ("xppd", &xppd)) xppd = dipStep_;
	/* step in processed partial images */

    if (!sf_getint ("yppn", &yppn)) yppn = sdipNum_;
	/* number of processed partial images */
    if (!sf_getfloat ("yppo", &yppo)) yppo = sdipStart_;
	/* first processed partial image */
    if (!sf_getfloat ("yppd", &yppd)) yppd = sdipStep_;
	/* step in processed partial images */

    // IMAGE PARAMS
    if (!sf_getint ("itn", &itn))        itn = tNum_;	
    /* number of imaged depth samples */
    if (!sf_getint ("ixn", &ixn))        ixn = xNum_;	
    /* number of imaged positions */
    if (!sf_getint ("iyn", &iyn))        iyn = yNum_;	
    /* number of imaged positions */
    if (!sf_getfloat ("ito", &ito))      ito = tStart_;
    /* first imaged time (in ms) */
    if (!sf_getfloat ("ixo", &ixo))      ixo = xStart_;
    /* first imaged position (in m) */
    if (!sf_getfloat ("iyo", &iyo))      iyo = yStart_;
    /* first imaged position (in m) */
    if (!sf_getfloat ("itd", &itd))      itd = tStep_;
    /* step in time (in ms) */
    if (!sf_getfloat ("ixd", &ixd))      ixd = xStep_;
    /* step in positions (in m) */
    if (!sf_getfloat ("iyd", &iyd))      iyd = yStep_;
    /* step in positions (in m) */

	// OUTPUT PARAMETERS
  	sf_putint (resFile, "n1", itn); 
  	sf_putint (resFile, "n2", ixn); 
  	sf_putint (resFile, "n3", iyn); 
  	sf_putint (resFile, "n4", xppn); 
  	sf_putint (resFile, "n5", yppn); 
  	
  	sf_putfloat (resFile, "d1", itd); 
  	sf_putfloat (resFile, "d2", ixd); 
  	sf_putfloat (resFile, "d3", iyd); 
  	sf_putfloat (resFile, "d4", xppd); 
	sf_putfloat (resFile, "d5", yppd); 

	sf_putfloat (resFile, "o1", ito); 
  	sf_putfloat (resFile, "o2", ixo); 
  	sf_putfloat (resFile, "o3", iyo); 
  	sf_putfloat (resFile, "o4", xppo); 
  	sf_putfloat (resFile, "o5", yppo); 

  	sf_putint (sembFile, "n1", itn); 
  	sf_putint (sembFile, "n2", ixn); 
  	sf_putint (sembFile, "n3", iyn); 
  	sf_putint (sembFile, "n4", xppn); 
  	sf_putint (sembFile, "n5", yppn); 
  	
  	sf_putfloat (sembFile, "d1", itd); 
  	sf_putfloat (sembFile, "d2", ixd); 
  	sf_putfloat (sembFile, "d3", iyd); 
  	sf_putfloat (sembFile, "d4", xppd); 
	sf_putfloat (sembFile, "d5", yppd); 

	sf_putfloat (sembFile, "o1", ito); 
  	sf_putfloat (sembFile, "o2", ixo); 
  	sf_putfloat (sembFile, "o3", iyo); 
  	sf_putfloat (sembFile, "o4", xppo); 
  	sf_putfloat (sembFile, "o5", yppo); 



// main part

	const size_t inSize = yNum_ * xNum_ * tNum_;
	float* partImage  = sf_floatalloc (inSize);
	const int outSize = iyn * ixn * itn;
	float* dPartImage = sf_floatalloc (outSize);
	float* sembMap  = sf_floatalloc (outSize);

	const int velSize = v_tNum_ * v_xNum_ * v_yNum_;
	float* vel = sf_floatalloc (velSize);
	sf_seek (velFile, 0, SEEK_SET);
	sf_floatread (vel, velSize, velFile);
	if (!isVelKMS) {
		float* ptrVel = vel;
		for (int iv = 0; iv < velSize; ++iv, ++ptrVel)
			*ptrVel *= 0.001f;
	}
	
	const int migNum = yppn * xppn;
	for (int idy = 0; idy < yppn; ++idy) {
		const float curDipY = yppo + idy * yppd;
		const int idpy = (curDipY - sdipStart_) / sdipStep_;
		for (int idx = 0; idx < xppn; ++idx) {
			const float curDipX = xppo + idx * xppd;
			const int idpx = (curDipX - dipStart_) / dipStep_;

			sf_warning ("dip %d of %d;", idy * xppn + idx + 1, migNum);	

			memset ( partImage,  0, inSize * sizeof (float) );
			memset ( dPartImage, 0, outSize * sizeof (float) );
			memset ( sembMap,    0, outSize * sizeof (float) );

			// read partial image
			const size_t startPos = (idpy * dipNum_ + idpx) * inSize * sizeof (float);
			sf_seek (piFile, startPos, SEEK_SET);		
			sf_floatread (partImage, inSize, piFile);

			// process partial image
			processPartImage (curDipY, curDipX, partImage, dPartImage, sembMap, vel);

			// write down the result
			sf_floatwrite (dPartImage, outSize, resFile);
			sf_floatwrite (sembMap, outSize, sembFile);
		}	
	}
	sf_warning (".");

	sf_fileclose (piFile);
	sf_fileclose (resFile);

	free (partImage);
	free (dPartImage);
	free (sembMap);

	return 0;
}
