/* Diffraction velocity analysis

   Input:
   dataFile_.rsf - migrated dip-angle gathers

   Output:
   sembFile_.rsf - semblance spectrum
*/

/*
  Copyright (C) 2011 University of Texas at Austin
  
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

// VARIABLES

// files
sf_file dataFile_;          // dip-angle gathers file - stack in the scattering-angle direction
sf_file velFile_;           // output file - semblance
sf_file sembFile_;          // output file - semblance

// data
float* pData_;               //  --"--  --"--  square gathers data 
float* pVel_;
float* pSembPanel_;          //  --"--  --"--  result - crs-based semblance

// dip-angle gathers dimensions:
int   tNum_;                 
float tStart_;
float tStep_;

int   dipNum_;
float dipStart_;
float dipStep_;

int   xNum_;
float xStart_;
float xStep_;

// velocity model dimensions:
int   v_tNum_;                 
float v_tStart_;
float v_tStep_;

int   v_xNum_;
float v_xStart_;
float v_xStep_;

// scan parameters
int   gammaNum_;
float gammaStep_;
float gammaStart_;

int   distNum_;
float distStart_;

int   coher_;
float dlim_;
bool  isSemb_;

bool  isVelKMS;

// FUNCTIONS

void processGatherStack (float* velTrace, int cigNum, float distShift) 
{
    int ig, it, velInd, idist, id, tIndBase, ic, tInd, ind;
    float curGamma, curT, velMig, curStack, curDist, curDip, curDipRad, gammaSin, aux, zd, xi, t, val;

    float* ptrSemb = pSembPanel_;

    for (ig = 0; ig < gammaNum_; ++ig) {
	curGamma = gammaStart_ + ig * gammaStep_;
	for (it = 0; it < tNum_; ++it, ++ptrSemb) {			
	    curT = tStart_ + it * tStep_;

	    velInd = (curT - v_tStart_) / v_tStep_;
	    velMig = *(velTrace + velInd);

	    curStack = 0.f;
	    for (idist = 0; idist < cigNum; ++idist) {
		curDist = distStart_ + idist * xStep_ + distShift;
		for (id = 0; id < dipNum_; ++id) {
		    curDip = dipStart_ + id * dipStep_;
		    if (fabs (curDip) > dlim_) continue;
		    curDipRad = curDip * SF_PI / 180.f;

		    gammaSin = curGamma * sin (curDipRad);
		    aux = 1 - pow (gammaSin, 2);
		    zd = velMig * curT / (2 * curGamma) + 1e-6;
		    xi = curDist / zd;
		    
		    t = curT * cos (curDipRad) * ( xi * gammaSin + sqrt (xi*xi + aux ) ) / aux;
		    
		    tIndBase = (t - tStart_) / tStep_;
	
		    for (ic = 0; ic < coher_; ++ic) {
			tInd = tIndBase + ic - coher_ / 2;
			if (tInd < 0 || tInd >= tNum_)
			    continue;
			ind = idist * tNum_ * dipNum_ + id * tNum_ + tInd;
			val = *(pData_ + ind);

			curStack += val;
		    }
		}
	    }	
	    *ptrSemb = fabs (curStack);
	}
    }

    return;
}

void processGatherSemb (float* velTrace, int cigNum, float distShift) 
{
    int ig, it, velInd, count, idist, id, tIndBase, ic, tInd, ind;
    float curGamma,  curT, velMig, curDist, curDip, curDipRad, gammaSin, aux, zd, xi, t, val;
    float fullInput, fullOutput, curSemb;

    float* ptrSemb = pSembPanel_;

    float* energyInput  = sf_floatalloc (coher_);
    float* energyOutput = sf_floatalloc (coher_);

    for (ig = 0; ig < gammaNum_; ++ig) {
	curGamma = gammaStart_ + ig * gammaStep_;
	for (it = 0; it < tNum_; ++it, ++ptrSemb) {			
	    curT = tStart_ + it * tStep_;

	    velInd = (curT - v_tStart_) / v_tStep_;
	    velMig = *(velTrace + velInd);

	    memset ( energyInput,  0, coher_ * sizeof (float) );   
	    memset ( energyOutput, 0, coher_ * sizeof (float) );   

	    count = 0;

	    for (idist = 0; idist < cigNum; ++idist) {
		curDist = distStart_ + idist * xStep_ + distShift;
		for (id = 0; id < dipNum_; ++id) {
		    curDip = dipStart_ + id * dipStep_;
		    if (fabs (curDip) > dlim_) continue;
		    curDipRad = curDip * SF_PI / 180.f;

		    gammaSin = curGamma * sin (curDipRad);
		    aux = 1 - pow (gammaSin, 2);
		    zd = velMig * curT / (2 * curGamma) + 1e-6;
		    xi = curDist / zd;

		    t = curT * cos (curDipRad) * ( xi * gammaSin + sqrt (xi*xi + aux ) ) / aux;

		    tIndBase = (t - tStart_) / tStep_;
	
		    for (ic = 0; ic < coher_; ++ic) {
			tInd = tIndBase + ic - coher_ / 2;
			if (tInd < 0 || tInd >= tNum_)
			    continue;
			ind = idist * tNum_ * dipNum_ + id * tNum_ + tInd;
			val = *(pData_ + ind);

			energyOutput [ic] += val;
			energyInput  [ic] += val * val;
		    }
		    ++count;
		}
	    }	
	    fullInput  = 0.f;
	    fullOutput = 0.f;	
	    for (ic = 0; ic < coher_; ++ic) {
		fullInput  += energyInput [ic];
		fullOutput += pow (energyOutput [ic], 2);
	    }
	    curSemb = fullInput ? fullOutput / (count * fullInput) : 0.f;
	    *ptrSemb = curSemb;
	}
    }

    free (energyInput);
    free (energyOutput);

    return;
}

int main (int argc, char* argv[]) 
{
    float firstvel, geoX, distShift;
    // time - in ms
    char* corUnit;
    char* unit;
    int dagSize, sembSize, velSize, halfXApp, iv, ix, xNum, temp, curPanelSize, offset, velInd;
    float *ptrVel, *velTrace;
    size_t startInd;

// Initialize RSF 
    sf_init (argc, argv);
// Input files
    dataFile_   = sf_input("in");

// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

    if ( NULL != sf_getstring("vel") ) {
	/* velocity model file (velocity in km/s) */ 
	velFile_  = sf_input ("vel");
	sf_floatread (&firstvel, 1, velFile_);
	isVelKMS = true;		
	if (firstvel > 15.f) {
	    sf_warning ("it seems that velocity is in m/s - will be divided by 1000");	
	    isVelKMS = false;					
	}			
    } else { sf_error ("Need input: velocity model"); }

// Output file
    sembFile_ = sf_output("out");

// Depth/time axis 
    if ( !sf_histint   (dataFile_, "n1", &tNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile_, "d1", &tStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile_, "o1", &tStart_) ) sf_error ("Need o1= in input");
// Dip angle axis 
    if ( !sf_histint   (dataFile_, "n2", &dipNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile_, "d2", &dipStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile_, "o2", &dipStart_) ) sf_error ("Need o2= in input");
// x axis 
    if ( !sf_histint   (dataFile_, "n3", &xNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile_, "d3", &xStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile_, "o3", &xStart_) )   sf_error ("Need o3= in input");

    // time - in ms
    corUnit = (char*) "ms"; unit = sf_histstring (dataFile_, "unit1"); if (!unit) sf_error ("unit1 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { tStep_ *= 1000; tStart_ *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile_, "unit3"); if (!unit) sf_error ("unit3 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { xStep_ *= 1000; xStart_ *= 1000; }

    // VELOCITY MODEL PARAMS

    if ( !sf_histint   (velFile_, "n1", &v_tNum_) )   sf_error ("Need n1= in velocity file");
    if ( !sf_histfloat (velFile_, "d1", &v_tStep_) )  sf_error ("Need d1= in velocity file");
    if ( !sf_histfloat (velFile_, "o1", &v_tStart_) ) sf_error ("Need o1= in velocity file");

    if ( !sf_histint   (velFile_, "n2", &v_xNum_) )   sf_error ("Need n2= in velocity file");
    if ( !sf_histfloat (velFile_, "d2", &v_xStep_) )  sf_error ("Need d2= in velocity file");
    if ( !sf_histfloat (velFile_, "o2", &v_xStart_) ) sf_error ("Need o2= in velocity file");

    corUnit = (char*) "ms"; unit = sf_histstring (velFile_, "unit1"); if (!unit) sf_error ("unit1 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_tStep_ *= 1000; v_tStart_ *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile_, "unit2"); if (!unit) sf_error ("unit2 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { v_xStep_ *= 1000; v_xStart_ *= 1000; }


    if ( !sf_getint ("gn",    &gammaNum_) ) gammaNum_ = 1;
    /* number of scanned Vm/V values  */
    if (!gammaNum_) {sf_warning ("gn value is changed to 1"); gammaNum_ = 1;}

    if ( !sf_getfloat ("go",    &gammaStart_) ) gammaStart_ = 1.0;
    /* start of Vm/V parameter */
    if (!gammaStart_) {sf_warning ("gn value is changed to 1.0"); gammaStart_ = 1.0;}

    if ( !sf_getfloat ("gd",    &gammaStep_) ) gammaStep_ = 1;
    /* increment of Vm/V parameter */
    if (!gammaStep_) {sf_warning ("gd value is changed to 0.01"); gammaStep_ = 0.01;}

    if ( !sf_getint ("coher",   &coher_) )   coher_ = 11;
    /* height of a vertical window for semblance calculation */
    if (!coher_) {sf_warning ("coher value is changed to 1"); coher_ = 1;}

    if ( !sf_getint ("cigNum",   &distNum_) ) distNum_ = 1;
    /* height of a vertical window for semblance calculation */
    if (!coher_) {sf_warning ("cigNum value is changed to 1"); distNum_ = 1;}

    if ( !sf_getfloat ("dlim",   &dlim_) )   dlim_ = fabs (dipStart_);
    /* defines dip-angle-window for the analysis */
    if (!dlim_) {sf_warning ("coher value is changed to 1"); coher_ = fabs (dipStart_);}

    if (!sf_getbool ( "isSemb", &isSemb_) ) isSemb_ = true;
    /* y - output is semblance; n - stack power */

// OUTPUT FILE

    sf_putint (sembFile_, "n1", tNum_); 
    sf_putint (sembFile_, "n2", gammaNum_); 
    sf_putint (sembFile_, "n3", xNum_); 
    sf_putfloat (sembFile_, "d1", tStep_); 
    sf_putfloat (sembFile_, "d2", gammaStep_); 
    sf_putfloat (sembFile_, "d3", xStep_); 
    sf_putfloat (sembFile_, "o1", tStart_); 
    sf_putfloat (sembFile_, "o2", gammaStart_); 
    sf_putfloat (sembFile_, "o3", xStart_); 
    sf_putstring (sembFile_, "label1", "time"); 
    sf_putstring (sembFile_, "label2", "Vm/V");
    sf_putstring (sembFile_, "label3", "inline");
    sf_putstring (sembFile_, "unit1", "s"); 
    sf_putstring (sembFile_, "unit2", ""); 
    sf_putstring (sembFile_, "unit3", "m"); 

    dagSize  = tNum_ * dipNum_;
    sembSize = tNum_ * gammaNum_;
    velSize = v_tNum_ * v_xNum_;

    halfXApp = distNum_ / 2;
    distStart_ = -halfXApp * xStep_;

    pSembPanel_ = sf_floatalloc (sembSize);

    pVel_ = sf_floatalloc (velSize);
    sf_seek (velFile_, 0, SEEK_SET);
    sf_floatread (pVel_, velSize, velFile_);
	
    if (!isVelKMS) {
	ptrVel = pVel_;
	for (iv = 0; iv < velSize; ++iv, ++ptrVel)
	    *ptrVel /= 1000;
    }

    for (ix = 0; ix < xNum_; ++ix) {

	sf_warning ("scanning: CIG %d of %d;", ix + 1, xNum_);

        xNum = distNum_;       
        startInd = ix - halfXApp;

	distShift = 0.f;

	// boundary checking
        temp = ix - halfXApp;
	if (temp < 0) { xNum += temp; startInd -= temp; distShift = -temp * xStep_; }
	temp = xNum_ - (ix + halfXApp) - 1;
	if (temp < 0) xNum += temp;

	// memory allocation
	curPanelSize = dagSize * xNum; 
	pData_ = sf_floatalloc (curPanelSize);

	// read data
	offset = dagSize * startInd * sizeof (float);
	sf_seek (dataFile_, offset, SEEK_SET);
	sf_floatread (pData_, curPanelSize, dataFile_);	

	geoX = xStart_ + ix * xStep_;
	velInd = (geoX - v_xStart_) / v_xStep_;
	velTrace = pVel_ + velInd;

	if (isSemb_)
	    processGatherSemb (velTrace, xNum, distShift);
	else
	    processGatherStack (velTrace, xNum, distShift);

	sf_floatwrite (pSembPanel_, sembSize, sembFile_);
    }

    sf_warning (".");

    free (pSembPanel_);
    free (pData_);		
    free (pVel_);		

    sf_fileclose (dataFile_);
    sf_fileclose (sembFile_);
    sf_fileclose (velFile_);

    return 0;
}
