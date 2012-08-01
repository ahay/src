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
sf_file sembFile_;          // output file - semblance

// data
float* pData_;               //  --"--  --"--  square gathers data 
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

// scan parameters
int   gammaNum_;
float gammaStep_;
float gammaStart_;

int   coher_;
float dlim_;
bool  isSemb_;

// FUNCTIONS

void processGatherStack () {
	
	float* ptrSemb = pSembPanel_;

	for (int ig = 0; ig < gammaNum_; ++ig) {
		const float curGamma = gammaStart_ + ig * gammaStep_;
		for (int it = 0; it < tNum_; ++it, ++ptrSemb) {			
			const float curT = tStart_ + it * tStep_;
			float curStack = 0.f;
			for (int id = 0; id < dipNum_; ++id) {
				const float curDip = dipStart_ + id * dipStep_;
				if (fabs (curDip) > dlim_) continue;
				const float curDipRad = curDip * SF_PI / 180.f;
	
				const float t = curT * cos (curDipRad) / sqrt (1 - pow (curGamma * sin (curDipRad), 2) );
				const int tIndBase = (t - tStart_) / tStep_;

				for (int ic = 0; ic < coher_; ++ic) {
					const int tInd = tIndBase + ic - coher_ / 2;
					if (tInd < 0 || tInd >= tNum_)
						continue;
					const int ind = id * tNum_ + tInd;
					const float val = *(pData_ + ind);

					curStack += val;
				}
			}
			*ptrSemb = fabs (curStack);
		}
	}

	return;
}

void processGatherSemb () {
	
	float* ptrSemb = pSembPanel_;

	float* energyInput  = sf_floatalloc (coher_);
	float* energyOutput = sf_floatalloc (coher_);

	for (int ig = 0; ig < gammaNum_; ++ig) {
		const float curGamma = gammaStart_ + ig * gammaStep_;
		for (int it = 0; it < tNum_; ++it, ++ptrSemb) {			
			const float curT = tStart_ + it * tStep_;
		    memset ( energyInput,  0, coher_ * sizeof (float) );   
    		memset ( energyOutput, 0, coher_ * sizeof (float) );   
			int count = 0;
			for (int id = 0; id < dipNum_; ++id) {
				const float curDip = dipStart_ + id * dipStep_;
				if (fabs (curDip) > dlim_) continue;
				const float curDipRad = curDip * SF_PI / 180.f;
	
				const float t = curT * cos (curDipRad) / sqrt (1 - pow (curGamma * sin (curDipRad), 2) );
				const int tIndBase = (t - tStart_) / tStep_;

				for (int ic = 0; ic < coher_; ++ic) {
					const int tInd = tIndBase + ic - coher_ / 2;
					if (tInd < 0 || tInd >= tNum_)
						continue;
					const int ind = id * tNum_ + tInd;
					const float val = *(pData_ + ind);

					energyOutput [ic] += val;
					energyInput  [ic] += val * val;
				}
				++count;
			}
			float fullInput  = 0.f;
			float fullOutput = 0.f;	
			for (int ic = 0; ic < coher_; ++ic) {
				fullInput  += energyInput [ic];
				fullOutput += pow (energyOutput [ic], 2);
			}
			
			*ptrSemb = fullInput ? fullOutput / (count * fullInput) : 0.f;
		}
	}

	free (energyInput);
	free (energyOutput);

	return;
}

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc, argv);
// Input files
    dataFile_   = sf_input("in");

// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

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

	const int dagSize  = tNum_ * dipNum_;
	const int sembSize = tNum_ * gammaNum_;

	pData_      = sf_floatalloc (dagSize);
	pSembPanel_ = sf_floatalloc (sembSize);

	for (int ix = 0; ix < xNum_; ++ix) {
		sf_warning ("CIG %d of %d;", ix + 1, xNum_);

		const size_t startPos = (size_t) ix * dagSize * sizeof (float);
		sf_seek (dataFile_, startPos, SEEK_SET);
		sf_floatread (pData_, dagSize, dataFile_);

		if (isSemb_)
			processGatherSemb ();
		else
			processGatherStack ();

		sf_floatwrite (pSembPanel_, sembSize, sembFile_);
	}

	sf_warning (".");

	free (pSembPanel_);
	free (pData_);		

	sf_fileclose (dataFile_);
	sf_fileclose (sembFile_);

	return 0;
}
