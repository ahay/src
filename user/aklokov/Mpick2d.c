/* 2D picking

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
sf_file paramFile_;           // output file - semblance
sf_file outFile_;           // output file - semblance

// data
float* pData_;               //  --"--  --"--  square gathers data 
float* panel_;          //  --"--  --"--  result - crs-based semblance
float* maxTrace_;       

// dip-angle gathers dimensions:
int   tNum_;                 
float tStart_;
float tStep_;

int   paramNum_;
float paramStart_;
float paramStep_;

int   xNum_;
float xStart_;
float xStep_;

// scan parameters
int   gammaNum_;
float gammaStep_;
float gammaStart_;

int   coher_;
float eps_;

// FUNCTIONS

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc, argv);
// Input files
    dataFile_   = sf_input("in");

// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */

// Output file
	char* paramTag = "maxForPicking.rsf";
    paramFile_ = sf_output (paramTag);
	char* sembTag = "sembForPicking.rsf";
	sembFile_ = sf_output (sembTag);
    outFile_ = sf_output ("out");

// Depth/time axis 
    if ( !sf_histint   (dataFile_, "n1", &tNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile_, "d1", &tStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile_, "o1", &tStart_) ) sf_error ("Need o1= in input");
// Dip angle axis 
    if ( !sf_histint   (dataFile_, "n2", &paramNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile_, "d2", &paramStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile_, "o2", &paramStart_) ) sf_error ("Need o2= in input");
// x axis 
    if ( !sf_histint   (dataFile_, "n3", &xNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile_, "d3", &xStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile_, "o3", &xStart_) )   sf_error ("Need o3= in input");

    if ( !sf_getfloat ("eps", &eps_) ) eps_ = 1;
    /* epsilon  */

// OUTPUT FILE

    sf_putint    (sembFile_, "n1", tNum_);     	 sf_putint    (sembFile_, "n2", xNum_); 
    sf_putfloat  (sembFile_, "d1", tStep_); 	 sf_putfloat  (sembFile_, "d2", xStep_); 
    sf_putfloat  (sembFile_, "o1", tStart_);	 sf_putfloat  (sembFile_, "o2", xStart_); 
    sf_putstring (sembFile_, "label1", "time");  sf_putstring (sembFile_, "label2", "inline");
    sf_putstring (sembFile_, "unit1", "s");      sf_putstring (sembFile_, "unit2", "m"); 

    sf_putint    (paramFile_, "n1", tNum_);     	 sf_putint    (paramFile_, "n2", xNum_); 
    sf_putfloat  (paramFile_, "d1", tStep_);   	 sf_putfloat  (paramFile_, "d2", xStep_); 
    sf_putfloat  (paramFile_, "o1", tStart_);	     sf_putfloat  (paramFile_, "o2", xStart_); 
    sf_putstring (paramFile_, "label1", "time");   sf_putstring (paramFile_, "label2", "inline");
    sf_putstring (paramFile_, "unit1", "s");       sf_putstring (paramFile_, "unit2", "m"); 


	const int panelSize = tNum_ * paramNum_;	
	panel_ = sf_floatalloc (panelSize);
	memset ( panel_, 0, panelSize * sizeof (float) );   

	float* maxTrace = sf_floatalloc (tNum_);
	float* paramTrace = sf_floatalloc (tNum_);

	// compute max

	for (int ix = 0; ix < xNum_; ++ix) {
		sf_warning ("CIG %d of %d;", ix + 1, xNum_);

		const size_t startPos = (size_t) ix * panelSize * sizeof (float);
		sf_seek (dataFile_, startPos, SEEK_SET);
		sf_floatread (panel_, panelSize, dataFile_);

		memset ( maxTrace, 0, tNum_ * sizeof (float) );   
	
		float* pPanel = panel_;

		for (int ip = 0; ip < paramNum_; ++ip) {
			float* pMax   = maxTrace;			
			float* pParam = paramTrace;			
			const float curParam = paramStart_ + ip * paramStep_;
			for (int it = 0; it < tNum_; ++it, ++pMax, ++pPanel, ++pParam) {
				if (*pMax < *pPanel) { *pMax = *pPanel; *pParam = curParam; }							
			}
		}
		sf_floatwrite (maxTrace,   tNum_, sembFile_);
		sf_floatwrite (paramTrace, tNum_, paramFile_);
	}

	sf_warning (".");

	free (pData_);
	free (panel_);
	free (maxTrace);		
	free (paramTrace);		

	sf_fileclose (dataFile_);
	sf_fileclose (paramFile_);
	sf_fileclose (sembFile_);

	// regularization

	sembFile_  = sf_input (sembTag);
	paramFile_ = sf_input (paramTag);

	for (int ix = 0; ix < xNum_; ++ix) {
		sf_warning ("CIG %d of %d;", ix + 1, xNum_);

/*        xNum = xApp;       
        startInd = it - halfXApp;

		// boundary checking
        int temp = it - halfXApp;
		if (temp < 0) {
	    	xNum += temp;
	    	startInd -= temp;
		}	
        temp = tracesNum - (it + halfXApp) - 1;
		if (temp < 0)
		    xNum += temp;
*/
	}

	sf_fileclose (paramFile_);
	sf_fileclose (sembFile_);
	sf_fileclose (outFile_);

	return 0;
}
