/* CRS-based semblance
Several CIGs are used simultaneously. Dip-angle sections corresponding to the same 
dip-angle compose a subvolume. The subvolume allows calculating semblance in the
scattering-angle direction along reflection boundaries.

Input:
	inDags_.rsf   - dip-angle gathers - stack in the scattering-angle direction
	InDagsSq_.rsf - stack of amplitude squares in the scattering-angle direction

Output:
	sembFile_.rsf - crs-based semblance file; has the same dimensions as the input files
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
sf_file inDags_;            // dip-angle gathers file - stack in the scattering-angle direction
sf_file inDagsSq_;          // stack of amplitude squares in the scattering-angle direction
sf_file sembFile_;          // output file - crs-base semblance

// data
float* ptrToDags_;          // pointer to the dip-angle gathers data
float* ptrToDagsSq_;        //  --"--  --"--  square gathers data 
float* ptrToSembPanel_;     //  --"--  --"--  result - crs-based semblance

// dip-angle gathers dimensions:
int zNum_;                 
float zStart_;
float zStep_;

int dipNum_;
float dipStart_;
float dipStep_;

int xNum_;
float xStart_;
float xStep_;

int dagSize_;               // size (in samples) of one dip-angle gather
int scatnum_;               // shows how many traces were stacked in the scattering angle direction
int coher_;                 // height (in samples) of a vertical window for semblance calculation
int halfCoher_;             // half of the "coher"
int xapp_;                  // number of CIGs in the inline-direction processed simultaneously
int halfXapp_;              // half of the "xapp"

// FUNCTIONS

// Check if the point with its surounding is inside data volume;
// if not - correct analysis aperture
void checkBoundary (int* xPos, int* curXapp) {

    // checking in X-dimension

    // left edge
    if (*xPos < 0) {
        int diff = -1 * *xPos;
		*curXapp -= diff;
		*xPos = 0;
	}

    // right edge
    if (*xPos + *curXapp - 1 >= xNum_) {
        int diff = *xPos + *curXapp - xNum_; 
        *curXapp -= diff;
    }

	return;
}

void readBlockAroundPoint (int xPos, int halfXapp, int* curXapp, int* leftShift) {

    int startX = xPos - halfXapp;
    
    // check if the apperture is adequate... if not - correct it and the start point
    checkBoundary (&startX, curXapp);
	*leftShift = xPos - startX;
	const int pointsNumToRead = dagSize_ * (*curXapp);
	ptrToDags_   = sf_floatalloc (pointsNumToRead);
	ptrToDagsSq_ = sf_floatalloc (pointsNumToRead);
	memset (ptrToDags_,   0, pointsNumToRead * sizeof (float));
	memset (ptrToDagsSq_, 0, pointsNumToRead * sizeof (float));

	const int startPos = startX * dagSize_ * sizeof (float);

	sf_seek (inDags_,   startPos, SEEK_SET);
	sf_seek (inDagsSq_, startPos, SEEK_SET);

	sf_floatread (ptrToDags_,   pointsNumToRead, inDags_);
	sf_floatread (ptrToDagsSq_, pointsNumToRead, inDagsSq_);
	
	return;
}

void getSemblanceForTrace (int gathersNum, int tracesNum, float* stack, float* stackSq, float* semb) {
   
    double sumOutput, sumInput;
    int im, it, j;

    const int targetZNum     = zNum_; 
	const int fullZNumber    = halfCoher_ + zNum_;

	float stackVal   = 0.f;
	float stackSqVal = 0.f;

	float* traceSumOutput = sf_floatalloc (fullZNumber);
	float* traceSumInput  = sf_floatalloc (fullZNumber);

	for (it = 0; it < fullZNumber; ++it) {
		traceSumOutput[it] = traceSumInput[it] = 0.f;
        for (im = 0; im < gathersNum; ++im){
			const int ind = im * zNum_ + it;
			stackVal   = stack[ind];
			stackSqVal = stackSq[ind];

		    traceSumOutput[it] += stackVal;
		    traceSumInput[it]  += stackSqVal;
		}
		traceSumOutput[it] *= traceSumOutput[it];
	}
 
    for (int ind = 0, it = 0; ind < targetZNum; ++it, ++ind) {
        sumOutput = 0.f;
        sumInput  = 0.f;
		const int temp = it + coher_;
		for (j = it; j < temp; ++j) {
		    sumOutput += traceSumOutput[j];
		    sumInput  += traceSumInput[j];
		}
		semb[ind] = sumInput ? sumOutput / (tracesNum * sumInput) : 0.f;
    }

    free (traceSumOutput);
    free (traceSumInput);

    return;
}

void buildFilter (int curxapp, int leftShift, float* ptrToSembPanel) {

	const int fullSampNumber = 2 * halfCoher_ + zNum_;
	const int stackGridSize = curxapp * fullSampNumber;
	const float CONVRATIO = M_PI / 180.f;
	
	float* ptrToSembTrace = ptrToSembPanel;

	for (int id = 0; id < dipNum_; ++id, ptrToSembTrace += zNum_) {

		const float curDip = dipStart_ + id * dipStep_;
		const float curDipInRad = curDip * CONVRATIO;
		const float tanDipInRad = tan (curDipInRad);
	
		float* stackGrid   = sf_floatalloc (stackGridSize);
		float* stackSqGrid = sf_floatalloc (stackGridSize);
		
	    int tracesNum = 0;
			
		for (int ix = 0; ix < curxapp; ++ix) {		
			const float xShift = (ix - leftShift) * xStep_;
   		    const float depthShift = xShift * tanDipInRad;
			const int depthShiftSamp = depthShift / zStep_;
			const int dataShift = (id + ix * dipNum_) * zNum_;
		
			float* ptrDataStack   = ptrToDags_ + dataShift;
			float* ptrDataStackSq = ptrToDagsSq_ + dataShift;

			const int stackShift = tracesNum * zNum_ + halfCoher_;

			float* ptrStackGrid   = stackGrid + stackShift;
			float* ptrStackSqGrid = stackSqGrid + stackShift;

			int zInd = -depthShiftSamp;
				
			for (int iz = 0; iz < zNum_; ++iz, ++ptrStackGrid, ++ptrStackSqGrid, ++zInd) {
				if (zInd < 0 || zInd >= zNum_) continue;
				*ptrStackGrid = *(ptrDataStack + zInd);
				*ptrStackSqGrid = *(ptrDataStackSq + zInd);
			}		
			++tracesNum;
		}

		getSemblanceForTrace (curxapp, curxapp * scatnum_, stackGrid, stackSqGrid, ptrToSembTrace);  

		free (stackGrid);
		free (stackSqGrid);
	}

	return;
}

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc,argv);
// Input files
    inDags_   = sf_input("in");
    inDagsSq_ = sf_input("dataSq");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (inDags_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */
    if ( SF_FLOAT != sf_gettype (inDagsSq_) ) sf_error ("Need float input: dip-angle gathers in squares");
    /* stacks of amplitude squares in the scattering-angle direction */

// Output file
    sembFile_ = sf_output("out");

// Depth/time axis 
    if ( !sf_histint   (inDags_, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (inDags_, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (inDags_, "o1", &zStart_) ) sf_error ("Need o1= in input");
// Dip angle axis 
    if ( !sf_histint   (inDags_, "n2", &dipNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (inDags_, "d2", &dipStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (inDags_, "o2", &dipStart_) ) sf_error ("Need o2= in input");
// x axis 
    if ( !sf_histint   (inDags_, "n3", &xNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (inDags_, "d3", &xStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (inDags_, "o3", &xStart_) )   sf_error ("Need o3= in input");
	
    if ( !sf_getint ("xapp",    &xapp_) )    xapp_ = 1;
    /* number of CIGs in the inline-direction processed simultaneously */
	if (!xapp_) {sf_warning ("xapp value is changed to 1"); xapp_ = 1;}

    if ( !sf_getint ("coher",   &coher_) )   coher_ = 11;
	/* height of a vertical window for semblance calculation */
	if (!coher_) {sf_warning ("coher value is changed to 1"); coher_ = 1;}

    if ( !sf_getint ("scatnum", &scatnum_) ) scatnum_ = 1;
	/* shows how many traces were stacked in the scattering angle direction; 
	   if the stack was normalized use the default value 
	*/ 

	dagSize_ = zNum_ * dipNum_;
	halfCoher_ = coher_ / 2;    // yes - this is the integer division	
	halfXapp_  = xapp_ / 2;     // this is the integer division too	

	for (int ix = 0; ix < xNum_; ++ix) {
		
		sf_warning ("CIG %d of %d;", ix + 1, xNum_);
		
		// xapp for the currect core CIG; it can be changed by the checkBoundary ()
		int curxapp = xapp_; 
        // distance between the core gather and the left side of the aperture
		int leftShift = 0;	

		ptrToSembPanel_ = sf_floatalloc (dagSize_);
		memset (ptrToSembPanel_, 0, dagSize_ * sizeof (float));

		readBlockAroundPoint (ix, halfXapp_, &curxapp, &leftShift);
		buildFilter (curxapp, leftShift, ptrToSembPanel_);

		free (ptrToDags_);
		free (ptrToDagsSq_);

		sf_floatwrite (ptrToSembPanel_, dagSize_, sembFile_);
		
		free (ptrToSembPanel_);
	}

	sf_warning (".");

	sf_fileclose (inDags_);
	sf_fileclose (inDagsSq_);
	sf_fileclose (sembFile_);

	return 0;
}
