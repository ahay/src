/* CRS-based semblance for 3D
Several CIGs are used simultaneously. Dip-angle sections corresponding to the same 
dip-angle compose a subvolume. The subvolume allows calculating semblance in the
scattering-angle direction along reflection boundaries.

Input:
	inDags_.rsf   - 3D dip-angle gathers - stack in the scattering-angle direction
	inDagsSq_.rsf - stack of amplitude squares in the scattering-angle direction

	Working with just dip-angle gathers use default value of "scatnum" parameter

Output:
	sembFile_.rsf - crs-based semblance file; has the same dimensions as the input files
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
#ifdef _OPENMP
#include <omp.h>
#endif

// VARIABLES

// files
sf_file inDags_;            // dip-angle gathers file - stack in the scattering-angle direction
sf_file inDagsSq_;          // stack of amplitude squares in the scattering-angle direction
sf_file sembFile_;          // output file - crs-base semblance

// data
float* ptrToDags_;          // pointer to the dip-angle gathers data
float* ptrToDagsSq_;        //  --"--  --"--  square gathers data 
float* ptrToSembPanel_;     //  --"--  --"--  result - crs-based semblance

float* ptrToData_;          // pointer to the dip-angle gathers data
float* ptrToDataSq_;        //  --"--  --"--  square gathers data 


// dip-angle gathers dimensions:
int   zNum_;                 
float zStart_;
float zStep_;

int   dipxNum_;
float dipxStart_;
float dipxStep_;

int   dipyNum_;
float dipyStart_;
float dipyStep_;

int   xNum_;
float xStart_;
float xStep_;

int   yNum_;
float yStart_;
float yStep_;

int dagSize_;               // size (in samples) of one dip-angle gather
int scatnum_;               // shows how many traces were stacked in the scattering angle direction
int coher_;                 // height (in samples) of a vertical window for semblance calculation
int halfCoher_;             // half of the "coher"

int xapp_;                  // number of CIGs in the inline-direction processed simultaneously
int halfXapp_;              // half of the "xapp"
int xdipapp_;               // number of traces in the x-dip direction processed simultaneously

int yapp_;                  // number of CIGs in the crossline-direction processed simultaneously
int halfYapp_;              // half of the "yapp"
int ydipapp_;               // number of traces in the y-dip direction processed simultaneously

// FUNCTIONS

// Check if the point with its surounding is inside data volume;
// if not - correct analysis aperture
void checkBoundary (int* yPos, int* curYapp, int* xPos, int* curXapp) {

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

    // checking in Y-dimension

    // bottom edge
    if (*yPos < 0) {
        int diff = -1 * *yPos;
		*curYapp -= diff;
		*yPos = 0;
	}

    // top edge
    if (*yPos + *curYapp - 1 >= yNum_) {
        int diff = *yPos + *curYapp - yNum_; 
        *curYapp -= diff;
    }

	return;
}

void readBlockAroundPoint (int yPos, int xPos, int halfYapp, int halfXapp, int* curYapp, int*bottomShift, 
						   int* curXapp, int* leftShift) {

    int startY = yPos - halfYapp;
    int startX = xPos - halfXapp;
    
    // check if the apperture is adequate... if not - correct it and the start point
    checkBoundary (&startY, curYapp, &startX, curXapp);
	*leftShift = xPos - startX;
	*bottomShift = yPos - startY;

	const int volSize = dagSize_ * (*curYapp) * (*curXapp);
	const int pointsNumToRead = dagSize_ * (*curXapp);

	ptrToDags_   = sf_floatalloc (volSize);
	ptrToDagsSq_ = sf_floatalloc (volSize);

	memset (ptrToDags_,   0, volSize * sizeof (float));
	memset (ptrToDagsSq_, 0, volSize * sizeof (float));

	ptrToData_   = sf_floatalloc (volSize);
	ptrToDataSq_ = sf_floatalloc (volSize);
	memset (ptrToData_,   0, volSize * sizeof (float));
	memset (ptrToDataSq_, 0, volSize * sizeof (float));			

	for (int iy = 0; iy < *curYapp; ++iy) {
		
		const size_t startPos = (size_t) ((iy + startY) * xNum_ + startX) * dagSize_ * sizeof (float);
			
		sf_seek (inDags_,   startPos, SEEK_SET);
		sf_seek (inDagsSq_, startPos, SEEK_SET);

		const size_t shift = iy * pointsNumToRead;

		sf_floatread (ptrToData_   + shift, pointsNumToRead, inDags_);
		sf_floatread (ptrToDataSq_ + shift, pointsNumToRead, inDagsSq_);
	}

	// substacking in the y-dip and x-dip directions
	for (int iy = 0; iy < *curYapp; ++iy) {
		const int yShift = iy * (*curXapp) * dipyNum_ * dipxNum_ * zNum_;
		for (int ix = 0; ix < *curXapp; ++ix) {
			const int xShift = ix * dipyNum_ * dipxNum_ * zNum_;
			const size_t shiftToGather = yShift + xShift;

			for (int idy = 0; idy < dipyNum_; ++idy) {
				for (int idx = 0; idx < dipxNum_; ++idx) {
					const size_t tracesShift = shiftToGather + (idy * dipxNum_ + idx) * zNum_;
		
					float* ptrToTrace   = ptrToDags_   + tracesShift;
					float* ptrToTraceSq = ptrToDagsSq_ + tracesShift;

					for (int idya = 0; idya < ydipapp_; ++idya) {		
						int indy = idy - ydipapp_ / 2 + idya;
						if (indy < 0 || indy >= dipyNum_) continue;		
						for (int idxa = 0; idxa < xdipapp_; ++idxa) {		
							int indx = idx - xdipapp_ / 2 + idxa;
							if (indx < 0 || indx >= dipxNum_) continue;		
			
							const size_t dataShift = shiftToGather + (indy * dipxNum_ + indx) * zNum_;
							float* ptrFrom   = ptrToData_   + dataShift;
							float* ptrSqFrom = ptrToDataSq_ + dataShift;

							float* ptrTo     = ptrToTrace;
							float* ptrSqTo   = ptrToTraceSq;

							for (int iz = 0; iz < zNum_; ++iz, ++ptrTo, ++ptrSqTo, ++ptrFrom, ++ptrSqFrom) {
								*ptrTo += *ptrFrom;
								*ptrSqTo += *ptrSqFrom;
							}
						}
					}
				}
			}
		}
	}

	return;
}

void getSemblanceForTrace (int gathersNum, float* stack, float* stackSq, float* semb) {
   
    double sumOutput, sumInput;
    int im, it;

    const int targetZNum     = zNum_; 
	const int zNumFull    = 2 * halfCoher_ + zNum_;

	float stackVal   = 0.f;
	float stackSqVal = 0.f;

	float* traceSumOutput = sf_floatalloc (zNumFull);
	float* traceSumInput  = sf_floatalloc (zNumFull);
    memset (traceSumOutput, 0, zNumFull * sizeof (float));   
    memset (traceSumInput,  0, zNumFull * sizeof (float));   

#ifdef _OPENMP 
#pragma omp parallel for
#endif
	for (it = 0; it < targetZNum; ++it) {
		const int trInd = it + halfCoher_;
        for (im = 0; im < gathersNum; ++im){
			const int ind = im * zNum_ + it;
			stackVal   = stack[ind];
			stackSqVal = stackSq[ind];

		    traceSumOutput[trInd] += stackVal;
		    traceSumInput[trInd]  += stackSqVal;
		}
		traceSumOutput[trInd] *= traceSumOutput[trInd];
	}
#ifdef _OPENMP 
#pragma omp parallel for
#endif
    for (int it = 0; it < targetZNum; ++it) {
        sumOutput = 0.f;
        sumInput  = 0.f;
		const int temp = it + coher_;
		for (int j = it; j < temp; ++j) {
		    sumOutput += traceSumOutput[j];
		    sumInput  += traceSumInput[j];
		}
		semb[it] = sumInput ? sumOutput / (gathersNum * scatnum_* sumInput) : 0.f;
    }

    free (traceSumOutput);
    free (traceSumInput);

    return;
}

void buildFilter (int curyapp, int bottomShift, int curxapp, int leftShift, float* ptrToSembPanel) {

	const int fullSampNumber = 2 * halfCoher_ + zNum_;
	const int dipNum = dipyNum_ * dipxNum_;
	const int dagSize = dipNum * zNum_;
	const int stackGridSize = curyapp * curxapp * fullSampNumber;
	const float CONVRATIO = M_PI / 180.f;
	
	float* ptrToSembTrace = ptrToSembPanel;

	for (int idy = 0; idy < dipyNum_; ++idy) {
		const float curyDip = dipyStart_ + idy * dipyStep_;
		const float curyDipInRad = curyDip * CONVRATIO;
		const float tanyDipInRad = tan (curyDipInRad);
		for (int idx = 0; idx < dipxNum_; ++idx, ptrToSembTrace += zNum_) {
			const float curxDip = dipxStart_ + idx * dipxStep_;
			const float curxDipInRad = curxDip * CONVRATIO;
			const float tanxDipInRad = tan (curxDipInRad);
	
			float* stackGrid   = sf_floatalloc (stackGridSize);
			float* stackSqGrid = sf_floatalloc (stackGridSize);
			memset (stackGrid,   0, stackGridSize * sizeof (float));
			memset (stackSqGrid, 0, stackGridSize * sizeof (float));

			const size_t shiftInDag = (idy * dipxNum_ + idx) * zNum_;
		
		    int tracesNum = 0;
			for (int iy = 0; iy < curyapp; ++iy) {					
				const float yShift = (iy - bottomShift) * yStep_;
	   		    const float yDepthShift = yShift * tanyDipInRad;
				const int yDepthShiftSamp = yDepthShift / zStep_;
				for (int ix = 0; ix < curxapp; ++ix) {		
					const float xShift = (ix - leftShift) * xStep_;
		   		    const float xDepthShift = xShift * tanxDipInRad;
					const int xDepthShiftSamp = xDepthShift / zStep_;

					const size_t dataShift  = shiftInDag + (iy * curxapp + ix) * dagSize;

					float* ptrDataStack   = ptrToDags_   + dataShift;
					float* ptrDataStackSq = ptrToDagsSq_ + dataShift;

					const int stackShift = tracesNum * zNum_;
		
					float* ptrStackGrid   = stackGrid + stackShift;
					float* ptrStackSqGrid = stackSqGrid + stackShift;

					int zInd = -yDepthShiftSamp - xDepthShiftSamp;
				
					for (int iz = 0; iz < zNum_; ++iz, ++ptrStackGrid, ++ptrStackSqGrid, ++zInd) {
						if (zInd < 0 || zInd >= zNum_) continue;
						*ptrStackGrid = *(ptrDataStack + zInd);
						*ptrStackSqGrid = *(ptrDataStackSq + zInd);
					}		
					++tracesNum;
				}
			}
			getSemblanceForTrace (tracesNum, stackGrid, stackSqGrid, ptrToSembTrace);  

			free (stackGrid);
			free (stackSqGrid);
		}
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
    if ( !sf_histint   (inDags_, "n1", &zNum_) )      sf_error ("Need n1= in input");
    if ( !sf_histfloat (inDags_, "d1", &zStep_) )     sf_error ("Need d1= in input");
    if ( !sf_histfloat (inDags_, "o1", &zStart_) )    sf_error ("Need o1= in input");
// x-dip angle axis 
    if ( !sf_histint   (inDags_, "n2", &dipxNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (inDags_, "d2", &dipxStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (inDags_, "o2", &dipxStart_) ) sf_error ("Need o2= in input");
// y-dip angle axis 
    if ( !sf_histint   (inDags_, "n3", &dipyNum_) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (inDags_, "d3", &dipyStep_) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (inDags_, "o3", &dipyStart_) ) sf_error ("Need o3= in input");
// x axis 
    if ( !sf_histint   (inDags_, "n4", &xNum_) )      sf_error ("Need n4= in input");
    if ( !sf_histfloat (inDags_, "d4", &xStep_) )     sf_error ("Need d4= in input");
    if ( !sf_histfloat (inDags_, "o4", &xStart_) )    sf_error ("Need o4= in input");
// y axis 
    if ( !sf_histint   (inDags_, "n5", &yNum_) )      sf_error ("Need n5= in input");
    if ( !sf_histfloat (inDags_, "d5", &yStep_) )     sf_error ("Need d5= in input");
    if ( !sf_histfloat (inDags_, "o5", &yStart_) )    sf_error ("Need o5= in input");
	
    if ( !sf_getint ("xapp",    &xapp_) )    xapp_ = 1;
    /* number of CIGs in the inline-direction processed simultaneously */
	if (!xapp_) {sf_warning ("xapp value is changed to 1"); xapp_ = 1;}

    if ( !sf_getint ("yapp",    &yapp_) )    yapp_ = 1;
    /* number of CIGs in the crossline-direction processed simultaneously */
	if (!yapp_) {sf_warning ("yapp value is changed to 1"); yapp_ = 1;}

    if ( !sf_getint ("dipappx",    &xdipapp_) ) xdipapp_ = 11;
    /* number of traces in the x-dip direction processed simultaneously */
	if (!xdipapp_) {sf_warning ("dipapp value is changed to 11"); xdipapp_ = 11;}

    if ( !sf_getint ("dipappy",    &ydipapp_) ) ydipapp_ = 11;
    /* number of traces in the y-dip direction processed simultaneously */
	if (!xdipapp_) {sf_warning ("dipappy value is changed to 11"); ydipapp_ = 11;}

    if ( !sf_getint ("coher",   &coher_) ) coher_ = 11;
	/* height of a vertical window for semblance calculation */
	if (!coher_) {sf_warning ("coher value is changed to 11"); coher_ = 11;}

    if ( !sf_getint ("scatnum", &scatnum_) ) scatnum_ = 1;
	/* shows how many traces were stacked in the scattering angle direction; 
	   if the stack was normalized use the default value 
	*/ 

	dagSize_ = zNum_ * dipxNum_ * dipyNum_;
	halfCoher_ = coher_ / 2;    // yes - this is the integer division	
	halfXapp_  = xapp_ / 2;     // this is the integer division too	
	halfYapp_  = yapp_ / 2;     // this is the integer division too	

	const int cigNum = xNum_ * yNum_;
	int cigInd = 1;

	for (int iy = 0; iy < yNum_; ++iy) {
		for (int ix = 0; ix < xNum_; ++ix, ++cigInd) {
		
			sf_warning ("CIG %d of %d;", cigInd, cigNum);
		
			// xapp for the currect core CIG; it can be changed by the checkBoundary ()
			int curxapp = xapp_; 
	        // distance between the core gather and the left side of the aperture
			int leftShift = 0;	
			// yapp for the currect core CIG; it can be changed by the checkBoundary ()
			int curyapp = yapp_; 
	        // distance between the core gather and the bottom side of the aperture
			int bottomShift = 0;	

			ptrToSembPanel_ = sf_floatalloc (dagSize_);
			memset (ptrToSembPanel_, 0, dagSize_ * sizeof (float));

			readBlockAroundPoint (iy, ix, halfYapp_, halfXapp_, &curyapp, &bottomShift, &curxapp, &leftShift);
			buildFilter (curyapp, bottomShift, curxapp, leftShift, ptrToSembPanel_);

			free (ptrToDags_);
			free (ptrToDagsSq_);

			free (ptrToData_);
			free (ptrToDataSq_);

			sf_floatwrite (ptrToSembPanel_, dagSize_, sembFile_);
			free (ptrToSembPanel_);
		}
	}

	sf_warning (".");

	sf_fileclose (inDags_);
	sf_fileclose (inDagsSq_);
	sf_fileclose (sembFile_);

	return 0;
}
