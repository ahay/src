/* CRS-based semblance

   Several CIGs are used simultaneously. Dip-angle sections corresponding to the same 
   dip-angle compose a subvolume. The subvolume allows calculating semblance in the
   scattering-angle direction along reflection boundaries.

   Input:
   inDags_.rsf   - dip-angle gathers - stack in the scattering-angle direction
   InDagsSq_.rsf - stack of amplitude squares in the scattering-angle direction

   Working with just dip-angle gathers use default value of "scatnum" parameter

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

/* VARIABLES */

/* files */
sf_file inDags_;            /* dip-angle gathers file - stack in the scattering-angle direction */
sf_file inDagsSq_;          /* stack of amplitude squares in the scattering-angle direction */
sf_file sembFile_;          /* output file - crs-base semblance */

/* data */
float* ptrToDags_;          /* pointer to the dip-angle gathers data */
float* ptrToDagsSq_;        /*  --"--  --"--  square gathers data */
float* ptrToSembPanel_;     /*  --"--  --"--  result - crs-based semblance */

float* ptrToData_;          /* pointer to the dip-angle gathers data */
float* ptrToDataSq_;        /*  --"--  --"--  square gathers data */


/* dip-angle gathers dimensions: */
int zNum_;                 
float zStart_;
float zStep_;

int dipNum_;
float dipStart_;
float dipStep_;

int xNum_;
float xStart_;
float xStep_;

int dagSize_;               /* size (in samples) of one dip-angle gather */
int scatnum_;               /* shows how many traces were stacked in the scattering angle direction */
int coher_;                 /* height (in samples) of a vertical window for semblance calculation */
int halfCoher_;             /* half of the "coher" */
int xapp_;                  /* number of CIGs in the inline-direction processed simultaneously */
int halfXapp_;              /* half of the "xapp" */
int xdipapp_;               /* number of traces in the x-dip direction processed simultaneously */

bool makeWeight_;
float s1_;
float s2_;
float ds_;

/* FUNCTIONS */

/* Check if the point with its surounding is inside data volume; */
/* if not - correct analysis aperture */
void checkBoundary (int* xPos, int* curXapp) 
{
    int diff;

    /* checking in X-dimension */

    /* left edge */
    if (*xPos < 0) {
        diff = -1 * *xPos;
		*curXapp -= diff;
		*xPos = 0;
    }

    /* right edge */
    if (*xPos + *curXapp - 1 >= xNum_) {
        diff = *xPos + *curXapp - xNum_; 
        *curXapp -= diff;
    }

    return;
}

void readBlockAroundPoint (int xPos, int halfXapp, int* curXapp, int* leftShift) {

    int startX = xPos - halfXapp;
    int ix, id, ida, iz, ind, dataShift, tracesShift, pointsNumToRead;
    size_t startPos;
    float *ptrToTrace, *ptrToTraceSq, *ptrFrom, *ptrTo, *ptrSqFrom, *ptrSqTo;

    /* check if the apperture is adequate... if not - correct it and the start point */
    checkBoundary (&startX, curXapp);
    *leftShift = xPos - startX;
    pointsNumToRead = dagSize_ * (*curXapp);

    ptrToDags_   = sf_floatalloc (pointsNumToRead);
    ptrToDagsSq_ = sf_floatalloc (pointsNumToRead);
    memset (ptrToDags_,   0, pointsNumToRead * sizeof (float));
    memset (ptrToDagsSq_, 0, pointsNumToRead * sizeof (float));

    ptrToData_   = sf_floatalloc (pointsNumToRead);
    ptrToDataSq_ = sf_floatalloc (pointsNumToRead);
    memset (ptrToData_,   0, pointsNumToRead * sizeof (float));
    memset (ptrToDataSq_, 0, pointsNumToRead * sizeof (float));			

    startPos = (size_t) startX * dagSize_ * sizeof (float);

    sf_seek (inDags_,   startPos, SEEK_SET);
    sf_seek (inDagsSq_, startPos, SEEK_SET);

    sf_floatread (ptrToData_,   pointsNumToRead, inDags_);
    sf_floatread (ptrToDataSq_, pointsNumToRead, inDagsSq_);

    /* substacking in the x-dip direction */
		
    for (ix = 0; ix < *curXapp; ++ix) {
		for (id = 0; id < dipNum_; ++id) {
			
		    tracesShift = ix * dipNum_ + id;
		    ptrToTrace   = ptrToDags_ + tracesShift * zNum_;
		    ptrToTraceSq = ptrToDagsSq_ + tracesShift * zNum_;

		    for (ida = 0; ida < xdipapp_; ++ida) {		
				ind = id - xdipapp_ / 2 + ida;
				if (ind < 0 || ind >= dipNum_) continue;		
				
				dataShift = ix * dipNum_ + ind;
				ptrFrom   = ptrToData_   + dataShift * zNum_;
				ptrSqFrom = ptrToDataSq_ + dataShift * zNum_;

				ptrTo     = ptrToTrace;
				ptrSqTo   = ptrToTraceSq;

				for (iz = 0; iz < zNum_; ++iz, ++ptrTo, ++ptrSqTo, ++ptrFrom, ++ptrSqFrom) {
				    *ptrTo += *ptrFrom;
				    *ptrSqTo += *ptrSqFrom;
				}
	    	}
		}
    }

    return;
}

void getSemblanceForTrace (int tracesNum, float* stack, float* stackSq, float* semb) {
   
    double sumOutput, sumInput;
    int im, it, temp, j, trInd, ind;

    const int targetZNum = zNum_; 
    const int zNumFull   = 2 * halfCoher_ + zNum_;

    const int k = tracesNum * scatnum_ * xdipapp_;

    float stackVal   = 0.f;
    float stackSqVal = 0.f;

    float sembval;

    float* traceSumOutput; 
    float* traceSumInput;  

    traceSumOutput = sf_floatalloc (zNumFull);
    traceSumInput  = sf_floatalloc (zNumFull);

    memset (traceSumOutput, 0, zNumFull * sizeof (float));   
    memset (traceSumInput,  0, zNumFull * sizeof (float));   

    for (it = 0; it < targetZNum; ++it) {
		trInd = it + halfCoher_;
        for (im = 0; im < tracesNum; ++im){
		    ind = im * zNum_ + it;
		    stackVal   = stack[ind];
		    stackSqVal = stackSq[ind];
	
		    traceSumOutput[trInd] += stackVal;
		    traceSumInput[trInd]  += stackSqVal;
		}
		traceSumOutput[trInd] *= traceSumOutput[trInd];
    }
 
    for (it = 0; it < targetZNum; ++it) {
        sumOutput = 0.f;
        sumInput  = 0.f;
		temp = it + coher_;
		for (j = it; j < temp; ++j) {
		    sumOutput += traceSumOutput[j];
		    sumInput  += traceSumInput[j];
		}

		sembval = sumInput ? sumOutput / (k * sumInput) : 0.f;

		if (makeWeight_) {
			if (sembval > s2_) semb [it] = 1.0;
			else if (sembval < s1_) semb [it] = 0.0;
			else semb[it] = (sembval - s1_) / ds_;
		} else {
			semb[it] = sembval;
		}
    }

    free (traceSumOutput);
    free (traceSumInput);

    return;
}

void buildFilter (int curxapp, int leftShift, float* ptrToSembPanel) {

    int id, tracesNum, ix, depthShiftSamp, dataShift, stackShift;
    int zInd, iz;
    float curDip, curDipInRad, tanDipInRad, xShift,  depthShift;
    float *stackGrid, *stackSqGrid, *ptrDataStack, *ptrDataStackSq;
    float *ptrStackGrid, *ptrStackSqGrid;

    const int fullSampNumber = 2 * halfCoher_ + zNum_;
    const int stackGridSize = curxapp * fullSampNumber;
    const float CONVRATIO = M_PI / 180.f;
	
    float* ptrToSembTrace = ptrToSembPanel;

    for (id = 0; id < dipNum_; ++id, ptrToSembTrace += zNum_) {

		curDip = dipStart_ + id * dipStep_;
		curDipInRad = curDip * CONVRATIO;
		tanDipInRad = tan (curDipInRad);
	
		stackGrid   = sf_floatalloc (stackGridSize);
		stackSqGrid = sf_floatalloc (stackGridSize);
		memset (stackGrid,   0, stackGridSize * sizeof (float));
		memset (stackSqGrid, 0, stackGridSize * sizeof (float));
		
		tracesNum = 0;
			
		for (ix = 0; ix < curxapp; ++ix) {		
		    xShift = (ix - leftShift) * xStep_;
		    depthShift = xShift * tanDipInRad;
		    depthShiftSamp = depthShift / zStep_;
		    dataShift = (id + ix * dipNum_) * zNum_;
		
		    ptrDataStack   = ptrToDags_ + dataShift;
		    ptrDataStackSq = ptrToDagsSq_ + dataShift;

		    stackShift = tracesNum * zNum_;

		    ptrStackGrid   = stackGrid + stackShift;
		    ptrStackSqGrid = stackSqGrid + stackShift;

		    zInd = -depthShiftSamp;
				
		    for (iz = 0; iz < zNum_; ++iz, ++ptrStackGrid, ++ptrStackSqGrid, ++zInd) {
				if (zInd < 0 || zInd >= zNum_) continue;
				*ptrStackGrid = *(ptrDataStack + zInd);
				*ptrStackSqGrid = *(ptrDataStackSq + zInd);
		    }		
		    ++tracesNum;
		}

		getSemblanceForTrace (curxapp, stackGrid, stackSqGrid, ptrToSembTrace);  

		free (stackGrid);
		free (stackSqGrid);
    }

    return;
}

int main (int argc, char* argv[]) 
{
    int ix, curxapp, leftShift;
   
/* Initialize RSF */
    sf_init (argc,argv);
/* Input files */
    inDags_   = sf_input("in");
    inDagsSq_ = sf_input("dataSq");

/* check that the input is float */
    if ( SF_FLOAT != sf_gettype (inDags_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */
    if ( SF_FLOAT != sf_gettype (inDagsSq_) ) sf_error ("Need float input: dip-angle gathers in squares");
    /* stacks of amplitude squares in the scattering-angle direction */

/* Output file */
    sembFile_ = sf_output("out");

/* Depth/time axis */
    if ( !sf_histint   (inDags_, "n1", &zNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (inDags_, "d1", &zStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (inDags_, "o1", &zStart_) ) sf_error ("Need o1= in input");
/* Dip angle axis */
    if ( !sf_histint   (inDags_, "n2", &dipNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (inDags_, "d2", &dipStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (inDags_, "o2", &dipStart_) ) sf_error ("Need o2= in input");
/* x axis */
    if ( !sf_histint   (inDags_, "n3", &xNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (inDags_, "d3", &xStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (inDags_, "o3", &xStart_) )   sf_error ("Need o3= in input");
	
    if ( !sf_getint ("xapp",    &xapp_) )    xapp_ = 1;
    /* number of CIGs in the inline-direction processed simultaneously */
    if (!xapp_) {sf_warning ("xapp value is changed to 1"); xapp_ = 1;}

    if ( !sf_getint ("dipapp",    &xdipapp_) ) xdipapp_ = 11;
    /* number of traces in the x-dip direction processed simultaneously */
    if (!xdipapp_) {sf_warning ("dipapp value is changed to 11"); xdipapp_ = 11;}

    if ( !sf_getint ("coher",   &coher_) )   coher_ = 11;
    /* height of a vertical window for semblance calculation */
    if (!coher_) {sf_warning ("coher value is changed to 1"); coher_ = 1;}

    if ( !sf_getint ("scatnum", &scatnum_) ) scatnum_ = 1;
    /* shows how many traces were stacked in the scattering angle direction; 
       if the stack was normalized use the default value 
    */ 

	makeWeight_ = true;
    if ( !sf_getfloat ("s1", &s1_) ) {s1_ = -1.f; makeWeight_ = false; }
    /* minimum semblance value */ 
    if ( !sf_getfloat ("s2", &s2_) ) {s2_ = -1.f; makeWeight_ = false; }
    /* maximum semblance value */ 
	ds_ = s2_ - s1_;

    dagSize_ = zNum_ * dipNum_;
    halfCoher_ = coher_ / 2;    /* yes - this is the integer division */	
    halfXapp_  = xapp_ / 2;     /* this is the integer division too	*/

    for (ix = 0; ix < xNum_; ++ix) {
		
		sf_warning ("CIG %d of %d;", ix + 1, xNum_);
		
		/* xapp for the currect core CIG; it can be changed by the checkBoundary () */
		curxapp = xapp_; 
		/* distance between the core gather and the left side of the aperture */
		leftShift = 0;	

		ptrToSembPanel_ = sf_floatalloc (dagSize_);
		memset (ptrToSembPanel_, 0, dagSize_ * sizeof (float));

		readBlockAroundPoint (ix, halfXapp_, &curxapp, &leftShift);
		buildFilter (curxapp, leftShift, ptrToSembPanel_);

		free (ptrToDags_);
		free (ptrToDagsSq_);

		free (ptrToData_);
		free (ptrToDataSq_);

		sf_floatwrite (ptrToSembPanel_, dagSize_, sembFile_);
		free (ptrToSembPanel_);
    }

    sf_warning (".");

    exit(0);
}
