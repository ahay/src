#include "depthBackfitMigrator2D.hh"
#include <rsf.hh>
#include "support.hh"
#include "iTracer2D.hh"

DepthBackfitMigrator2D::DepthBackfitMigrator2D () {
}

DepthBackfitMigrator2D::~DepthBackfitMigrator2D () {
}

void DepthBackfitMigrator2D::init (int zNum, float zStart, float zStep, 
			 					   int pNum, float pStart, float pStep,
							 	   int xNum, float xStart, float xStep) {
	
	zNum_   = zNum;
	zStep_  = zStep;
	zStart_ = zStart;

	pNum_   = pNum;
	pStep_  = pStep;
	pStart_ = pStart;

	xNum_   = xNum;
	xStep_  = xStep;
	xStart_ = xStart;	

	return;
}

void DepthBackfitMigrator2D::processParialImage (float* piData, float curP, float* xVol, float* tVol, float* piImage) {

	ITracer2D iTracer;
	iTracer.init (zNum_, zStart_, zStep_, 
  			      pNum_, pStart_, pStep_,
			      xNum_, xStart_, xStep_);

	float* xRes = sf_floatalloc (pNum_);
	float* zRes = sf_floatalloc (pNum_);

	for (int ix = 0; ix < xNum_; ++ix) {
		const float curX = xStart_ + ix * xStep_;
		for (int iz = 0; iz < zNum_; ++iz) {
			const float curZ = zStart_ + iz * zStep_;
	
			memset ( xRes, 0, pNum_ * sizeof (float) );
			memset ( zRes, 0, pNum_ * sizeof (float) );

			iTracer.traceImage (xVol, tVol, curX, curZ, curP, xRes, zRes);
		}
	}

	// FINISH

	free (xRes);
	free (zRes);

	return;
}
