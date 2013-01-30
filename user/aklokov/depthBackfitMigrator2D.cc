#include "depthBackfitMigrator2D.hh"
#include <rsf.hh>
#include "support.hh"

DepthBackfitMigrator2D::DepthBackfitMigrator2D () {
}

DepthBackfitMigrator2D::~DepthBackfitMigrator2D () {
}

void DepthBackfitMigrator2D::init (int zNum, float zStart, float zStep, 
			 					   int pNum, float pStart, float pStep,
							 	   int xNum, float xStart, float xStep,
							 	   int rNum, float rStart, float rStep,
  								   float* xVol, float* tVol) {
	
	zNum_   = zNum;
	zStep_  = zStep;
	zStart_ = zStart;

	pNum_   = pNum;
	pStep_  = pStep;
	pStart_ = pStart;

	xNum_   = xNum;
	xStep_  = xStep;
	xStart_ = xStart;	

	rNum_   = rNum;
	rStep_  = rStep;
	rStart_ = rStart;	

	xVol_ = xVol;
	tVol_ = tVol;

	return;
}

void DepthBackfitMigrator2D::getImageSample (float* piData, float curX, float curZ, float curP, float* sample) {

	float* xRes = sf_floatalloc (rNum_);
	float* zRes = sf_floatalloc (rNum_);

	memset ( xRes, 0, rNum_ * sizeof (float) );
	memset ( zRes, 0, rNum_ * sizeof (float) );

	iTracer_.traceImage (xVol_, tVol_, curX, curZ, curP, xRes, zRes);

	*sample = 0.f;

	// loop over depth-line
	for (int ir = 0; ir < rNum_; ++ir) {			
		const float lz = zRes[ir];
		if (lz <= 0) continue; // bad point
		const float lx = xRes[ir];

		float curSample (0.f);
		bool goodSample = this->getSample (piData, lx, lz, curP, curSample);

		if (!goodSample) continue; // bad sample
		*sample += curSample;
	}

	free (xRes);
	free (zRes);

	return;
}

void DepthBackfitMigrator2D::processPartialImage (float* piData, float curP, float* piImage) {

	iTracer_.init (zNum_, zStart_, zStep_, 
  			       rNum_, rStart_, rStep_,
			       xNum_, xStart_, xStep_);

	for (int ix = 0; ix < xNum_; ++ix) {
		const float curX = xStart_ + ix * xStep_;
		sf_warning ("%g", curX);
		float* iTrace = piImage + ix * zNum_;
#pragma omp parallel for
		for (int iz = 0; iz < zNum_; ++iz) {
			const float curZ = zStart_ + iz * zStep_;
			this->getImageSample (piData, curX, curZ, curP, iTrace + iz);
		}
	}

	return;
}

bool DepthBackfitMigrator2D::getSample (float* data, const float curX, const float curZ,
										const float curP, float &sample) {

	// limits checking	
	const int izMiddle = (int) ((curZ - zStart_) / zStep_);
	if (izMiddle < 0 || izMiddle >= zNum_) return false;

	const int xSamp = (int) ((curX - xStart_) / xStep_);
	if (xSamp < 0 || xSamp >= xNum_) return false;

	float* trace = data + xSamp * zNum_;

	// middle (main) sample
    
    const float befMiddle = (curZ - zStart_) * 1.f / zStep_ - izMiddle;
    const float aftMiddle = 1.f - befMiddle;

	const float sampleMiddle = aftMiddle * trace[izMiddle] + befMiddle * trace[izMiddle + 1];
    
	sample = sampleMiddle;

	return true;
}
