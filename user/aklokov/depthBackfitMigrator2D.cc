#include "depthBackfitMigrator2D.hh"
#include <rsf.hh>
#include "support.hh"
#include <list>

static bool deleteAll (ImagePoint2D* p) { delete p; return true; }

DepthBackfitMigrator2D::DepthBackfitMigrator2D () {
}

DepthBackfitMigrator2D::~DepthBackfitMigrator2D () {
}

//  Causal integration of a trace[n]
void DepthBackfitMigrator2D::applyCasualIntegration (float *trace, int n) {

    for (int i = 1; i < n; ++i)
        trace[i] += trace[i - 1];
        
    return;      
}

// Anticausal integrations of a trace[n]
void DepthBackfitMigrator2D::applyAnticasualIntegration (float *trace, int n) {

    for (int i = n - 2; i >= 0; --i)
        trace[i] += trace[i + 1];

    return;
}

void DepthBackfitMigrator2D::init (int zNum, float zStart, float zStep, 
			 					   int pNum, float pStart, float pStep,
							 	   int xNum, float xStart, float xStep,
							 	   int rNum, float rStart, float rStep,
								   float dx, float dt,
  								   float* xVol, float* tVol, bool isAA) {
	
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
	
	dx_   = dx;
	dt_   = dt;

	xVol_ = xVol;
	tVol_ = tVol;

	isAA_ = isAA;

	return;
}

bool DepthBackfitMigrator2D::getSample (float* piData, const float& curX, const float& curZ, const float& curP, float &sample) {

	const int xInd = (curX - xStart_) / xStep_;
	
	float befX = curX - (xStart_ + xInd * xStep_);	
	float aftX = xStart_ + (xInd + 1) * xStep_ - curX;

	float bef = befX / xStep_;
	float aft = 1.f - bef;

	float x1 = curX - befX;
	float z1 = curZ + befX * curP;

	float x2 = curX + aftX;
	float z2 = curZ - aftX * curP;

	float sample1, sample2;
	bool goodSample = this->getSampleFromImage (piData, x1, z1, curP, sample1);
	if (!goodSample) return false;
	goodSample = this->getSampleFromImage (piData, x2, z2, curP, sample2);
	if (!goodSample) return false;
	
	sample = bef * sample2 + aft * sample1;

	return true;
}

void DepthBackfitMigrator2D::getImageSample (float* piData, float curX, float curZ, float curP, float* sample) {

	float* xRes = sf_floatalloc (rNum_);
	float* zRes = sf_floatalloc (rNum_);

	memset ( xRes, 0, rNum_ * sizeof (float) );
	memset ( zRes, 0, rNum_ * sizeof (float) );

	iTracer_.traceImage (xVol_, tVol_, curX, curZ, curP, xRes, zRes);

	// filter points
	std::list<ImagePoint2D*> goodPoints;
	for (int ir = 0; ir < rNum_; ++ir) {
		const float lz = zRes[ir];
		if (lz <= 0) continue; // bad point
	    ImagePoint2D* p = new ImagePoint2D (xRes[ir], lz, 0, 0);
		goodPoints.push_back (p);
	}				

	const int listSize = goodPoints.size ();

	std::list<ImagePoint2D*>::iterator iter = goodPoints.begin ();
	*sample = 0.f;
	int count = 0;
	int pcount = 0;
	// loop over depth-line
	for (iter = goodPoints.begin (); iter != goodPoints.end(); ++iter, ++pcount) {
		ImagePoint2D* point = *iter;
		float px = point->x_;
		float pz = point->z_;
		float curP = 0.f;
		if (pcount && pcount != (listSize - 1)) {
			--iter;
			ImagePoint2D* prevPoint = *iter;			
			++iter; ++iter;
			ImagePoint2D* nextPoint = *iter;			
			--iter;
			const float dx = nextPoint->x_ - prevPoint->x_;
			const float dz = nextPoint->z_ - prevPoint->z_;
			curP = -dz / dx;
		}

		float curSample (0.f);
		bool goodSample = this->getSample (piData, px, pz, curP, curSample);			

		if (!goodSample) continue; // bad sample
		*sample += curSample;
		++count;
	}

	if (count)
		*sample /= count;

	// FINISH

	goodPoints.remove_if (deleteAll);

	free (xRes);
	free (zRes);

	return;
}

void DepthBackfitMigrator2D::processData (float* piData) {

	float* ptrData = piData;
	for (int ix = 0; ix < xNum_; ++ix, ptrData += zNum_) {
	    applyCasualIntegration (ptrData, zNum_);
	    applyAnticasualIntegration (ptrData, zNum_);
	}

	return;
}

void DepthBackfitMigrator2D::processPartialImage (float* piData, float curP, float* piImage) {

	iTracer_.init (zNum_, zStart_, zStep_, 
  			       rNum_, rStart_, rStep_,
			       xNum_, xStart_, xStep_,
				   dx_, dt_);

	if (isAA_)
		this->processData (piData); 
/*
	for (int ix = 0; ix < 1; ++ix) {
		const float curX = 6750;
		float* iTrace = piImage + ix * zNum_;
		for (int iz = 0; iz < 1; ++iz) {
			float curZ = 1500;
			this->getImageSample (piData, curX, curZ, curP, iTrace + iz);
		}
	}
*/
	for (int ix = 0; ix < xNum_; ++ix) {
		const float curX = xStart_ + ix * xStep_;
		float* iTrace = piImage + ix * zNum_;
#pragma omp parallel for
		for (int iz = 0; iz < zNum_; ++iz) {
			const float curZ = zStart_ + iz * zStep_;
			this->getImageSample (piData, curX, curZ, curP, iTrace + iz);
		}
	}

	return;
}

bool DepthBackfitMigrator2D::getSampleFromImage (float* data, const float curX, const float curZ,
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

	if (!isAA_) {
		sample = sampleMiddle;
		return true;
	}

	const float p = fabs ( curP );
	const float filterLength = p * xStep_ + zStep_;
    
  	// left sample
 
 	const float zLeft = curZ - filterLength;
 	const int  izLeft = (int) ((zLeft - zStart_) / zStep_); 
	
	if (izLeft < 0) return false;

    const float befLeft = (zLeft - zStart_) * 1.f / zStep_ - izLeft;
    const float aftLeft = 1.f - befLeft;

	const float sampleLeft = aftLeft   * trace[izLeft]   + befLeft   * trace[izLeft   + 1];

	// right sample
 
 	const float zRight = curZ + filterLength;
 	const int  izRight = (int) ((zRight - zStart_) / zStep_); 

	if (izRight >= zNum_ - 1) return false;

    const float befRight = (zRight - zStart_) * 1.f / zStep_ - izRight;
    const float aftRight = 1.f - befRight;

	const float sampleRight = aftRight  * trace[izRight]  + befRight  * trace[izRight  + 1];

	// norm
 
    float imp = zStep_ / (zStep_ + zRight - zLeft);
    imp *= imp;
    
	sample = (2.f * sampleMiddle - sampleLeft - sampleRight) * imp;
		
	return true;
}
