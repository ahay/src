#include <string.h>
#include "tmigrator2D.h"
#include "support.h"
#include "curveDefinerBase.h"
#include "curveDefinerDipOffset.h"


TimeMigrator2D::TimeMigrator2D () {
}

TimeMigrator2D::~TimeMigrator2D () {
}

void TimeMigrator2D::processGather (Point2D& curGatherCoords, float curOffset, const float* const velTrace, const bool isAzDip,
									float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq) {
   
    const int   tNum     = ip_->zNum;
    const float tStart   = ip_->zStart;
    const float tStep    = ip_->zStep;

    const int   curX     = (int) curGatherCoords.getX ();    
	const float xCIG     = ip_->xStart + curX * ip_->xStep;

    const int   dipNum   = gp_->dipNum;
    const float dipStart = gp_->dipStart;
    const float dipStep  = gp_->dipStep;

	const float curAz    = gp_->sdipStart;

	const float dummy    = 0.f;

    float*     ptrGather = curoffsetGather;

	curveDefiner_->curOffset_ = curOffset;
	curOffset_ = curOffset;

    // compose gather
	for (int id = 0; id < dipNum; ++id) {
		const float curDip = dipStart + id * dipStep;
    
//		int*   ptrMutingMask = mutingMask;					
	    float* ptrImage      = curoffsetImage;
		float* ptrImageSq    = curoffsetImageSq;
		
	    for (int it = 0; it < tNum; ++it, ++ptrGather, ++ptrImage, ++ptrImageSq) { //, ++ptrMutingMask) {
	
//			if (!(*ptrMutingMask)) continue; // the sample is muted

  		    const float curTime = tStart + it * tStep;
			const float migVel = this->getMigVel (velTrace, curTime);
			if (migVel < 1e-3) continue;
	    	
			float sample (0.f);
    		int badRes = this->getSampleByBeam (dummy, xCIG, curTime, curDip, curAz, migVel, isAzDip, sample);
    		if (badRes)
    			sample = this->getSampleByRay (dummy, xCIG, curTime, curDip, curAz, migVel, isAzDip, dummy, dummy);

			*ptrGather += sample;
			*ptrImage += sample;
			*ptrImageSq += sample * sample;
	    }
	}
    
    return;
}
	
int TimeMigrator2D::getSampleByBeam (const float yCIG, const float xCIG, const float curZeroTime, 
	  							     const float curDip, const float curAz, const float migVel, const bool isAzDip, 
									 float& sample) {

	const float xStep   = dp_->xStep;
	const float dipStep = gp_->dipStep;

// the first point
	float time_1 (0.f); float xRayExit_1(0.f);	float yRayExit_1 (0.f);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, -curDip - dipStep * 0.5, curAz, migVel, time_1, xRayExit_1, yRayExit_1);

// the second point
	float time_2 (0.f); float xRayExit_2(0.f);	float yRayExit_2 (0.f);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, -curDip + dipStep * 0.5, curAz, migVel, time_2, xRayExit_2, yRayExit_2);

	int curXSamp = xRayExit_1 / xStep;
	if (curXSamp * xStep < xRayExit_1) ++curXSamp;
	
	int count (0);
	float curXpos = curXSamp * xStep;
	while (curXpos < xRayExit_2 + 1e-6) {
		if (curXpos < dataXMin_ || curXpos > dataXMax_) {curXpos += xStep; continue; }
		float dx = curXpos - xCIG;
		float p (0.f);
		float curTime = curveDefiner_->getEscapeTime (yCIG, dx, curZeroTime, migVel, p);
		sample += this->getSampleFromData (yCIG, curXpos, curTime, p);
		curXpos += xStep;
		++count;
	}
	if (!count) return -1;

	sample /= count;

	return 0;
}

float TimeMigrator2D::getSampleByRay (const float yCIG, const float xCIG, const float curZeroTime, 
									   const float curDip, const float curAz, const float migVel, const bool isAzDip,
									   const float yEmAngle, const float xEmAngle) {

	const float xStep   = dp_->xStep;

	float time (0.f); float xRayExit (0.f);	float yRayExit (0.f);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, -curDip, curAz, migVel, time, xRayExit, yRayExit);

	// check if the ray exits outside the acquisition surface
	if (xRayExit - dataXMin_ < -1e-4 || xRayExit - dataXMax_ > 1e-4) return 0.f;
	if (time < dataZMin_ || time > dataZMax_) return 0.f;

	int curXSamp = xRayExit / xStep;
	float curXpos = curXSamp * xStep;
	if (curXpos < dataXMin_ || curXpos > dataXMax_) {return 0.f;}
	float dx = curXpos - xCIG;

	float p (0.f);	

	float curTime1 = curveDefiner_->getEscapeTime (yCIG, dx, curZeroTime, migVel, p);
	float sample1  = this->getSampleFromData (yCIG, curXpos, curTime1, p);	

	curXpos += xStep;
	if (curXpos < dataXMin_ || curXpos > dataXMax_) {return 0.f;}
	dx = curXpos - xCIG;

	float curTime2 = curveDefiner_->getEscapeTime (yCIG, dx, curZeroTime, migVel, p);
	float sample2  = this->getSampleFromData (yCIG, curXpos, curTime2, p);		
	
	float bef = (xRayExit - curXSamp * xStep ) / xStep;
	float aft = 1.f - bef;
	
	float sample = bef * sample2 + aft * sample1;
	
	return sample;
}

float TimeMigrator2D::getSampleFromData (const float geoY, const float geoX1, const float ti, const float p) {
	
	int zNum_ = dp_->zNum;

	float zStep_ = dp_->zStep;
	float xStep_ = dp_->xStep;

	float zStart_ = dp_->zStart;
	float xStart_ = dp_->xStart;

	float geoX = isCMP_ ? geoX1 - curOffset_ * 0.5 : geoX1;

	const int itMiddle = (ti - zStart_) / zStep_;
	if (itMiddle < 0 || itMiddle >= zNum_) return 0.f;

	const int xSamp = (geoX - xStart_) / xStep_;
	if (xSamp < 0 || xSamp >= zNum_) return 0.f;

	float* const trace = ptrToData_ + xSamp * zNum_;

	// middle (main) sample
    
    const float befMiddle = (ti - zStart_) * 1.f / zStep_ - itMiddle;
    const float aftMiddle = 1.f - befMiddle;

	const float sampleMiddle = aftMiddle * trace[itMiddle] + befMiddle * trace[itMiddle + 1];
    
	if (!isAA_) 
		return sampleMiddle;

	const float filterLength = p * xStep_ + zStep_;
    
  	// left sample
 
 	const float timeLeft = ti - filterLength;
 	const int     itLeft = (timeLeft - zStart_) / zStep_; 
	
	if (itLeft < 0) return 0.f;

    const float befLeft = (timeLeft - zStart_) * 1.f / zStep_ - itLeft;
    const float aftLeft = 1.f - befLeft;

	const float sampleLeft = aftLeft   * trace[itLeft]   + befLeft   * trace[itLeft   + 1];

	// right sample
 
 	const float timeRight = ti + filterLength;
 	const int     itRight = (timeRight - zStart_) / zStep_; 

	if (itRight >= zNum_ - 1) return 0.f;

    const float befRight = (timeRight - zStart_) * 1.f / zStep_ - itRight;
    const float aftRight = 1.f - befRight;

	const float sampleRight = aftRight  * trace[itRight]  + befRight  * trace[itRight  + 1];

	// norm
 
    float imp = zStep_ / (zStep_ + timeRight - timeLeft);
    imp *= imp;
    
	const float aaSample = (2.f * sampleMiddle - sampleLeft - sampleRight) * imp;
		
	return aaSample;
}
