#include "tmigrator3D.hh"
#include "support.hh"
#include <string.h>
#include <algorithm>
#include <math.h>
#include "curveDefinerBase.hh"
#include "curveDefinerDipOffset.hh"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

TimeMigrator3D::TimeMigrator3D () {
}

TimeMigrator3D::~TimeMigrator3D () {
}

void TimeMigrator3D::processGather (Point2D& curGatherCoords, float curOffset, const float* const velTrace,
      							    const bool isAzDip, float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq) {


    const int   tNum     = ip_->zNum;
    const float tStart   = ip_->zStart;
    const float tStep    = ip_->zStep;

	const int   curX     = (int) curGatherCoords.getX ();    
	const int   curY     = (int) curGatherCoords.getY ();    

	const float xCIG     = ip_->xStart + curX * ip_->xStep;
	const float yCIG     = ip_->yStart + curY * ip_->yStep;

    const int   dipNum   = gp_->dipNum;
    const float dipStart = gp_->dipStart;
    const float dipStep  = gp_->dipStep;

    const int   sdipNum   = gp_->sdipNum;
    const float sdipStart = gp_->sdipStart;
    const float sdipStep  = gp_->sdipStep;

	curveDefiner_->curOffset_ = curOffset;

 //   const int gatherSize = tNum * dipNum * dipNum;

    float*     ptrGather = curoffsetGather;
	float*     pTaper    = stackTaper_;

    // compose gather
    for (int idy = 0; idy < sdipNum; ++idy) {
        const float curSDip = sdipStart + idy * sdipStep;
	    for (int idx = 0; idx < dipNum; ++idx, ++pTaper) {
    	    const float curDip = dipStart + idx * dipStep;

			float* ptrImage  = curoffsetImage;
			float* ptrImageSq = curoffsetImageSq;

			const float dummy = 0.f;	
			const float w = *pTaper; // stacking weight

#ifdef _OPENMP 
#pragma omp parallel for
#endif
		    for (int it = 0; it < tNum; ++it) {
				const float curTime = tStart + it * tStep;
				const float migVel = this->getMigVel (velTrace, curTime);
				if (migVel < 1e-3) continue;

				float sample (0.f);
	    		int badRes = this->getSampleByBeam (yCIG, xCIG, curTime, curSDip, curDip, migVel, isAzDip, sample);
	    		if (badRes)
	    			sample = this->getSampleByRay (yCIG, xCIG, curTime, curSDip, curDip, migVel, isAzDip, dummy, dummy);
	
				const int gInd = it + (idy * dipNum + idx) * tNum;
				sample *= w;

				ptrGather[gInd] += sample;

				ptrImage[it] += sample;
				ptrImageSq[it] += sample * sample;
		    }
		}
	}	

	return;
}

float TimeMigrator3D::getSampleFromData (const float geoY, const float geoX, const float ti, const float p) {
	
	int zNum_ = dp_->zNum;
	int xNum_ = dp_->xNum;
	int yNum_ = dp_->yNum;

	float zStep_ = dp_->zStep;
	float xStep_ = dp_->xStep;
	float yStep_ = dp_->yStep;

	float zStart_ = dp_->zStart;
	float xStart_ = dp_->xStart;
	float yStart_ = dp_->yStart;


	const int itMiddle = (int) ((ti - zStart_) / zStep_);
	if (itMiddle < 0 || itMiddle >= zNum_) return 0.f;
	const int xSamp = (int) ((geoX - xStart_) / xStep_);
	if (xSamp < 0 || xSamp >= xNum_) return 0.f;
	const int ySamp = (int) ((geoY - yStart_) / yStep_);
	if (ySamp < 0 || ySamp >= yNum_) return 0.f;


	float* const trace = ptrToData_ + (ySamp * xNum_ + xSamp) * zNum_;

	// middle (main) sample
    
    const float befMiddle = (ti - zStart_) * 1.f / zStep_ - itMiddle;
    const float aftMiddle = 1.f - befMiddle;

	const float sampleMiddle = aftMiddle * trace[itMiddle] + befMiddle * trace[itMiddle + 1];
    
//	if (!rp_.isAA) 
//		return sampleMiddle;

	const float filterLength = p * xStep_ + zStep_;
    
  	// left sample
 
 	const float timeLeft = ti - filterLength;
 	const int     itLeft = (int) ((timeLeft - zStart_) / zStep_); 
	
	if (itLeft < 0) return 0.f;

    const float befLeft = (timeLeft - zStart_) * 1.f / zStep_ - itLeft;
    const float aftLeft = 1.f - befLeft;

	const float sampleLeft = aftLeft   * trace[itLeft]   + befLeft   * trace[itLeft   + 1];

	// right sample
 
 	const float timeRight = ti + filterLength;
 	const int     itRight = (int) ((timeRight - zStart_) / zStep_); 

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

	
int TimeMigrator3D::getSampleByBeam (const float yCIG, const float xCIG, const float curZeroTime, 
									 const float curDipY, const float curDipX, const float migVel, const bool isAzDip,
									 float& sample) {

	const float xStart   = dp_->xStart;
	const float yStart   = dp_->yStart;

	const float xStep   = dp_->xStep;
	const float yStep   = dp_->yStep;

	const float dipStep  = gp_->dipStep;
	const float sdipStep = isAzDip ? 1.f : gp_->sdipStep; // for dip/azimuth mode 
	

	// escape points

	float xRayExit [4];
	float yRayExit [4];
	float timeExit [4];

	float curDip (0.f);
	float curAz (0.f);
													 // if it is not correct to use half of azimuth step; temporarily it is set to dip step 
	if (isAzDip) { curDip = curDipX - dipStep * 0.5; curAz  = curDipY - dipStep * 0.5; }
	else this->getAzDipFromXY (curDipY - sdipStep * 0.5, curDipX - dipStep * 0.5, curDip, curAz);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, curDip, curAz, migVel, timeExit [0], xRayExit [0], yRayExit [0]);

	if (isAzDip) { curDip = curDipX - dipStep * 0.5; curAz  = curDipY + dipStep * 0.5; }
	else this->getAzDipFromXY (curDipY + sdipStep * 0.5, curDipX - dipStep * 0.5, curDip, curAz);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, curDip, curAz, migVel, timeExit [1], xRayExit [1], yRayExit [1]);

	if (isAzDip) { curDip = curDipX + dipStep * 0.5; curAz  = curDipY + dipStep * 0.5; }
	else this->getAzDipFromXY (curDipY + sdipStep * 0.5, curDipX + dipStep * 0.5, curDip, curAz);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, curDip, curAz, migVel, timeExit [2], xRayExit [2], yRayExit [2]);

	if (isAzDip) { curDip = curDipX + dipStep * 0.5; curAz  = curDipY - dipStep * 0.5; }
	else this->getAzDipFromXY (curDipY - sdipStep * 0.5, curDipX + dipStep * 0.5, curDip, curAz);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, curDip, curAz, migVel, timeExit [3], xRayExit [3], yRayExit [3]);

	Point2D poly [4];
	for (int ic = 0; ic < 4; ++ic) {
		poly[ic].setX (xRayExit [ic]);
		poly[ic].setY (yRayExit [ic]);
	}

	float minX = *min_element (xRayExit, xRayExit + 3);
	float maxX = *max_element (xRayExit, xRayExit + 3);
	float minY = *min_element (yRayExit, yRayExit + 3);
	float maxY = *max_element (yRayExit, yRayExit + 3);

	int curXSamp = (int) ((minX - xStart) / xStep);
	if (curXSamp * xStep - xStart < minX) ++curXSamp;

	int curYSamp = (int) ((minY - yStart) / yStep);
	if (curYSamp * yStep - yStart < minY) ++curYSamp;

	float curYpos = curYSamp * yStep + yStart;
	float curXposStart = curXSamp * xStep + xStart;

	int count (0);

	while (curYpos < maxY + 1e-6) {
		if (curYpos < dataYMin_ || curYpos > dataYMax_) {curYpos += yStep; continue; }
		float dy = curYpos - yCIG;
		float curXpos = curXposStart;
		bool isInside = false;	

		while (curXpos < maxX + 1e-6) {
			if (curXpos < dataXMin_ || curXpos > dataXMax_) {curXpos += xStep; continue; }		
		
			Point2D p0 (curXpos, curYpos);

    		if ( !this->isPointInsidePoly (poly, 4, p0) ) {
				if (isInside) curXpos = 2 * maxX;					
				else curXpos += xStep;			
				continue; 
			}		

			isInside = true;

			float dx = curXpos - xCIG;
			float p (0.f);
			float curTime = curveDefiner_->getEscapeTime (dy, dx, curZeroTime, migVel, p);
			sample += this->getSampleFromData (curYpos, curXpos, curTime, p);
			curXpos += xStep;
			++count;
		}	
		curYpos += yStep;
	}
	if (!count) return -1;

	sample /= count;

	return 0;
}

float TimeMigrator3D::getSampleByRay (const float yCIG, const float xCIG, const float curZeroTime, 
									  const float curDipY, const float curDipX, const float migVel, const bool isAzDip,
									  const float yEmAngle, const float xEmAngle) {

//	int zNum_ = dp_->zNum;
	int xNum_ = dp_->xNum;
	int yNum_ = dp_->yNum;

//	float zStep_ = dp_->zStep;
	float xStep_ = dp_->xStep;
	float yStep_ = dp_->yStep;

//	float zStart_ = dp_->zStart;
	float xStart_ = dp_->xStart;
	float yStart_ = dp_->yStart;

	float curDip (0.f); float curAz (0.f);
	if (isAzDip) { curDip = curDipX; curAz  = curDipY; }
	else this->getAzDipFromXY (curDipY, curDipX, curDip, curAz);
	float timeExit (0.f); float xRayExit (0.f);	float yRayExit (0.f);
	curveDefiner_->getEscapePoint (yCIG, xCIG, curZeroTime, curDip, curAz, migVel, timeExit, xRayExit, yRayExit);

	const int xSamp = (int) ((xRayExit - xStart_) / xStep_);
	if (xSamp < 1 || xSamp >= xNum_ - 1) return 0.f;
	const int ySamp = (int) ((yRayExit - yStart_) / yStep_);
	if (ySamp < 1 || ySamp >= yNum_ - 1) return 0.f;

	// get p
//	float p = 0;
//	float curTime = curveDefiner_->getEscapeTime (yRayExit - yCIG, xRayExit - xCIG, curZeroTime, migVel, p);

	// which triangle contains the escape point

	Point2D p0 (xRayExit, yRayExit);

	Point2D poly [3];	
	poly[0].setX (xStart_ + xSamp * xStep_);
	poly[0].setY (yStart_ + ySamp * yStep_);

	poly[1].setX (xStart_ + (xSamp + 1) * xStep_);
	poly[1].setY (yStart_ + ySamp * yStep_);

	poly[2].setX (xStart_ + xSamp * xStep_);
	poly[2].setY (yStart_ + (ySamp + 1) * yStep_);

	if (!this->isPointInsidePoly (poly, 3, p0) ) {
		poly[0].setX (xStart_ + (xSamp + 1) * xStep_);
		poly[0].setY (yStart_ + (ySamp + 1) * yStep_);
	}

	float samples [3];
//	float weights [3];

//	float dummy = 0.f;

	for (int ic = 0; ic < 3; ++ic) {
		float curP = 0.f;
		float curTime = curveDefiner_->getEscapeTime (poly[ic].getY () - yCIG, poly[ic].getX () - xCIG, curZeroTime, migVel, curP);
		samples[ic] = this->getSampleFromData (poly[ic].getY (), poly[ic].getX (), curTime, curP);	
	}

	const float x13 = poly[0].getX () - poly[2].getX ();
	const float x32 = poly[2].getX () - poly[1].getX ();

	const float y23 = poly[1].getY () - poly[2].getY ();
	const float y13 = poly[0].getY () - poly[2].getY ();

	const float denom = y23 * x13 + x32 * y13; 

	float w1 = ( y23 * (xRayExit - poly[2].getX())  + x32 * (yRayExit - poly[2].getY()) ) / denom;
	float w2 = ( -y13 * (xRayExit - poly[2].getX()) + x13 * (yRayExit - poly[2].getY()) ) / denom;
	float w3 = 1 - w1 - w2;

	float sample = w1 * samples[0] + w2 * samples[1] + w3 * samples[2];

	return sample;
}

void TimeMigrator3D::getStackTaper (const float edgeTaper, const bool isDipAz) {

    const int   dipNum   = gp_->dipNum;
    const float dipStart = gp_->dipStart;
    const float dipStep  = gp_->dipStep;

    const int   sdipNum   = gp_->sdipNum;
    const float sdipStart = gp_->sdipStart;
    const float sdipStep  = gp_->sdipStep;

	const float dipMax = dipStep * dipNum / 2;
	const float edgeDip = dipMax - edgeTaper;

	const int taperSize = sdipNum * dipNum;
	stackTaper_ = new float [taperSize];

	float* pTaper = stackTaper_;
	
	if (isDipAz) {
	    for (int it = 0; it < taperSize; ++it, ++pTaper)
			*pTaper = 1.f;
		return;
	}

    for (int idy = 0; idy < sdipNum; ++idy) {
        const float curSDip = sdipStart + idy * sdipStep;
		const float sdip2 = curSDip * curSDip;
	    for (int idx = 0; idx < dipNum; ++idx, ++pTaper) {
    	    const float curDip = dipStart + idx * dipStep;
			const float dip2 = curDip * curDip;

			// stacking taper	
			const float dip = sqrt (sdip2 + dip2);
			float w = 1.f;		
			if (dip > edgeDip) {
				if (dip > dipMax) w = 0.f;
				else w = 1.f - (dip - edgeDip) / edgeTaper;
			}	
			*pTaper = w;
		}
	}

	return;
}
