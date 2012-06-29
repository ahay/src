#include <rsf.hh>
#include <string.h>
#include "dmigrator2D.hh"
#include "support.hh"
#include "curveDefinerBase.hh"
#include "curveDefinerDipOffset.hh"

#include <fstream>

#include <iostream>
using namespace std;

DepthMigrator2D::DepthMigrator2D () {
}

DepthMigrator2D::~DepthMigrator2D () {
}

void DepthMigrator2D::processGather (Point2D& curGatherCoords, const float* const data, float* dag, float* aCig) {

	// CONSTANTS

	// inline position
	const float xCIG = curGatherCoords.getX ();
	// depth	
    const int   zNum   = gp_->zNum;
	const float zStart = gp_->zStart;
	const float zStep  = gp_->zStep;
	// dip-angle	
	const int   dipNum   = gp_->dipNum;
	const float dipStart = gp_->dipStart;
	const float dipStep  = gp_->dipStep;
	// scattering-angle
	const int   scatNum   = gp_->scatNum;
	const float scatStart = gp_->scatStart;
	const float scatStep  = gp_->scatStep;
	// size of dip-angle gather
	const int dagSize = zNum * dipNum;
	// travel-times rays number
	const int ttNum = wavefrontTracer_.wp_.rNum;

	// ACTION

	// dip-angle gather for current scattering angle
	float* curDag = new float [dagSize];
	memset ( curDag, 0, dagSize * sizeof (float) );

	// loop over depth samples
    for (int iz = 0; iz < zNum; ++iz) {	
		const float curZ = zStart + iz * zStep;		
		travelTimes_ = new EscapePoint [ttNum];		
		this->calcTravelTimes (curZ, xCIG, travelTimes_);
		// loop over scattering-angle
		for (int is = 0; is < scatNum; ++is) {
			const float curScatAngle = scatStart + is * scatStep;
			// loop over dip-angle
			for (int id = 0; id < dipNum; ++id)	{
				const float curDipAngle = -1 * (dipStart + id * dipStep); // "-1" is to consist with an agreement
				
				float sample (0.f);
				int isBad = this->getSampleByBeam (curScatAngle, curDipAngle, sample);	 
				if (isBad)
					sample = 0.f; // some function may be here

				const int ind = id * zNum + iz;
				curDag [ind] += sample;
				aCig [is*zNum + iz] += sample;
			}
			// add current scattering-angle migration to the main result
			float* pTo = dag;
			float* pFrom = curDag;
			for (int id = 0; id < dagSize; ++id, ++pTo, ++pFrom)
				*pTo += *pFrom;
		}

		delete [] travelTimes_;

		sf_warning ("iz %d\n", iz);
	}

	delete [] curDag;

	return;
}

int DepthMigrator2D::getSampleByBeam (float curScatAngle, float curDipAngle, float& sample) {

	const int   hNum   = dp_->hNum;
	const float hStart = dp_->hStart; 
	const float hStep  = dp_->hStep; 
	
	const int   xNum   = dp_->xNum;
	const float xStart = dp_->xStart; 
	const float xStep  = dp_->xStep; 

	const float dipStep = gp_->dipStep;

	const int ttNum = wavefrontTracer_.wp_.rNum;

	float baseDir = curDipAngle + 0.5 * curScatAngle;
	const float shiftDir = 0.5 * dipStep;

	// calc recs lane
	float dir1 = baseDir - shiftDir;
	EscapePoint* escPointRec1 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttNum, dir1, *escPointRec1);

	float dir2 = baseDir + shiftDir;
	EscapePoint* escPointRec2 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttNum, dir2, *escPointRec2);
	
	float recLaneLeft  = escPointRec1->x;
	float recLaneRight = escPointRec2->x;
	if (recLaneLeft > recLaneRight) { float temp = recLaneRight; recLaneRight = recLaneLeft; recLaneLeft = temp; }
	delete escPointRec1;
	delete escPointRec2;

	// calc sources lane
	baseDir = curDipAngle - 0.5 * curScatAngle;
	float dir3 = baseDir - shiftDir;
	EscapePoint* escPointRec3 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttNum, dir3, *escPointRec3);

	float dir4 = baseDir + shiftDir;
	EscapePoint* escPointRec4 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttNum, dir4, *escPointRec4);			

	float srcLaneLeft  = escPointRec3->x;
	float srcLaneRight = escPointRec4->x;
	if (srcLaneLeft > srcLaneRight) { float temp = srcLaneRight; srcLaneRight = srcLaneLeft; srcLaneLeft = temp; }
	delete escPointRec3;
	delete escPointRec4;

	// loop over receivers
	int count (0);

	float curRecPos = xStart + xStep * (xNum - 1);
	while (curRecPos > recLaneRight) curRecPos -= xStep;

	while (curRecPos > recLaneLeft) {
		// loop over offsets
		for (int ih = 0; ih < hNum; ++ih) {
			const float curOffset = hStart + ih * hStep;
			const float curSrcPos = curRecPos - curOffset;
			if (curSrcPos < srcLaneLeft || curSrcPos > srcLaneRight)
				continue;
			// get escape point
			float timeToRec (0.f);
			float recAbsP (0.f);
			bool onSurf;
			this->getRayToPoint (curRecPos, dir1, dir2, timeToRec, recAbsP, onSurf);
			if (!onSurf) continue; // the ray does not reach the surface
			float timeToSrc (0.f);
			float srcAbsP (0.f);
			this->getRayToPoint (curSrcPos, dir3, dir4, timeToSrc, srcAbsP, onSurf);
			if (!onSurf) continue; // the ray does not reach the surface
			float curTime = (timeToRec + timeToSrc) * 1000; // transform time to "ms"
			sample += this->getSampleFromData (curOffset, 0, curRecPos, curTime, recAbsP);
			++count;
		}
		curRecPos -= xStep;
	}
			
	if (!count) return -1; // no sample is returned
	sample /= count;

	return 0;
}

void DepthMigrator2D::getRayToPoint (float curRecPos, float dir1, float dir2, float& timeToPoint, float& pointAbsP, bool& isSurf) {

	//

	int size = wavefrontTracer_.wp_.rNum;
	float basedir = dir1 > dir2 ? dir1 : dir2;

	//

	EscapePoint* pTT = travelTimes_;
	
	float ttDir = pTT->startDir;

	int count (0);
	
	while (ttDir > basedir && count < size) {
		++pTT; 
		ttDir = pTT->startDir; 
		++count;
	}
	if (count == size) return;
	--pTT;

	float ttX = pTT->x;
	while (ttX > curRecPos && count < size) { ++pTT; ttX = pTT->x; ++count; } 	
	if (count == size) return;

	// get two basic points around the target point
	EscapePoint* rightPoint = pTT - 1;
	EscapePoint* leftPoint  = pTT;				

	if (!leftPoint->isSurf && !rightPoint->isSurf) {
		isSurf = false; return;
	} else if (!leftPoint->isSurf) {
		timeToPoint = rightPoint->t;		
		pointAbsP   = rightPoint->p;	
		isSurf = true;	
		return;
	} else if (!rightPoint->isSurf) {
		timeToPoint = leftPoint->t;		
		pointAbsP   = leftPoint->p;		
		isSurf = true;	
		return;
	} else {
		float bef = ( rightPoint->x - curRecPos ) / (rightPoint->x - leftPoint->x);
		float aft = 1 - bef;

		timeToPoint = leftPoint->t * bef + rightPoint->t * aft;
		pointAbsP   = leftPoint->p * bef + rightPoint->p * aft;
		isSurf = true;	
		return;
	}	
}

void DepthMigrator2D::getEscPointByDirection (EscapePoint* escPoints, int size, float targetStartDir, EscapePoint &resEscPoint) {


	if (targetStartDir < startDirMin_) {resEscPoint.isSurf = false; return;}
	if (targetStartDir > startDirMax_) {resEscPoint.isSurf = false; return;}
		

	EscapePoint* pEscPoint = escPoints;
	float curStartDir = pEscPoint->startDir;
	if (curStartDir < targetStartDir) {
		resEscPoint = *pEscPoint;
		return;
	}

	int count (0);

	while (curStartDir > targetStartDir && count <= size) {
		++pEscPoint;
		curStartDir = pEscPoint->startDir;
		++count;
	}

	if (count == size) {
		resEscPoint = *(pEscPoint - 1);
		return;
	}

//	!!! this interolation is NOT the best solution

	float bef = ( pEscPoint->startDir - targetStartDir ) / ( pEscPoint->startDir - (pEscPoint - 1)->startDir );
	float aft = 1.f - bef;

	resEscPoint.x = pEscPoint->x * aft + (pEscPoint - 1)->x * bef;
	resEscPoint.p = pEscPoint->p * aft + (pEscPoint - 1)->p * bef;
	resEscPoint.t = pEscPoint->t * aft + (pEscPoint - 1)->t * bef;

	return;
}

void DepthMigrator2D::calcTravelTimes (float curZ, float curX, EscapePoint* escPoints) {

	const float xCIG = ip_->xStart + curX * ip_->xStep;
	wavefrontTracer_.getEscapePoints (xCIG, curZ, escPoints);
		
	startDirMin_ = wavefrontTracer_.wp_.rStart - 180.f;
	startDirMax_ = wavefrontTracer_.wp_.rStart +  (wavefrontTracer_.wp_.rNum - 1) * wavefrontTracer_.wp_.rStep - 180.f;


	return;
}

void DepthMigrator2D::processGatherOLD (Point2D& curGatherCoords, float curOffset, const float* const velTrace, const bool isAzDip,
									float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq) {
   
    const int   zNum     = ip_->zNum;
    const float zStart   = ip_->zStart;
    const float zStep    = ip_->zStep;

    const int   curX     = (int) curGatherCoords.getX ();    
	const float xCIG     = ip_->xStart + curX * ip_->xStep;

    const int   dipNum   = gp_->dipNum;
    const float dipStart = gp_->dipStart;
    const float dipStep  = gp_->dipStep;

	const float curAz    = gp_->sdipStart;




//wavefrontTracer_.setParams (raysNum, raysStep, raysStart);

//	const float dummy    = 0.f;

//    float*     ptrGather = curoffsetGather;

//	curveDefiner_->curOffset_ = curOffset;
//	curOffset_ = (int) curOffset;

	// calculate traveltimes

//	escPoints = new EscapePoint [raysNum_ * zNum];
//	EscapePoint* pEscPoints = escPoints;
		
//    for (int iz = 0; iz < zNum; ++iz, pEscPoints += raysNum) {
//        const float curZ = ip_->zStart + iz * ip_->zStep;	
//		this->getEscapePoints (xCIG, curZ, pEscPoints);
//    }

	
//			if (!(*ptrMutingMask)) continue; // the sample is muted

//  		    const float curTime = zStart + iz * zStep;
//			const float migVel = this->getMigVel (velTrace, curTime);
//			if (migVel < 1e-3) continue;
	    	
//			float sample (0.f);
//    		int badRes = this->getSampleByBeam (dummy, xCIG, curTime, curDip, curAz, migVel, isAzDip, sample);
//    		if (badRes)
//    			sample = this->getSampleByRay (dummy, xCIG, curTime, curDip, curAz, migVel, isAzDip, dummy, dummy);

//			*ptrGather += sample;
//			*ptrImage += sample;
//			*ptrImageSq += sample * sample;
//	    }
//	}

   
    return;
}

void DepthMigrator2D::getSampleByRay (EscapePoint& escPoint, float& sample) {

	if (!escPoint.isSurf) {
		sample = 0.f;
		return; 
	}

	float xStartData_ = dp_->xStart;
	float xStepData_ = dp_->xStep;

	float dataXMin_ = dp_->xStart;
	float dataXMax_ = dataXMax_ + (dp_->xNum - 1) * dp_->xStep;
	float dataTMin_ = dp_->zStart;
	float dataTMax_ = dataTMin_ + (dp_->zNum - 1) * dp_->zStep;


	if (escPoint.x - dataXMin_ < -1e-4 || escPoint.x - dataXMax_ > 1e-4) return;
	if (escPoint.t < dataTMin_ || escPoint.t > dataTMax_) return;

	int curXSamp = (escPoint.x - xStartData_) / xStepData_;
	float posLeftTrace = curXSamp * xStepData_ + xStartData_;
	float posRightTrace = (curXSamp + 1) * xStepData_ + xStartData_;
	
	float p = escPoint.p;
	float timeLeft  = escPoint.t - (escPoint.x - posLeftTrace) * p * 0.001;
	float timeRight = escPoint.t + (posRightTrace - escPoint.x) * p * 0.001;	

	float absP = fabs (p);
	float sampleLeft  = this->getSampleFromData (escPoint.offset, 0, posLeftTrace,  1000 * timeLeft, absP);
	float sampleRight = this->getSampleFromData (escPoint.offset, 0, posRightTrace, 1000 * timeRight, absP);

	float bef = (escPoint.x - posLeftTrace) / xStepData_;
	float aft = 1.f - bef;
	
	sample = bef * sampleRight + aft * sampleLeft;
	
	return;
}

float DepthMigrator2D::getSampleFromData (const float h, const float geoY, const float geoX1, const float ti, const float p) {
	
	int tNum = dp_->zNum;
	int xNum = dp_->xNum;

	float tStep = dp_->zStep;
	float xStep = dp_->xStep;

	float tStart = dp_->zStart;
	float xStart = dp_->xStart;

	float geoX = isCMP_ ? geoX1 - curOffset_ * 0.5 : geoX1;

	const int itMiddle = (int) ((ti - tStart) / tStep);
	if (itMiddle < 0 || itMiddle >= tNum) return 0.f;

	const int xSamp = (int) ((geoX - xStart) / xStep);
	if (xSamp < 0 || xSamp >= xNum) return 0.f;

	// offset

	const int offsetInd = h / dp_->hStep;
//	sf_warning ("offset %g %d", h, offsetInd);
	if (offsetInd >= dp_->hNum || offsetInd < 0) return 0.f;

	float* const trace = ptrToData_ + xSamp * tNum + xNum * tNum * offsetInd;

	// middle (main) sample
    
    const float befMiddle = (ti - tStart) * 1.f / tStep - itMiddle;
    const float aftMiddle = 1.f - befMiddle;

	const float sampleMiddle = aftMiddle * trace[itMiddle] + befMiddle * trace[itMiddle + 1];
    
	if (!isAA_) 
		return sampleMiddle;

	const float filterLength = p * xStep + tStep;
    
  	// left sample
 
 	const float timeLeft = ti - filterLength;
 	const int     itLeft = (int) ((timeLeft - tStart) / tStep); 
	
	if (itLeft < 0) return 0.f;

    const float befLeft = (timeLeft - tStart) * 1.f / tStep - itLeft;
    const float aftLeft = 1.f - befLeft;

	const float sampleLeft = aftLeft   * trace[itLeft]   + befLeft   * trace[itLeft   + 1];

	// right sample
 
 	const float timeRight = ti + filterLength;
 	const int     itRight = (int) ((timeRight - tStart) / tStep); 

	if (itRight >= tNum - 1) return 0.f;

    const float befRight = (timeRight - tStart) * 1.f / tStep - itRight;
    const float aftRight = 1.f - befRight;

	const float sampleRight = aftRight  * trace[itRight]  + befRight  * trace[itRight  + 1];

	// norm
 
    float imp = tStep / (tStep + timeRight - timeLeft);
    imp *= imp;
    
	const float aaSample = (2.f * sampleMiddle - sampleLeft - sampleRight) * imp;
		
	return aaSample;
}
