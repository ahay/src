#include <rsf.hh>
#include <string.h>
#include "dmigrator2D.hh"
#include "support.hh"

DepthMigrator2D::DepthMigrator2D () {
}

DepthMigrator2D::~DepthMigrator2D () {
}

void DepthMigrator2D::processGather (Point2D& curGatherCoords, const float* const data, float* image, float* dag, float* aCig) {

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
	// size of internal scattering-angle gather
	const int scatSize = zNum * scatNum;
	// velocity model limits
	const float velModelDepthMin = vp_->zStart;
	const float velModelDepthMax = velModelDepthMin + (vp_->zNum - 1) * vp_->zStep;

	// ACTION

	// dip-angle gather for current scattering angle
	double* curDag = new double [dagSize];
	memset ( curDag, 0, dagSize * sizeof (double) );
	// internal scattering angle gather
	double* curCig = new double [scatSize];
	memset ( curCig, 0, scatSize * sizeof (double) );
	// internal image
	double* curImage = new double [zNum];
	memset ( curImage, 0, zNum * sizeof (double) );

	// mask for dip-angle gather for current scattering angle
	int* maskDag = new int [dagSize];
	memset ( maskDag, 0, dagSize * sizeof (int) );
	// mask for internal scattering angle gather
	int* maskCig = new int [scatSize];
	memset ( maskCig, 0, scatSize * sizeof (int) );
	// mask for internal image
	int* maskImage = new int [zNum];
	memset ( maskImage, 0, zNum * sizeof (int) );

	// loop over depth samples
    for (int iz = 0; iz < zNum; ++iz) {	
		const float curZ = zStart + iz * zStep;		
		if (curZ < velModelDepthMin || curZ > velModelDepthMax)
			continue;
		travelTimes_ = new EscapePoint [ttRayNum_];		
		this->calcTravelTimes (curZ, xCIG, travelTimes_);
		// loop over scattering-angle
		for (int is = 0; is < scatNum; ++is) {
			const float curScatAngle = scatStart + is * scatStep;
			// loop over dip-angle
			for (int id = 0; id < dipNum; ++id)	{
				const float curDipAngle = -1 * (dipStart + id * dipStep); // "-1" is to consist with an agreement
				
				float sample (0.f);
				bool isGood = this->getSampleByBeam (curScatAngle, curDipAngle, sample);	 
				if (!isGood && curScatAngle < 1e-6) // implemented for zero-offset only
					isGood = this->getSampleByRay (curDipAngle, sample);

				if (!isGood)
					continue;

				const int indDag = id * zNum + iz;
				curDag [indDag]  += sample;
				maskDag [indDag] += 1;
				const int indCig = is * zNum + iz;
				curCig [indCig] += sample;
				maskCig [indCig] += 1;
				curImage [iz] += sample;
				maskImage [iz] += 1;
			}
			// add current scattering-angle migration to the main result
			float*  pTo   = dag;
			double* pFrom = curDag;
			for (int id = 0; id < dagSize; ++id, ++pTo, ++pFrom)
				*pTo += *pFrom;
		}

		delete [] travelTimes_;

		sf_warning ("iz %d\n", iz);
	}

	// transfer data from internal angle gather (in double) to the external one (in float)
	int* pMask = maskCig; float* pTo = aCig; double* pFrom = curCig;
	for (int id = 0; id < scatSize; ++id, ++pTo, ++pFrom, ++pMask) {
		*pTo += *pFrom;
		if (*pMask) *pTo /= *pMask;
	}
	// output dip-angle gather normalization
	pMask = maskDag; pTo = dag; 
	for (int id = 0; id < dagSize; ++id, ++pTo, ++pMask) {
		if (*pMask) *pTo /= *pMask;
	}
	// image normalization
	pMask = maskImage; pFrom = curImage; pTo = image;
	for (int iz = 0; iz < zNum; ++iz, ++pTo, ++pFrom, ++pMask) {
		if (*pMask) *pTo = *pFrom / *pMask;		
	}

	delete [] curDag;
	delete [] curCig;
	delete [] curImage;

	delete [] maskDag;
	delete [] maskCig;
	delete [] maskImage;

	return;
}

bool DepthMigrator2D::getSampleByBeam (float curScatAngle, float curDipAngle, float& sample) {

	const int   hNum   = dp_->hNum;
	const float hStart = dp_->hStart; 
	const float hStep  = dp_->hStep; 
	
	const int   xNum   = dp_->xNum;
	const float xStart = dp_->xStart; 
	const float xStep  = dp_->xStep; 

	float baseDir = curDipAngle + 0.5 * curScatAngle;
	const float shiftDir = 0.25 * gp_->scatStep;

	// calc recs lane
	float dir1 = baseDir - shiftDir;
	EscapePoint* escPointRec1 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttRayNum_, dir1, *escPointRec1);

	float dir2 = baseDir + shiftDir;
	EscapePoint* escPointRec2 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttRayNum_, dir2, *escPointRec2);
	
	float recLaneLeft  = escPointRec1->x;
	float recLaneRight = escPointRec2->x;
	if (recLaneLeft > recLaneRight) { float temp = recLaneRight; recLaneRight = recLaneLeft; recLaneLeft = temp; }
	delete escPointRec1;
	delete escPointRec2;

	// calc sources lane
	baseDir = curDipAngle - 0.5 * curScatAngle;
	float dir3 = baseDir - shiftDir;
	EscapePoint* escPointRec3 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttRayNum_, dir3, *escPointRec3);

	float dir4 = baseDir + shiftDir;
	EscapePoint* escPointRec4 = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttRayNum_, dir4, *escPointRec4);			

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
			bool goodRes = this->getRayToPoint (curRecPos, dir1, dir2, timeToRec, recAbsP, onSurf);
			if (!goodRes) continue; // no ray was found
			if (!onSurf) continue; // the ray does not reach the surface
			float timeToSrc (0.f);
			float srcAbsP (0.f);
			goodRes = this->getRayToPoint (curSrcPos, dir3, dir4, timeToSrc, srcAbsP, onSurf);
			if (!goodRes) continue; // no ray was found
			if (!onSurf) continue; // the ray does not reach the surface
			float curTime = (timeToRec + timeToSrc) * 1000; // transform time to "ms"
			sample += this->getSampleFromData (curOffset, 0, curRecPos, curTime, recAbsP);
			++count;
		}
		curRecPos -= xStep;
	}
			
	if (!count) return false; // no sample is returned
	sample /= count;

	return true;
}

// get sample by only ray trace; implemented for zero-offset section only
bool DepthMigrator2D::getSampleByRay (float dipAngle, float& sample) {

	EscapePoint* escPoint = new EscapePoint ();
	this->getEscPointByDirection (travelTimes_, ttRayNum_, dipAngle, *escPoint);

	// check if the ray reaches the daylight surface
	if (!escPoint->isSurf) { return false; }

	// check if the sample is inside the data volume
	if (escPoint->x - dataXMin_ < -1e-4 || escPoint->x - dataXMax_ > 1e-4) { return false; }
	if (escPoint->t < dataTMin_ || escPoint->t > dataTMax_) { return false; }

	// the sample coordinates
	const float sampleX = escPoint->x;
	const float sampleT = escPoint->t;

	const float xStartData = dp_->xStart;
	const float xStepData  = dp_->xStep;

	const int   curXSamp = (sampleX - xStartData) / xStepData;
	const float posLeftTrace = curXSamp * xStepData + xStartData;
	const float posRightTrace = (curXSamp + 1) * xStepData + xStartData;
	
	const float p = escPoint->p;
	const float smallp = p * 0.001;
	const float timeLeft  = 2.f * (sampleT - (sampleX - posLeftTrace)  * smallp); // double time
	const float timeRight = 2.f * (sampleT + (posRightTrace - sampleX) * smallp);	

	const float absP = fabs (p);
	const float zeroOffset = 0.f;
	const float sampleLeft  = this->getSampleFromData (zeroOffset, 0, posLeftTrace,  1000 * timeLeft,  absP); // transform time to "ms"
	const float sampleRight = this->getSampleFromData (zeroOffset, 0, posRightTrace, 1000 * timeRight, absP);
	// interpolation parameters
	const float bef = (sampleX - posLeftTrace) / xStepData;
	const float aft = 1.f - bef;
	
	// get the sample by interpolation
	sample = bef * sampleRight + aft * sampleLeft;

	delete escPoint;
	
	return true;
}

// calculate ray touching the current receiver
bool DepthMigrator2D::getRayToPoint (float curRecPos, float dir1, float dir2, float& timeToPoint, float& pointAbsP, bool& isSurf) {

	const float basedir = dir1 > dir2 ? dir1 : dir2; // basedir "starts" the target beam

	EscapePoint* pTT = travelTimes_; // the first point in travelTimes has the most positive direction !
	float ttDir = pTT->startDir;

	int count (0);
	// get ray inside the target beam
	while (ttDir > basedir && count < ttRayNum_) {
		++pTT; 
		ttDir = pTT->startDir; 
		++count;
	}
	if (count == ttRayNum_) return false; // no ray was found

	float ttX = pTT->x;

	// get two basic points around the target point
	EscapePoint* rightPoint (NULL);
	EscapePoint* leftPoint	(NULL);	
	if (ttX > curRecPos) {
		while (ttX > curRecPos && count < ttRayNum_) {
			++pTT; 
			ttX = pTT->x; 
			++count;
	    }
		if (count == ttRayNum_) return false; // no ray was found
		rightPoint = pTT - 1;
		leftPoint  = pTT; 
	} else {
		while (ttX < curRecPos && count < ttRayNum_) {
			++pTT; 
			ttX = pTT->x; 
			++count;
	    } 
		if (count == ttRayNum_) return false; // no ray was found
		rightPoint = pTT;
		leftPoint  = pTT - 1; 
	}	
	// check if the basic points are valid
	if (!leftPoint->isSurf && !rightPoint->isSurf) {
		isSurf = false; return false;
	} else if (!leftPoint->isSurf) {
		timeToPoint = rightPoint->t;		
		pointAbsP   = rightPoint->p;	
		isSurf = true;	
		return false;
	} else if (!rightPoint->isSurf) {
		timeToPoint = leftPoint->t;		
		pointAbsP   = leftPoint->p;		
		isSurf = true;	
		return false;
	} else {
		const float bef = ( rightPoint->x - curRecPos ) / (rightPoint->x - leftPoint->x);
		const float aft = 1 - bef;
		timeToPoint = leftPoint->t * bef + rightPoint->t * aft;
		pointAbsP   = leftPoint->p * bef + rightPoint->p * aft;
		isSurf = true;	
		return true;
	}	
}

void DepthMigrator2D::getEscPointByDirection (EscapePoint* escPoints, int size, float targetStartDir, EscapePoint &resEscPoint) {

	if (targetStartDir < startDirMin_) {resEscPoint.isSurf = false; return;}
	if (targetStartDir > startDirMax_) {resEscPoint.isSurf = false; return;}

	EscapePoint* pEscPoint = escPoints;
	float curStartDir = pEscPoint->startDir;
	if (curStartDir < targetStartDir) {	resEscPoint = *pEscPoint; return; }

	int count (0);

	while (curStartDir > targetStartDir && count <= size) {
		++pEscPoint;
		curStartDir = pEscPoint->startDir;
		++count;
	}
	if (count == size) {resEscPoint = *(pEscPoint - 1);	return;	}

//	liniar interpolation - NOT the optimal solution

	const float bef = ( pEscPoint->startDir - targetStartDir ) / ( pEscPoint->startDir - (pEscPoint - 1)->startDir );
	const float aft = 1.f - bef;

	resEscPoint.x = pEscPoint->x * aft + (pEscPoint - 1)->x * bef;
	resEscPoint.p = pEscPoint->p * aft + (pEscPoint - 1)->p * bef;
	resEscPoint.t = pEscPoint->t * aft + (pEscPoint - 1)->t * bef;
	resEscPoint.z = pEscPoint->z * aft + (pEscPoint - 1)->z * bef;
	resEscPoint.startDir = pEscPoint->startDir * aft + (pEscPoint - 1)->startDir * bef;
	resEscPoint.offset = pEscPoint->offset * aft + (pEscPoint - 1)->offset * bef;
	resEscPoint.isSurf = resEscPoint.z > 0 ? false : true; // z > 0 - ray is below the daylight surface

	return;
}

void DepthMigrator2D::calcTravelTimes (float curZ, float curX, EscapePoint* escPoints) {

	const float xCIG = ip_->xStart + curX * ip_->xStep;
	wavefrontTracer_.getEscapePoints (xCIG, curZ, escPoints);
		
	return;
}

float DepthMigrator2D::getSampleFromData (const float h, const float geoY, const float geoX1, const float ti, const float p) {
	
	int tNum = dp_->zNum;
	int xNum = dp_->xNum;

	float tStep = dp_->zStep;
	float xStep = dp_->xStep;

	float tStart = dp_->zStart;
	float xStart = dp_->xStart;

//	float geoX = isCMP_ ? geoX1 - h * 0.5 : geoX1;
	float geoX = geoX1 - h;

	const int itMiddle = (int) ((ti - tStart) / tStep);
	if (itMiddle < 0 || itMiddle >= tNum) return 0.f;

	const int xSamp = (int) ((geoX - xStart) / xStep);
	if (xSamp < 0 || xSamp >= xNum) return 0.f;

	// offset

	const int offsetInd = h / dp_->hStep;
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
// transfer parameters to wavefrontTracer
void DepthMigrator2D::setWavefrontTracerParams (int ttRayNum, float ttRayStep, float ttRayStart, 
												int ttNum, float ttStep, float ttStart) {

	wavefrontTracer_.setParams (ttRayNum, ttRayStep, ttRayStart, ttNum, ttStep, ttStart);

	ttRayNum_ = ttRayNum;
	ttRayStep_ = ttRayStep;
	ttRayStart_ = ttRayStart;

	startDirMax_ = ttRayStart_ - 180.f;
	startDirMax_ *= -1; // "-1" is to consist with an agreement
	startDirMin_ = ttRayStart_ +  (ttRayNum_ - 1) * ttRayStep_ - 180.f;
	startDirMin_ *= -1; // "-1" is to consist with an agreement

	return;
}
