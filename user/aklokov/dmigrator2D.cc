#include <rsf.hh>
#include "dmigrator2D.hh"
#include "support.hh"

DepthMigrator2D::DepthMigrator2D () {
}

DepthMigrator2D::~DepthMigrator2D () {
}

void DepthMigrator2D::setWavefrontTracerAxes () {
    wavefrontTracer_.setAxes ();
}

void DepthMigrator2D::processGather (Point2D& curGatherCoords, const float* const data, float* image, float* dag, float* aCig, float* mCig,
				     float* xEsc, float* tEsc) {

    // CONSTANTS

    // inline position
    const float xInd = curGatherCoords.getX ();
    const float curX = ip_->xStart + xInd * ip_->xStep;
    // depth	
    const int   zNum   = gp_->zNum;
    const float zStart = gp_->zStart;
    const float zStep  = gp_->zStep;
    // dip-angle	
    const int   dipNum   = gp_->dipNum;
    const int   dagSize  = dipNum * zNum;
    // scattering-angle
    const int   scatNum  = gp_->scatNum;
    const int   scatSize = scatNum * zNum;
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
    // multi-gather
    memset ( mCig, 0, scatSize * dipNum * sizeof (float) );

    // loop over depth samples
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iz = 0; iz < zNum; ++iz) {	  
	const float curZ = zStart + iz * zStep;		
	if (curZ < velModelDepthMin || curZ > velModelDepthMax)
	    continue;
	this->processDepthSample (curX, curZ, data, curImage + iz, curDag + iz, curCig + iz, mCig + iz, xEsc + iz, tEsc + iz);
    }

    // transfer data from internal gathers (in double) to the external ones (in float)
    float* pTo = aCig; double* pFrom = curCig;
    for (int id = 0; id < scatSize; ++id, ++pTo, ++pFrom) {	
	*pTo = *pFrom;
    }		
    pTo = dag; pFrom = curDag;
    for (int id = 0; id < dagSize; ++id, ++pTo, ++pFrom) {
	*pTo = *pFrom;		
    }
    pTo = image; pFrom = curImage;
    for (int iz = 0; iz < zNum; ++iz, ++pTo, ++pFrom) {
	*pTo = *pFrom;		
    }

    delete [] curDag;
    delete [] curCig;
    delete [] curImage;

    return;
}

void DepthMigrator2D::processDepthSample (const float curX, const float curZ, const float* const data, 
					  double* curImage, double* curDag, double* curCig, float* mCig,
					  float* xEsc, float* tEsc) {

    // CONSTANTS

    const int   zNum   = gp_->zNum;
    // dip-angle	
    const int   dipNum   = gp_->dipNum;
    const float dipStart = gp_->dipStart;
    const float dipStep  = gp_->dipStep;
    // scattering-angle
    const int   scatNum   = gp_->scatNum;
    const float scatStart = gp_->scatStart;
    const float scatStep  = gp_->scatStep;

    const bool zeroStartOffset = dp_->hStart < 1e-6 ? true : false; // start offset in data 

    // ACTION

    // masks for illumination normalization
    int* maskDag = new int [dipNum];
    memset ( maskDag, 0, dipNum * sizeof (int) );
    int* maskCig = new int [scatNum];
    memset ( maskCig, 0, scatNum * sizeof (int) );
    int maskImage (0);

    const float velInPoint = this->getVel (curZ, curX);
    EscapePoint* travelTimes = new EscapePoint [ttRayNum_];		
    this->calcTravelTimes (curZ, curX, travelTimes);
    // travel-time
    for (int it = 0; it < ttRayNum_; ++it) {
	*(xEsc + it * zNum) = travelTimes[it].x;
	*(tEsc + it * zNum) = travelTimes[it].t;
    }

    for (int is = 0; is < scatNum; ++is) {

	const float curScatAngle = scatStart + is * scatStep;
	const float H = 2 * cos (curScatAngle * SF_PI / 360.f) / velInPoint; // 0.5 * SF_PI / 180.f;		

	// loop over dip-angle
	for (int id = 0; id < dipNum; ++id)	{
	    const float curDipAngle = -1 * (dipStart + id * dipStep); // "-1" is to consist with an agreement
	
	    float sample (0.f);
	    bool isGood = this->getSampleByBeam (travelTimes, curScatAngle, curDipAngle, sample);	 
	    if (zeroStartOffset && !isGood && curScatAngle < 1e-6) // implemented for zero-offset only
		isGood = this->getSampleByRay (travelTimes, curDipAngle, sample);
	    if (!isGood)
		continue;
	
	    const float hSample = sample * H;

	    const int dagInd = id * zNum;
	    curDag   [dagInd] += hSample;
	    maskDag [id] += 1;
	    const int scatInd = is * zNum;
	    curCig   [scatInd] += hSample;
	    maskCig [is] += 1;
	    *curImage += hSample;
	    maskImage += 1;
	    const int mInd = (is * dipNum + id) * zNum;
	    mCig [mInd] += hSample;
	}
    }
    // illumination normalization
    double* pRes = curDag; int* pMask = maskDag;
    for (int id = 0; id < dipNum; ++id, pRes += zNum, ++pMask)
	if (*pMask) *pRes /= *pMask;
    pRes = curCig; pMask = maskCig;
    for (int is = 0; is < scatNum; ++is, pRes += zNum, ++pMask)
	if (*pMask) *pRes /= *pMask;
    if (maskImage) *curImage /= maskImage;

    delete [] maskDag;
    delete [] maskCig; 

    delete [] travelTimes;

    return;
}

bool DepthMigrator2D::getSampleByBeam (EscapePoint* travelTimes, float curScatAngle, float curDipAngle, float& sample) {

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
    this->getEscPointByDirection (travelTimes, dir1, *escPointRec1);
    if (!escPointRec1->isSurf) { delete escPointRec1; return false; }

    float dir2 = baseDir + shiftDir;
    EscapePoint* escPointRec2 = new EscapePoint ();
    this->getEscPointByDirection (travelTimes, dir2, *escPointRec2);
    if (!escPointRec2->isSurf) { delete escPointRec1; delete escPointRec2; return false; }
	
    float recLaneLeft  = escPointRec1->x;
    float recLaneRight = escPointRec2->x;
    if (recLaneLeft > recLaneRight) { float temp = recLaneRight; recLaneRight = recLaneLeft; recLaneLeft = temp; }
    delete escPointRec1;
    delete escPointRec2;

    // calc sources lane
    baseDir = curDipAngle - 0.5 * curScatAngle;
    float dir3 = baseDir - shiftDir;
    EscapePoint* escPointRec3 = new EscapePoint ();
    this->getEscPointByDirection (travelTimes, dir3, *escPointRec3);
    if (!escPointRec3->isSurf) { delete escPointRec3; return false; }

    float dir4 = baseDir + shiftDir;
    EscapePoint* escPointRec4 = new EscapePoint ();
    this->getEscPointByDirection (travelTimes, dir4, *escPointRec4);			
    if (!escPointRec4->isSurf)  { delete escPointRec3; delete escPointRec4; return false; }

    float srcLaneLeft  = escPointRec3->x;
    float srcLaneRight = escPointRec4->x;
    if (srcLaneLeft > srcLaneRight) { float temp = srcLaneRight; srcLaneRight = srcLaneLeft; srcLaneLeft = temp; }
    delete escPointRec3;
    delete escPointRec4;

    // loop over receivers
    int count (0);

    float curRecPos = xStart + xStep * (xNum - 1);
    while (curRecPos > recLaneRight) curRecPos -= xStep; // starting from the right

    while (curRecPos > recLaneLeft) {
	// loop over offsets
	for (int ih = 0; ih < hNum; ++ih) {
	    const float curOffset = hStart + ih * hStep;
	    const float curSrcPos = curRecPos - curOffset;
	    if (curSrcPos < srcLaneLeft || curSrcPos > srcLaneRight)
		continue;
	    // get escape point
	    float timeToRec (0.f);
	    float recAbsP   (0.f);
	    bool onSurf;
	    bool goodRes = this->getRayToPoint (travelTimes, curRecPos, dir1, dir2, timeToRec, recAbsP, onSurf);
	    if (!goodRes) continue; // no ray was found
	    if (!onSurf)  continue; // the ray does not reach the surface
	    float timeToSrc (0.f);
	    float srcAbsP   (0.f);
	    goodRes = this->getRayToPoint (travelTimes, curSrcPos, dir3, dir4, timeToSrc, srcAbsP, onSurf);
	    if (!goodRes) continue; // no ray was found
	    if (!onSurf)  continue; // the ray does not reach the surface
	    float curTime = (timeToRec + timeToSrc) * 1000; // transform time to "ms"
	    float curSample (0.f);
	    goodRes = this->getSampleFromData (curOffset, 0, curRecPos, curTime, recAbsP, curSample);
	    if (!goodRes) continue;
	    sample += curSample;
	    ++count;
	}
	curRecPos -= xStep;
    }
			
    if (!count) return false; // no sample is returned

    sample /= count;
    return true;
}

// get sample by only ray trace; implemented for zero-offset section only
bool DepthMigrator2D::getSampleByRay (EscapePoint* travelTimes, float dipAngle, float& sample) {

    EscapePoint* escPoint = new EscapePoint ();
    this->getEscPointByDirection (travelTimes, dipAngle, *escPoint);

    // check if the ray reaches the daylight surface
    if (!escPoint->isSurf) { delete escPoint; return false; }

    // check if the sample is inside the data volume
    if (escPoint->x - dataXMin_ < -1e-4 || escPoint->x - dataXMax_ > 1e-4) { delete escPoint; return false; }
    if (escPoint->t < dataTMin_ || escPoint->t > dataTMax_) { delete escPoint; return false; }

    // the sample coordinates
    const float sampleX = escPoint->x;
    const float sampleT = escPoint->t;

    const float xStartData = dp_->xStart;
    const float xStepData  = dp_->xStep;

    const int   curXSamp = (int) ((sampleX - xStartData) / xStepData);
    const float posLeftTrace = curXSamp * xStepData + xStartData;
    const float posRightTrace = (curXSamp + 1) * xStepData + xStartData;
	
    const float p = escPoint->p;
    delete escPoint;	
    const float smallp = p * 0.001;
    const float timeLeft  = 2.f * (sampleT - (sampleX - posLeftTrace)  * smallp); // double time
    const float timeRight = 2.f * (sampleT + (posRightTrace - sampleX) * smallp);	

    const float absP = fabs (p);
    const float zeroOffset = 0.f;
    float sampleLeft (0.f);
    bool goodRes = this->getSampleFromData (zeroOffset, 0, posLeftTrace,  1000 * timeLeft, absP, sampleLeft); // transform time to "ms"
    if (!goodRes) return false; // bad result
    float sampleRight (0.f);
    goodRes = this->getSampleFromData (zeroOffset, 0, posRightTrace, 1000 * timeRight, absP, sampleRight);
    if (!goodRes) return false; // bad result
    // interpolation parameters
    const float bef = (sampleX - posLeftTrace) / xStepData;
    const float aft = 1.f - bef;
	
    // get the sample by interpolation
    sample = bef * sampleRight + aft * sampleLeft;
	
    return true;
}

// calculate ray touching the current receiver
bool DepthMigrator2D::getRayToPoint (EscapePoint* travelTimes, float curRecPos, float dir1, float dir2, float& timeToPoint, float& pointAbsP, bool& isSurf) {

    const float basedir = dir1 > dir2 ? dir1 : dir2; // basedir "starts" the target beam

    EscapePoint* pTT = travelTimes; // the first point in travelTimes has the most positive direction !
    float ttDir = pTT->startDir;

    int count (0);
    // get ray inside the target beam
    while (ttDir > basedir && count < ttRayNum_) {
	++pTT; ttDir = pTT->startDir; ++count;
    }
    if (count == ttRayNum_) return false; // no ray was found

    // get two basic points around the target point
    EscapePoint* rightPoint (NULL);
    EscapePoint* leftPoint	(NULL);	
    float ttX = pTT->x;
    if (ttX > curRecPos) { // !!! assumtion of regular behaviour inside the beam; may cause errors
	while (ttX > curRecPos && count < ttRayNum_) {
	    ++pTT; ttX = pTT->x; ++count;
	}
	if (count == ttRayNum_) return false; // no ray was found
	rightPoint = pTT - 1;
	leftPoint  = pTT; 
    } else {
	while (ttX < curRecPos && count < ttRayNum_) {
	    ++pTT; ttX = pTT->x; ++count;
	} 
	if (count == ttRayNum_) return false; // no ray was found
	rightPoint = pTT;
	leftPoint  = pTT - 1; 
    }	
    // check if the basic points are valid
    if (!leftPoint->isSurf || !rightPoint->isSurf) {
	isSurf = false; return false;
    }
    const float bef = ( rightPoint->x - curRecPos ) / (rightPoint->x - leftPoint->x);
    const float aft = 1.f - bef;
    timeToPoint = leftPoint->t * bef + rightPoint->t * aft;
    pointAbsP   = leftPoint->p * bef + rightPoint->p * aft;
    isSurf = true;	

    return true;
}

void DepthMigrator2D::getEscPointByDirection (EscapePoint* travelTimes, const float targetStartDir, EscapePoint &resEscPoint) {

    if (targetStartDir < startDirMin_) { resEscPoint.isSurf = false; return; }
    if (targetStartDir > startDirMax_) { resEscPoint.isSurf = false; return; }

    EscapePoint* pEscPoint = travelTimes;
    float curStartDir = pEscPoint->startDir;
    if (curStartDir < targetStartDir)  { resEscPoint.isSurf = false; return; }
    else if (curStartDir == targetStartDir) { // "==" for float ?!
	resEscPoint.x = pEscPoint->x;
	resEscPoint.p = pEscPoint->p;
	resEscPoint.t = pEscPoint->t;
	resEscPoint.z = pEscPoint->z;
	resEscPoint.startDir = pEscPoint->startDir;
	resEscPoint.offset = pEscPoint->offset;
	resEscPoint.isSurf = pEscPoint->isSurf;
    } else {
	int count (0);

	while (curStartDir > targetStartDir && count < ttRayNum_) {
	    ++pEscPoint; curStartDir = pEscPoint->startDir;	++count;
	}
	if (count == ttRayNum_) { resEscPoint.isSurf = false; return; }

	//	liniar interpolation - NOT the optimal solution
	const float bef = ( pEscPoint->startDir - targetStartDir ) / ( pEscPoint->startDir - (pEscPoint - 1)->startDir );
	const float aft = 1.f - bef;

	const EscapePoint* prevPoint = pEscPoint - 1;
	resEscPoint.x = pEscPoint->x * aft + prevPoint->x * bef;
	resEscPoint.p = pEscPoint->p * aft + prevPoint->p * bef;
	resEscPoint.t = pEscPoint->t * aft + prevPoint->t * bef;
	resEscPoint.z = pEscPoint->z * aft + prevPoint->z * bef;
	resEscPoint.startDir = pEscPoint->startDir * aft + prevPoint->startDir * bef;
	resEscPoint.offset = pEscPoint->offset * aft + prevPoint->offset * bef;
	resEscPoint.isSurf = resEscPoint.z > 0 ? false : true; // z > 0 - ray is below the daylight surface
    }

    return;
}

void DepthMigrator2D::calcTravelTimes (float curZ, float xCIG, EscapePoint* escPoints) {

    wavefrontTracer_.getEscapePoints (xCIG, curZ, escPoints);
    return;
}

bool DepthMigrator2D::getSampleFromData (const float h, const float geoY, const float geoX1, const float ti, 
					 const float p, float &sample) {
	
    int tNum = dp_->zNum;
    int xNum = dp_->xNum;

    float tStep = dp_->zStep;
    float xStep = dp_->xStep;

    float tStart = dp_->zStart;
    float xStart = dp_->xStart;

    float geoX (0.f);
    if      (0 == axis2label_) geoX = geoX1 - h;       // if axis2 is "shot"
    else if (1 == axis2label_) geoX = geoX1 - h * 0.5; // if axis2 is "cmp"
    else if (2 == axis2label_) geoX = geoX1;           // if axis2 is "receiver"

    const int itMiddle = (int) ((ti - tStart) / tStep);
    if (itMiddle < 0 || itMiddle >= tNum) return false;

    const int xSamp = (int) ((geoX - xStart) / xStep);
    if (xSamp < 0 || xSamp >= xNum) return false;

    // offset

    const int offsetInd = (int) ((h - dp_->hStart) / dp_->hStep);
    if (offsetInd >= dp_->hNum || offsetInd < 0) return false;

    float* const trace = ptrToData_ + xSamp * tNum + xNum * tNum * offsetInd;

    // middle (main) sample
    
    const float befMiddle = (ti - tStart) * 1.f / tStep - itMiddle;
    const float aftMiddle = 1.f - befMiddle;

    const float sampleMiddle = aftMiddle * trace[itMiddle] + befMiddle * trace[itMiddle + 1];
    
    if (!isAA_) {
	sample = sampleMiddle;
	return true;
    }

    const float filterLength = p * xStep + tStep;
    
    // left sample
 
    const float timeLeft = ti - filterLength;
    const int     itLeft = (int) ((timeLeft - tStart) / tStep); 
	
    if (itLeft < 0) return false;

    const float befLeft = (timeLeft - tStart) * 1.f / tStep - itLeft;
    const float aftLeft = 1.f - befLeft;

    const float sampleLeft = aftLeft   * trace[itLeft]   + befLeft   * trace[itLeft   + 1];

    // right sample
 
    const float timeRight = ti + filterLength;
    const int     itRight = (int) ((timeRight - tStart) / tStep); 

    if (itRight >= tNum - 1) return false;

    const float befRight = (timeRight - tStart) * 1.f / tStep - itRight;
    const float aftRight = 1.f - befRight;

    const float sampleRight = aftRight  * trace[itRight]  + befRight  * trace[itRight  + 1];

    // norm
 
    float imp = tStep / (tStep + timeRight - timeLeft);
    imp *= imp;
    
    sample = (2.f * sampleMiddle - sampleLeft - sampleRight) * imp;
		
    return true;
}
// transfer parameters to wavefrontTracer
void DepthMigrator2D::setWavefrontTracerParams (int ttRayNum, float ttRayStep, float ttRayStart, 
						int ttNum, float ttStep, float ttStart) {

    wavefrontTracer_.setParams (ttRayNum, ttRayStep, ttRayStart + 180.f, ttNum, ttStep, ttStart);
    // "+180.f" is because in wavefrontTracer the vertical direction is 180 degree

    ttRayNum_ = ttRayNum;
    ttRayStep_ = ttRayStep;
    ttRayStart_ = ttRayStart;

    startDirMax_ = ttRayStart_;
    startDirMax_ *= -1; // "-1" is to consist with an agreement
    startDirMin_ = ttRayStart_ +  (ttRayNum_ - 1) * ttRayStep_;
    startDirMin_ *= -1; // "-1" is to consist with an agreement

    return;
}

float DepthMigrator2D::getVel (float curZ, float xCIG) {

    int zInd = (int) ((curZ - vp_->zStart) / vp_->zStep);	
    if (zInd < 0 || zInd >= vp_->zNum) return 1.f;	

    int xInd = (int) ((xCIG - vp_->xStart) / vp_->xStep);
    if (xInd < 0 || xInd >= vp_->xNum) return 1.f;	
		
    return velField_[xInd][zInd];
}
