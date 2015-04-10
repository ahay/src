#include <string.h>
#include "vsptmigrator2D.hh"
#include "support.hh"
#include "curveDefinerBase.hh"
#include "curveDefinerDipOffset.hh"
#include <rsf.hh>
#ifdef _OPENMP
#include <omp.h>
#endif

VSPTimeMigrator2D::VSPTimeMigrator2D () {
}

VSPTimeMigrator2D::~VSPTimeMigrator2D () {
}

void VSPTimeMigrator2D::processGather (Point2D& curGatherCoords, const float* const data, float* image, float* dag, float* cig, float* mCig, 
				       float* xEsc, float* tEsc) {
   
    const int   tNum     = ip_->zNum;
    const float tStart   = ip_->zStart;
    const float tStep    = ip_->zStep;

    const int   curX     = (int) curGatherCoords.getX ();    
    const float xCIG     = ip_->xStart + curX * ip_->xStep;

#pragma omp parallel for
    for (int it = 0; it < tNum; ++it) {	  
	const float curT = tStart + it * tStep;		
	this->processTimeSample (xCIG, curT, data, image + it, dag + it, cig + it);
    }
   
    return;
}

float VSPTimeMigrator2D::getSampleByRay (const float yCIG, const float xCIG, const float curZeroTime, 
					 const float curDip, const float curAz, const float migVel, const bool isAzDip,
					 const float yEmAngle, const float xEmAngle) {
    
    return 0.f;
}

bool VSPTimeMigrator2D::getSampleFromData (const float geoS, const float geoX, const float ti, 
					   const float p, float &sample) {
	
    int tNum = dp_->zNum;
    int xNum = dp_->xNum;
    int hNum = dp_->hNum;

    float tStep = dp_->zStep;
    float xStep = dp_->xStep;
    float hStep = dp_->hStep;

    float tStart = dp_->zStart;
    float xStart = dp_->xStart;
    float hStart = dp_->hStart;


    const int itMiddle = (int) ((ti - tStart) / tStep); 
    if (itMiddle < 0 || itMiddle >= tNum) return false;

    const int xSamp = (int) ((geoX - xStart) / xStep);
    if (xSamp < 0 || xSamp >= xNum) return false;

    const int sSamp = (int) ((geoS - hStart) / hStep);
    if (sSamp < 0 || sSamp >= hNum) return false;

    float* const trace = ptrToData_ + (xSamp + sSamp * xNum) * tNum;

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

float VSPTimeMigrator2D::getVel (float curZ, float xCIG) {

    int zInd = (int) ((curZ - vp_->zStart) / vp_->zStep);	
    if (zInd < 0 || zInd >= vp_->zNum) return 1.f;	

    int xInd = (int) ((xCIG - vp_->xStart) / vp_->xStep);
    if (xInd < 0 || xInd >= vp_->xNum) return 1.f;	
		
    return velField_[xInd][zInd];
}

void VSPTimeMigrator2D::processTimeSample (const float curX, const float curT, const float* const data, 
					   float* curImage, float* curDag, float* curCig) {

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

    // ACTION

    // masks for illumination normalization
    int* maskDag = new int [dipNum];
    memset ( maskDag, 0, dipNum * sizeof (int) );
    int* maskCig = new int [scatNum];
    memset ( maskCig, 0, scatNum * sizeof (int) );
    int maskImage (0);

    const float velInPoint = this->getVel (curT, curX);

    const float curZ = 0.5 * curT * velInPoint / 1000; // imaged point position

    for (int is = 0; is < scatNum; ++is) {

	const float curScatAngle = scatStart + is * scatStep;
	const float H = 2 * cos (curScatAngle * SF_PI / 360.f) / velInPoint; // 0.5 * SF_PI / 180.f;		

	// loop over dip-angle
	for (int id = 0; id < dipNum; ++id)	{
	    const float curDipAngle = dipStart + id * dipStep;
	
	    float sample (0.f);
	    bool isGood = this->getSampleByBeam (curX, curZ, velInPoint, curScatAngle, curDipAngle, sample);	 
//			if (zeroStartOffset && !isGood && curScatAngle < 1e-6) // implemented for zero-offset only
//				isGood = this->getSampleByRay (travelTimes, curDipAngle, sample);
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
	}
    }
    // illumination normalization
    float* pRes = curDag; int* pMask = maskDag;
    for (int id = 0; id < dipNum; ++id, pRes += zNum, ++pMask)
	if (*pMask) *pRes /= *pMask;
    pRes = curCig; pMask = maskCig;
    for (int is = 0; is < scatNum; ++is, pRes += zNum, ++pMask)
	if (*pMask) *pRes /= *pMask;
    if (maskImage) *curImage /= maskImage;

    delete [] maskDag;
    delete [] maskCig; 

    return;
}

bool VSPTimeMigrator2D::getSampleByBeam (float curX, float curZ, float vel, float curScatAngle, float curDipAngle, float& sample) {

    const int   hNum   = dp_->hNum;
    const float hStart = dp_->hStart; 
    const float hStep  = dp_->hStep; 
	
    const int   xNum   = dp_->xNum;
    const float xStart = dp_->xStart; 
    const float xStep  = dp_->xStep; 

    float baseDir = curDipAngle - 0.5 * curScatAngle;
    const float shiftDir = 0.25 * gp_->scatStep;

    // calc recs lane
    float dir1 = baseDir - shiftDir;
    float r1 = 0.f;
    bool res = this->getReceiverByDirection (curX, curZ, dir1, r1);
    if (!res) return false;

    float dir2 = baseDir + shiftDir;
    float r2 = 0.f;
    res = this->getReceiverByDirection (curX, curZ, dir2, r2);
    if (!res) return false;

    float recLaneTop    = r1;
    float recLaneBottom = r2;
    if (recLaneTop > recLaneBottom) { float temp = recLaneTop; recLaneTop = recLaneBottom; recLaneBottom = temp; }
	
    // calc sources lane
    baseDir = curDipAngle + 0.5 * curScatAngle;
    float dir3 = baseDir - shiftDir;
    float s1 = 0.f;
    res = this->getSourceByDirection (curX, curZ, dir3, s1);	
    if (!res) return false;

    float dir4 = baseDir + shiftDir;
    float s2 = 0.f;
    res = this->getSourceByDirection (curX, curZ, dir4, s2);		
    if (!res) return false;

    float srcLaneLeft  = s1;
    float srcLaneRight = s2;
    if (srcLaneLeft > srcLaneRight) { float temp = srcLaneRight; srcLaneRight = srcLaneLeft; srcLaneLeft = temp; }

    // loop over receivers
    int count (0);

    float curRecPos = xStart + xStep * (xNum - 1);
    while (curRecPos > recLaneBottom) curRecPos -= xStep; // starting from the bottom

    while (curRecPos > recLaneTop) {
	// loop over offsets
	for (int ih = 0; ih < hNum; ++ih) {
	    const float curOffset = hStart + ih * hStep;
	    const float curSrcPos = curOffset; // position of the well is zero
	    if (curSrcPos < srcLaneLeft || curSrcPos > srcLaneRight)
		continue;
	    // get escape point
	    float timeToRec = sqrt ( pow (curX, 2) + pow (curZ - curRecPos, 2) ) / vel;

	    // for anti-aliasing
	    float dz = 1;
	    float time1 = sqrt ( pow (curX, 2) + pow (curZ - curRecPos + dz, 2) ) / vel;
	    float time2 = sqrt ( pow (curX, 2) + pow (curZ - curRecPos - dz, 2) ) / vel;
	    float p = fabs ( (time1 - time2) / (2 * dz) );

	    float timeToSrc = sqrt ( pow (curX - curSrcPos, 2) + pow (curZ, 2) ) / vel;

	    float curTime = (timeToRec + timeToSrc) * 1000; // transform time to "ms"
	    float curSample (0.f);
	    bool goodRes = this->getSampleFromData (curSrcPos, curRecPos, curTime, p, curSample);

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

bool VSPTimeMigrator2D::getReceiverByDirection (float curX, float curZ, float dir2, float& rec) {

    if (dir2 == 0) return false; // vertical direction
    const float tang = tan ( dir2 * M_PI / 180.f );	
    if (tang > 0) return false; // wrong direction - not to the receivers
	
    // get receiver
    rec = curZ + curX / tang;
	
    // limit checking
    const float rTop = dp_->xStart;
    if (rec < rTop) return false;
    const float rBot = dp_->xStart + dp_->xStep * (dp_->xNum - 1);
    if (rec > rBot) return false;

    return true; // good receiver
}

bool VSPTimeMigrator2D::getSourceByDirection (float curX, float curZ, float dir2, float& src) {

    if (dir2 == 0) return false;

    src = curX + curZ * tan ( dir2 * M_PI / 180.f );	
	
    if (src < 0) return false; // temporary condition

    return true;
}
