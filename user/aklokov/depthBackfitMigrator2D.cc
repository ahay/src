#include "depthBackfitMigrator2D.hh"
#include <rsf.hh>
#include "support.hh"
#include <list>

#include <iostream>
using namespace std;

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
				   int sNum, float sStart, float sStep,
				   int izn,  float izo,    float izd,
				   int ixn,  float ixo,    float ixd,
				   float dx, float dt, float xlim, float xapert, int pj,
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

    sNum_   = sNum;
    sStep_  = sStep;
    sStart_ = sStart;	

    izn_    = izn;
    izo_    = izo;
    izd_    = izd;

    ixn_    = ixn;
    ixo_    = ixo;
    ixd_    = ixd;

    dx_   = dx;
    dt_   = dt;
    xlim_ = xlim;
    xapert_ = xapert;
	
    pj_   = pj;

    xVol_ = xVol;
    tVol_ = tVol;

    isAA_ = isAA;

    return;
}

bool DepthBackfitMigrator2D::getSample (float* piData, const float& curX, const float& curZ, const float& curP, float &sample) {

    const int xInd = (int) ((curX - xStart_) / xStep_);
	
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

    list<float> xpnts;
    list<float> zpnts;
	
    iTracer_.traceImage (xVol_, tVol_, curX, curZ, curP, sStart_, &xpnts, &zpnts);

    const int isize = xpnts.size ();
    std::list<float>::iterator iterx = xpnts.begin ();
    std::list<float>::iterator iterz = zpnts.begin ();

    // filter points
    std::list<ImagePoint2D*> goodPoints;
    for (int ir = 0; ir < isize; ++ir, ++iterx, ++iterz) {
	const float lz = *iterz;
	if (lz <= 0) continue; // bad point
	const float lx = *iterx;
	if (fabs (lx - curX) > xapert_) continue;

	if (ir % pj_) continue; // each pj-th point is used

	ImagePoint2D* p = new ImagePoint2D (lx, lz, 0, 0);
	goodPoints.push_back (p);
    }				

    std::list<ImagePoint2D*>::iterator iter = goodPoints.begin ();
    std::list<ImagePoint2D*>::iterator iterLast = goodPoints.end ();
    iterLast--;

    *sample = 0.f;
    int count = 0;
    for (iter = goodPoints.begin (); iter != iterLast;) {
	ImagePoint2D* point1 = *iter;
	++iter;		
	ImagePoint2D* point2 = *iter;
		
	float px1 = point1->x_;		
	float px2 = point2->x_;				

	float pz1 = point1->z_;		
	float pz2 = point2->z_;		

	if (px1 > px2) {
	    float temp (0.f);
	    temp = px2; px2 = px1; px1 = temp; 
	    temp = pz2;	pz2 = pz1; pz1 = temp;
	}

	// slope
	const float dx = px2 - px1;
	if (fabs (dx) > xlim_) continue;
	const float dz = pz2 - pz1;
	curP = -dz / dx;

	// integrate amps
	const int pind = (int) ((px1 - xStart_) / xStep_);
	float px = pind * xStep_ + xStart_;
	if (px - px1 < -1e-6) px += xStep_;

	while (px < px2) {
	    const float bef = (px - px1) / dx;
			
	    float pz = pz1 + bef * dz;

	    float curSample (0.f);
	    bool goodSample = this->getSampleFromImage (piData, px, pz, curP, curSample);			

	    if (!goodSample) {px += xStep_; continue;} // bad sample
			
	    *sample += curSample;
	    ++count;

	    px += xStep_;
	}
    }		

    // every point
    for (iter = goodPoints.begin (); iter != iterLast; ++iter) {
	ImagePoint2D* point1 = *iter;

	float px1 = point1->x_;		
	float pz1 = point1->z_;		

	float curSample (0.f);
	bool goodSample = this->getSample (piData, px1, pz1, curP, curSample);						
	if (goodSample) {
	    *sample += curSample;
	    ++count;
	}
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

    const int pNum = ixn_ * izn_;
#pragma omp parallel for
    for (int ip = 0; ip < pNum; ++ip) {
	const int ix = ip / izn_;
	const int iz = ip % izn_;

	const float curX = ixo_ + ix * ixd_;
	const float curZ = izo_ + iz * izd_;

//		if (curX != 18800 || curZ != 9000) continue;

	float* iPoint = piImage + (ix * izn_ + iz);
	this->getImageSample (piData, curX, curZ, curP, iPoint);
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
