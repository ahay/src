#include <rsf.hh>
#include "dmigratorBase.hh"
#include "hwt2d.hh"

DepthMigratorBase::DepthMigratorBase () {
}

DepthMigratorBase::~DepthMigratorBase () {
}

void DepthMigratorBase::processGather (Point2D& curGatherCoords, const float* const data, 
									   float* image, float* gather, float* aCig, float* mCig, float* xEsc, float* tEsc) {
	return;
}

void DepthMigratorBase::setWavefrontTracerParams (int ttRayNum, float ttRayStep, float ttRayStart,
												  int ttNum, float ttStep, float ttStart) {
	return;
}

void DepthMigratorBase::setVelModel (float** velField) {
	return;
}

void DepthMigratorBase::setWavefrontTracerAxes () {
	return;
}

void DepthMigratorBase::setDataLimits () {

	const int tNum = dp_->zNum;
	const int xNum = dp_->xNum;
	const int yNum = dp_->yNum;

    dataTMaxSamp_ = tNum - 1;
    dataXMaxSamp_ = xNum - 1;
    dataYMaxSamp_ = yNum - 1;

    dataTMin_ = dp_->zStart;
    dataXMin_ = dp_->xStart;
    dataYMin_ = dp_->yStart;

    dataTMax_ = dataTMin_ + dataTMaxSamp_ * dp_->zStep;
    dataXMax_ = dataXMin_ + dataXMaxSamp_ * dp_->xStep;
    dataYMax_ = dataYMin_ + dataYMaxSamp_ * dp_->yStep;

	return;
}

///////////////////////////
// Class WavefrontTracer //
///////////////////////////

WavefrontTracer::WavefrontTracer () {
}

WavefrontTracer::~WavefrontTracer () {
}

void WavefrontTracer::setAxes () {
	// Cartesian coordinates 
	// z axis 
    sf_axis az = sf_maxa (vp_.zNum, vp_.zStart, vp_.zStep);
    sf_setlabel (az, "z");
	// x axis 
    sf_axis ax = sf_maxa (vp_.xNum, vp_.xStart, vp_.xStep);
    sf_setlabel (ax, "x");
	// ray coordinates 
	// time axis 
    sf_axis at = sf_maxa (wp_.tNum, wp_.tStart, wp_.tStep);
    sf_setlabel (at, "t");
	// shooting angle axis 
    sf_axis ag = sf_maxa (wp_.rNum, wp_.rStart, wp_.rStep);
    sf_setlabel (ag, "g");

    hwt2d_init (pVelField_, az, ax, at, ag);

	sf_maxa_free (az);
	sf_maxa_free (ax);
	sf_maxa_free (at);	
	sf_maxa_free (ag);
}

void WavefrontTracer::getEscapePoints (float xSource, float zSource, EscapePoint* ep) {

	// CONSTANTS
	const int   rNum   = wp_.rNum;
	const float rStep  = wp_.rStep;
	const float rStart = wp_.rStart;

	const int   tNum   = wp_.tNum;
	const float tStep  = wp_.tStep;
	const float tStart = wp_.tStart;

    pt2d*    prevWF (NULL); // wavefront it-1 
    pt2d*    curWF  (NULL); // wavefront it   
    pt2d*    nextWF (NULL); // wavefront it+1 

    pt2d     pointPrevWF;            // point  on wft it - 1
    pt2d     pointNextWF;            // point  on wft it + 1 
	pt2d     pointCurWF_center; 

    // allocate wavefronts
    prevWF = pt2dalloc1 (rNum);
    curWF  = pt2dalloc1 (rNum);
    nextWF = pt2dalloc1 (rNum);

    // initialize wavefronts
    for (int ir = 0; ir < rNum; ++ir) {
		prevWF[ir].x = curWF[ir].x = nextWF[ir].x = 0.f;
		prevWF[ir].z = curWF[ir].z = nextWF[ir].z = 0.f;
		prevWF[ir].v = curWF[ir].v = nextWF[ir].v = 0.f;
    }

    // construct it = 0 wavefront
	pt2d* pWF = prevWF;
    int it = 0;
    for (int ir = 0; ir < rNum; ++ir, ++pWF) {
		pWF->x = xSource;
		pWF->z = zSource;
		pWF->v = hwt2d_getv ( *pWF );
    }

  	// construct it = 1 wavefront
	pWF = curWF;    
	it = 1;
    for (int ir = 0; ir < rNum; ++ir, ++pWF) {
		const double d =  wp_.tStep  * hwt2d_getv ( prevWF[ir] );
		const double g = (rStart + ir * rStep) * SF_PI / 180.f;
		pWF->x = xSource + d * sin(g);
		pWF->z = zSource + d * cos(g);
		pWF->v = hwt2d_getv ( *pWF );
    }

    // loop over time 
    for (it = 2; it < tNum; ++it) {
		const float curTime = tStart + (it - 2) * tStep;
	    for (int ir = 0; ir < rNum; ++ir) {
			pointCurWF_center = curWF [ir];  
			pointNextWF = prevWF[ir];
			pointPrevWF = hwt2d_raytr (pointNextWF, pointCurWF_center);
			nextWF[ir] = pointPrevWF;
    	}
		// checking if the ray has reached the daylight surface
		EscapePoint* pPoint = ep;
	    for (int ir = 0; ir < rNum; ++ir, ++pPoint) {
			const float prevZ = prevWF[ir].z; 							
			if (pPoint->isSurf) continue; // the ray has already reached the surface
			if (prevZ < 0) { // the ray has left the model
				const float bef = (1 - prevZ / (prevZ - pPoint->z));
				const float dx  = prevWF[ir].x - pPoint->x;
				pPoint->x += bef * dx;
				pPoint->t += bef * tStep;
				pPoint->z = 0.f;
				pPoint->isSurf = true;
			} else { // the ray is inside the model	- update escape point	 
				pPoint->x = prevWF[ir].x;				
				pPoint->z = prevZ;				
				pPoint->t = curTime;					
			}
		}
		// step in time
		pt2d* pPrev = prevWF; pt2d* pCur = curWF; pt2d* pNext = nextWF;
		for (int ir = 0; ir < rNum; ++ir, ++pPrev, ++pCur, ++pNext) {
			*pPrev = *pCur;
			*pCur = *pNext;
		}
    } 

	EscapePoint* pPoint = ep + 1; 
	const int rNumReduced = rNum - 1;   
	for (int ir = 1; ir < rNumReduced; ++ir, ++pPoint) {
		// start direction
		pPoint->startDir = rStart + ir * rStep - 180.f;
		pPoint->startDir *= -1;
		// p
		const float dx = (pPoint + 1)->x - (pPoint - 1)->x;
		if (fabs (dx) < 1e-6) { pPoint->p = 0.f; continue; } // vertical ray
		const float dt = (pPoint + 1)->t - (pPoint - 1)->t;
		pPoint->p = fabs (1000 * dt / dx);
	}
	// handle edge points			
	ep[0].p = ep[1].p;
	ep[rNum - 1].p = ep[rNum - 2].p;
 
	ep[0].startDir = rStart - 180.f;
	ep[0].startDir *= -1;
	ep[rNum - 1].startDir = rStart + (rNum - 1) * rStep - 180.f;
    ep[rNum - 1].startDir *= -1;

	pPoint = ep;
    for (int ir = 0; ir < rNum; ++ir, ++pPoint) {
		if (!pPoint->isSurf) pPoint->t = -1.f;
	}

	pt2dfree1 (prevWF);
	pt2dfree1 (curWF);
	pt2dfree1 (nextWF);

    return;
}
