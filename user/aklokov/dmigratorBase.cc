#include <math.h>
#include <rsf.hh>
#include "dmigratorBase.hh"
#include "curveDefinerBase.hh"
#include "curveDefinerDipOffset3D.hh"
#include "curveDefinerDipOffset.hh"
#include "hwt2d.hh"

DepthMigratorBase::DepthMigratorBase () {
}

DepthMigratorBase::~DepthMigratorBase () {
    delete curveDefiner_;
}

void DepthMigratorBase::processGather (Point2D& curGatherCoords, const float* const data, float* gather, float* aCig) {
}

void DepthMigratorBase::getAzDipFromXY (float curDipY, float curDipX, float& curDip, float&curAz) {

    const float curDipYRad = curDipY * CONVRATIO;
    const float curDipXRad = curDipX * CONVRATIO;

	const float tanSqY = pow (tan (curDipYRad), 2);
	const float tanSqX = pow (tan (curDipXRad), 2);

	curDip = atan ( sqrt (tanSqY + tanSqX) ) / CONVRATIO;
	if (tanSqY) {
		curAz  = atan ( tan (curDipXRad) / tan (curDipYRad) ) / CONVRATIO; 	
	} else if (curDipX < 0) {
		curAz = -90.f;
	} else 
		curAz =  90.f;

	if (curDipY < 0) curAz += 180;

	return;
}

float DepthMigratorBase::getMigVel (const float* const velTrace, const float curZ) {

	const int    zNum  = vp_->zNum;
	const float  zStep = vp_->zStep;
	const float zStart = vp_->zStart;
	
	const float zEnd = zStart + (zNum - 1) * zStep;

	if (curZ < zStart || curZ > zEnd)
		return -1;

	const int ind = (int) ((curZ - zStart) / zStep);

	return 0.001 * velTrace[ind];
}

void DepthMigratorBase::setDataLimits () {

	const int zNum = dp_->zNum;
	const int xNum = dp_->xNum;
	const int yNum = dp_->yNum;

    dataZMaxSamp_ = zNum - 1;
    dataXMaxSamp_ = xNum - 1;
    dataYMaxSamp_ = yNum - 1;

    dataZMin_ = dp_->zStart;
    dataXMin_ = dp_->xStart;
    dataYMin_ = dp_->yStart;

    dataZMax_ = dataZMin_ + dataZMaxSamp_ * dp_->zStep;
    dataXMax_ = dataXMin_ + dataXMaxSamp_ * dp_->xStep;
    dataYMax_ = dataYMin_ + dataYMaxSamp_ * dp_->yStep;

	return;
}

bool DepthMigratorBase::isPointInsidePoly (Point2D* poly, int nPoly, Point2D& p0) {

	int ip = 0;

	int firstRotation = (int) ( p0.getY() - poly[ip].getY() ) * (poly[ip+1].getX() - poly[ip].getX()) - 
						( p0.getX() - poly[ip].getX() ) * (poly[ip+1].getY() - poly[ip].getY());

	for (int ip = 1; ip < nPoly - 1; ++ip) {
	    int curRotation = (int) ( p0.getY() - poly[ip].getY() ) * (poly[ip+1].getX() - poly[ip].getX()) - 
			  		      ( p0.getX() - poly[ip].getX() ) * (poly[ip+1].getY() - poly[ip].getY());
		if (curRotation * firstRotation < 0)
			return false;
	}
	return true;
}

void DepthMigratorBase::initCurveDefiner (bool is3D) {
	
	if (is3D)
		curveDefiner_ = new CurveDefinerDipOffset3D ();
	else
		curveDefiner_ = new CurveDefinerDipOffset ();
	return;
}

WavefrontTracer::WavefrontTracer () {
    wp_.tNum = 1001;
    wp_.tStart = 0;
    wp_.tStep = 0.002;
}

WavefrontTracer::~WavefrontTracer () {
}

void WavefrontTracer::getEscapePoints (float xSource, float zSource, EscapePoint* ep) {

    bool rays = true;

    pt2d*    prevWF (NULL); // wavefront it-1 
    pt2d*    curWF  (NULL); // wavefront it   
    pt2d*    nextWF (NULL); // wavefront it+1 

    pt2d     pointPrevWF;            // point  on wft it - 1
    pt2d     pointNextWF;            // point  on wft it + 1 
	pt2d     pointCurWF_left;        // points on wft it  
	pt2d     pointCurWF_center; 
    pt2d     pointCurWF_right;      

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

    // allocate wavefronts
    prevWF = pt2dalloc1 (wp_.rNum);
    curWF  = pt2dalloc1 (wp_.rNum);
    nextWF = pt2dalloc1 (wp_.rNum);

    // initialize wavefronts
    for (int ir = 0; ir < wp_.rNum; ++ir) {
		prevWF[ir].x = curWF[ir].x = nextWF[ir].x = 0.f;
		prevWF[ir].z = curWF[ir].z = nextWF[ir].z = 0.f;
		prevWF[ir].v = curWF[ir].v = nextWF[ir].v = 0.f;
    }

    // init HWT 
    hwt2d_init (velField_, az, ax, at, ag);

    // construct it = 0 wavefront
    int it = 0;
    for (int ir = 0; ir < wp_.rNum; ++ir) {
		prevWF[ir].x = xSource;
		prevWF[ir].z = zSource;
		prevWF[ir].v = hwt2d_getv (prevWF[ir]);
    }

  	// construct it = 1 wavefront
    it = 1;
    for (int ig = 0; ig < wp_.rNum; ++ig) {
		double d, g;
		d =  wp_.tStep  * hwt2d_getv (prevWF[ig]);
		g = (wp_.rStart + ig * wp_.rStep) * SF_PI / 180;

		curWF[ig].x = xSource + d * sin(g);
		curWF[ig].z = zSource + d * cos(g);
		curWF[ig].v = hwt2d_getv (curWF [ig]);
    }

    // loop over time 
    for (it = 2; it < wp_.tNum; ++it) {
		float curTime = wp_.tStart + (it - 2) * wp_.tStep;

		if (wp_.rNum > 3 && !rays) {
	 		// boundaries 
		    int ig = 0;     
		    nextWF[ig] = hwt2d_raytr (prevWF[ig], curWF [ig]);
		    ig = wp_.rNum - 1; 
		    nextWF[ig] = hwt2d_raytr (prevWF[ig], curWF [ig]);

	  		for (int ig = 1; ig < wp_.rNum - 1; ++ig) {
		
				pointCurWF_left = curWF [ig - 1];
				pointCurWF_center = curWF [ig];
				pointNextWF = prevWF[ig];
				pointCurWF_right = curWF [ig + 1];
		
				if (hwt2d_cusp (pointNextWF, pointCurWF_left, 
								pointCurWF_center, pointCurWF_right)) {
				    pointPrevWF = hwt2d_raytr (pointNextWF,
				    							pointCurWF_center   );
				} else {
				    pointPrevWF = hwt2d_wfttr (pointNextWF, 
				    							pointCurWF_left, 
				    							pointCurWF_center, 
				    							pointCurWF_right);
				}
				nextWF[ig] = pointPrevWF;
			  }
		} else {
		    for (int ig = 0; ig < wp_.rNum; ++ig) {
				pointCurWF_center = curWF [ig];  
				pointNextWF = prevWF[ig];
				pointPrevWF = hwt2d_raytr (pointNextWF, pointCurWF_center);
				nextWF[ig] = pointPrevWF;
	    	}

			// checking if the ray has reached the daylight surface
	
		    for (int ig = 0; ig < wp_.rNum; ++ig) {
				if (ep[ig].isSurf) continue;
				else if (prevWF[ig].z >= 0) {
					ep[ig].x = prevWF[ig].x;				
					ep[ig].z = prevWF[ig].z;				
					ep[ig].t = curTime;					
				} else if (prevWF[ig].z < 0 && !ep[ig].isSurf) {
								
					float bef = (1 - prevWF[ig].z / (prevWF[ig].z - ep[ig].z));
					float dx = prevWF[ig].x - ep[ig].x;
					
					ep[ig].x += bef * dx;
					ep[ig].t += bef * wp_.tStep;
//					ep[ig].t *= 2; // because we need double time
							
					ep[ig].isSurf = true;
				}
				ep[ig].startDir = wp_.rStart + ig * wp_.rStep - 180.f;
				ep[ig].startDir *= -1;
			}
		}	

		// step in time
		for (int ir = 0; ir < wp_.rNum; ++ir) {
		    prevWF[ir] = curWF [ir];
		    curWF [ir] = nextWF[ir];
		}
    } 
    
    
    for (int ig = 1; ig < wp_.rNum - 1; ++ig) {
		const float dx = ep[ig + 1].x - ep[ig - 1].x;
		if (fabs (dx) < 1e-6) {ep[ig].p = 0.f; continue;}
		const float dt = ep[ig + 1].t - ep[ig - 1].t;
		ep[ig].p = fabs (1000 * dt / dx);
	}
		
	ep[0].p = ep[1].p;
	ep[wp_.rNum - 1].p = ep[wp_.rNum - 2].p;
       
    return;
}
