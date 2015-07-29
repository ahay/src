#include <math.h>
#include "curveDefinerDipOffset.hh"
#include "support.hh"

CurveDefinerDipOffset::CurveDefinerDipOffset () {
					
}

CurveDefinerDipOffset::~CurveDefinerDipOffset () {
}

void CurveDefinerDipOffset::getEscapePoint (const float curY, const float curX, const float curZeroTime, 
							 			   const float curDip, const float curAz, const float migVel,
										   float& time, float& xRayExit, float& yRayExit) {

	if (curZeroTime < 1e-6) {
		xRayExit = curX;
		yRayExit = curY;
		time     = 0.f;
		return;
	}
	
	if (fabs(curDip - 90) < 1e-3) {
		xRayExit = curX;
		yRayExit = curY;
		time     = curZeroTime;
		return;
	}
									   
  	const float dipInRad = curDip * CONVRATIO;
    const float azInRad  = curAz  * CONVRATIO;

	// escape position

	const float rho = 2 * curOffset_ / (migVel * curZeroTime);
	const float tanAlpha2 = pow (tan (dipInRad), 2);
	
	const float a = tanAlpha2 * rho;
	const float b = 2 * (tanAlpha2 + 1);
	const float c = -rho;
	
	const float D = b*b - 4*a*c;
	
	const float tanBeta = a ? ( -b + sqrt (D) ) / (2 * a) : -c / b;
	// const float tanBeta2 = a ? ( -b - sqrt (D) ) / (2 * a) : -c / b;

	const float beta = atan (tanBeta) * 180 / 3.1415;
	const float dirAngle = beta + curDip;

	if (fabs (dirAngle - 90) < 1e-3) {
		xRayExit = curX;
		yRayExit = curY;
		time = curZeroTime;
		return;	
	}

	const float dirAngleInRad = dirAngle * CONVRATIO;
		
	const float depth     = curZeroTime * migVel  / 2.f;	
	const float shiftToRec = depth * tan (dirAngleInRad);	

    xRayExit = curX + shiftToRec * sin (azInRad);
    yRayExit = curY + shiftToRec * cos (azInRad);

	// escape time
	const float distToRec = sqrt (pow (shiftToRec, 2) + depth*depth);	

	const float srcPos = xRayExit - curOffset_;
	const float distToSrc = sqrt (pow (srcPos - curX, 2) + depth*depth);

	time = (distToSrc + distToRec) / migVel;

	return;
}

float CurveDefinerDipOffset::getEscapeTime (const float dy, const float dx, const float t0, const float migVel, float& p) {

	const float z2 = pow (t0 * migVel * 0.5, 2);

	const float distToRec = sqrt (z2 + dx * dx);	
	const float distToSrc = sqrt (z2 + pow (curOffset_ - dx, 2));

	float time = (distToRec + distToSrc) / migVel;

	p = (distToRec && distToSrc) ? fabs ((dx / distToRec - (curOffset_ - dx) / distToSrc) / migVel) : 0.f;
	
	return time;
}
