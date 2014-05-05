#include <math.h>
#include "tmigratorBase.hh"
#include "curveDefinerBase.hh"
#include "curveDefinerDipOffset3D.hh"
#include "curveDefinerDipOffset.hh"

TimeMigratorBase::TimeMigratorBase () {
}

TimeMigratorBase::~TimeMigratorBase () {
    delete curveDefiner_;
	delete stackTaper_;
}

void TimeMigratorBase::processGather (Point2D& curGatherCoords, float curOffset, const float* const velTrace, const bool isAzDip,
								  float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq) {
}

void TimeMigratorBase::getAzDipFromXY (float curDipY, float curDipX, float& curDip, float&curAz) {

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

float TimeMigratorBase::getMigVel (const float* const velTrace, const float curZ) {

	const int    zNum  = vp_->zNum;
	const float  zStep = vp_->zStep;
	const float zStart = vp_->zStart;
	
	const float zEnd = zStart + (zNum - 1) * zStep;

	if (curZ < zStart || curZ > zEnd)
		return -1;

	const int ind = (int) ((curZ - zStart) / zStep);

	return 0.001 * velTrace[ind];
}

void TimeMigratorBase::setDataLimits () {

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

bool TimeMigratorBase::isPointInsidePoly (Point2D* poly, int nPoly, Point2D& p0) {

	int ip = 0;

	int firstRotation = (int) (( p0.getY() - poly[ip].getY() ) * (poly[ip+1].getX() - poly[ip].getX()) - 
				   ( p0.getX() - poly[ip].getX() ) * (poly[ip+1].getY() - poly[ip].getY()));

	for (int ip = 1; ip < nPoly - 1; ++ip) {
	    int curRotation = (int) (( p0.getY() - poly[ip].getY() ) * (poly[ip+1].getX() - poly[ip].getX()) - 
				     ( p0.getX() - poly[ip].getX() ) * (poly[ip+1].getY() - poly[ip].getY()));
		if (curRotation * firstRotation < 0)
			return false;
	}
	return true;
}

void TimeMigratorBase::initCurveDefiner (bool is3D) {
	
	if (is3D)
		curveDefiner_ = new CurveDefinerDipOffset3D ();
	else
		curveDefiner_ = new CurveDefinerDipOffset ();
	return;
}

void TimeMigratorBase::getStackTaper (const float edgeTaper, const bool isDipAz) {

	return;
}
