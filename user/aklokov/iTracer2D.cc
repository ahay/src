#include "iTracer2D.hh"
#include <rsf.hh>
#include "support.hh"
#include <algorithm>
#include <iostream>

static bool deleteAll (ImagePoint2D* p) { delete p; return true; }
static bool predPos   (const ImagePoint2D* lhs, const ImagePoint2D* rhs) { return lhs->x_ < rhs->x_; }
static bool predTime  (const ImagePoint2D* lhs, const ImagePoint2D* rhs) { return lhs->z_ < rhs->z_; }

// -- class ImagePoint2D ---
ImagePoint2D::ImagePoint2D () : x_ (0.f), z_ (0.f),
    ix_ (0), iz_ (0) {
}

ImagePoint2D::ImagePoint2D (float x, float z, int ix, int iz) : x_ (x), z_ (z), ix_ (ix), iz_ (iz) {
}

ImagePoint2D::~ImagePoint2D () {

}

ImagePoint2D& ImagePoint2D::operator= (const ImagePoint2D& point) {
    if (this == &point)
        return *this;
    x_ = point.x_;
    z_ = point.z_;
    ix_ = point.ix_;
    iz_ = point.iz_;

    return *this;
}

//bool ImagePoint2D::operator< (const ImagePoint2D& p) const {
//	return (iz_ < p.iz_);
//}

ITracer2D::ITracer2D () {
}

ITracer2D::~ITracer2D () {
}

void ITracer2D::init (int zNum, float zStart, float zStep, 
		      int pNum, float pStart, float pStep,
		      int xNum, float xStart, float xStep,
		      float dx, float dt) {
	
    zNum_   = zNum;
    zStep_  = zStep;
    zStart_ = zStart;

    pNum_   = pNum;
    pStep_  = pStep;
    pStart_ = pStart;

    xNum_   = xNum;
    xStep_  = xStep;
    xStart_ = xStart;	

    dx_     = dx;
    dt_     = dt;

    return;
}

bool ITracer2D::isPointInsideTriangle (float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {

    double x01 = x0 - x1;
    double y01 = y0 - y1;
    if (fabs (x01) < 1e-6 && fabs (y01) < 1e-6) return true;

    double x02 = x0 - x2;
    double y02 = y0 - y2;
    if (fabs (x02) < 1e-6 && fabs (y02) < 1e-6) return true;

    double x03 = x0 - x3;
    double y03 = y0 - y3;
    if (fabs (x03) < 1e-6 && fabs (y03) < 1e-6) return true;

    //
		
    double x31 = x3 - x1; 
    double y31 = y3 - y1;

    double x21 = x2 - x1;
    double y21 = y2 - y1;

    // dot products

    double dot00 = x31*x31 + y31*y31;
    double dot01 = x31*x21 + y31*y21;
    double dot02 = x31*x01 + y31*y01;
    double dot11 = x21*x21 + y21*y21;
    double dot12 = x21*x01 + y21*y01;

    // Compute barycentric coordinates
    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

void ITracer2D::traceImage (float* xVol, float* tVol, float x0, float z0, float p0, float sa, list<float>* xRes, list<float>* zRes) {

    const int zInd = (int) ((z0 - zStart_) / zStep_);
    // const int pInd = (p0 - pStart_) / pStep_;
    const int xInd = (int) ((x0 - xStart_) / xStep_);

    const int xRed = xNum_ - 1;
    const int zRed = zNum_ - 1;

    const float halfScatNum = 0.5 * sa;

    // escape line for the diffraction point
    list<ImagePoint2D*> escPoints;
    for (int ip = 0; ip < pNum_; ++ip) {
	const float curP = pStart_ + ip * pStep_;
	const float p1 = curP - halfScatNum;
	const int ip1 = (int) ((p1 - pStart_) / pStep_);
	if (ip1 < 0 || ip1 > pNum_ - 1) continue;
	const int pind1 = (xInd * pNum_ + ip1) * zNum_ + zInd;		
	const float t1 = tVol [pind1];
	if (t1 < 0) continue; // this is a "bad" escape point

	const float p2 = curP + halfScatNum;
	const int ip2 = (int) ((p2 - pStart_) / pStep_);
	if (ip2 < 0 || ip2 > pNum_ - 1) continue;
	const int pind2 = (xInd * pNum_ + ip2) * zNum_ + zInd;				
	const float t2 = tVol [pind2];
	if (t2 < 0) continue; // this is a "bad" escape point
		
	ImagePoint2D* p = new ImagePoint2D (xVol [pind2], t1 + t2, ip, 0); // two-way time
	escPoints.push_back (p);
    }

    // sort escape points by position
    escPoints.sort (predPos);

    // constant-dip panel extraction
    const int panelSize = zNum_ * xNum_;	
    float* xPanel = sf_floatalloc (panelSize);
    float* tPanel = sf_floatalloc (panelSize);
    memset ( xPanel, 0, panelSize * sizeof (float) );
    memset ( tPanel, 0, panelSize * sizeof (float) );

    list<ImagePoint2D*> allPoints;

    float* pXPanel = xPanel;
    float* pTPanel = tPanel;

    const float p1 = p0 - halfScatNum;
    const float p2 = p0 + halfScatNum;

    for (int ix = 0; ix < xNum_; ++ix) {
	for (int iz = 0; iz < zNum_; ++iz, ++pXPanel, ++pTPanel) {

	    const int ip1 = (int) ((p1 - pStart_) / pStep_);
	    if (ip1 < 0 || ip1 > pNum_ - 1) continue;
	    const int pind1 = (ix * pNum_ + ip1) * zNum_ + iz;		
	    float t1 = tVol [pind1];

	    const int ip2 = (int) ((p2 - pStart_) / pStep_);
	    if (ip2 < 0 || ip2 > pNum_ - 1) continue;
	    const int pind2 = (ix * pNum_ + ip2) * zNum_ + iz;		
	    float t2 = tVol [pind2];

	    *pXPanel = xVol [pind2];
	    *pTPanel = t1 + t2;

	    ImagePoint2D* p = new ImagePoint2D (*pXPanel, *pTPanel, ix, iz);
	    allPoints.push_back (p);
	}
    }

    // sort image points by position
    allPoints.sort (predPos);

    const float minX = allPoints.front()->x_;
    const float maxX = allPoints.back()->x_;	

    // loop over escape points
    list<ImagePoint2D*>::iterator iterCurX = allPoints.begin ();
    list<ImagePoint2D*>::iterator iterEP;
    for (iterEP = escPoints.begin (); iterEP != escPoints.end(); ++iterEP) {

	ImagePoint2D* escPoint = *iterEP;

	const float curX = escPoint->x_;
	if (curX < minX || curX > maxX) continue;		
	
	const float curT = escPoint->z_;

	const float x1 = curX - dx_;
	const float x2 = curX + dx_;

	list<ImagePoint2D*>::iterator iter = iterCurX;

	bool found1 (false);
	list<ImagePoint2D*>::iterator iterMin;
	while (!found1 && iter != allPoints.end()) {
	    ImagePoint2D* p1 = *iter;
	    if (p1->x_ > x1) { iterMin = iter; found1 = true; }
	    ++iter;
	}
	if (!found1) continue;
	iterCurX = iterMin;
	bool found2 = false;
	list<ImagePoint2D*>::iterator iterMax;
	while ( !found2 && iter != allPoints.end () ) {
	    ImagePoint2D* p2 = *iter;
	    if (p2->x_ > x2) {iterMax = iter; found2 = true; }
	    ++iter;
	}
	if (!found2) {
	    iterMax = iter;
	}

	// get good x-points
	list<ImagePoint2D*> goodXPoints;
	goodXPoints.insert (goodXPoints.end(), iterMin, iterMax);
		
	// narrow by times
	const float t1 = curT - dt_;
	const float t2 = curT + dt_;

	goodXPoints.sort (predTime);

	iter = goodXPoints.begin ();

	found1 = false;
	while (!found1 && iter != goodXPoints.end() ) {
	    ImagePoint2D* p1 = *iter;
	    if (p1->z_ > t1) { iterMin = iter; found1 = true; }
	    ++iter;
	}
	if (!found1) continue;

	found2 = false;
	while ( !found2 && iter != goodXPoints.end() ) {
	    ImagePoint2D* p2 = *iter;
	    if (p2->z_ > t2) {iterMax = iter; found2 = true; }
	    ++iter;
	}
	if (!found2) {
	    iterMax = iter;
	}

	// get good points
	list<ImagePoint2D*> goodPoints;
	goodPoints.insert (goodPoints.end(), iterMin, iterMax);
		
	bool isFound = false;
	for (iter = goodPoints.begin(); iter != goodPoints.end(); ++iter) {
	    ImagePoint2D* iPoint = *iter;			
	    int ix = iPoint->ix_;			
	    int iz = iPoint->iz_;			
	
	    float foundX (0.f);	
	    float foundZ (0.f);	
	    float foundDist (FLT_MAX);	
	
	    float x (0.f);
	    float z (0.f);
	    float dist (FLT_MAX);

	    if (ix && iz) {	
		const int mode = -1; // upper triangle
		isFound = this->checkTriangle (curX, curT, ix, iz, mode, xPanel, tPanel, &x, &z, dist);
		if (isFound) {
//					xRes->push_back (x);
//					zRes->push_back (z);
		}
		if (isFound && dist < foundDist) {
		    foundX = x;
		    foundZ = z;
		    foundDist = dist;
		}
	    }

	    if (ix < xRed && iz < zRed) {
		const int mode = 1; // lower triangle
		isFound = this->checkTriangle (curX, curT, ix, iz, mode, xPanel, tPanel, &x, &z, dist);
		if (isFound) {
//					xRes->push_back (x);
//					zRes->push_back (z);
		}
		if (isFound && dist < foundDist) {
		    foundX = x;
		    foundZ = z;
		    foundDist = dist;
		}
	    }

	    if (foundDist < FLT_MAX) {
		xRes->push_back (foundX);
		zRes->push_back (foundZ);
	    }
	}
    }		

    // remove duplicate values 
    xRes->unique ();
    zRes->unique ();	

    // the following part of the code was written to remove dublicate values from the UNsorted lists
    // However, we may see that the source point dublicates only and, moreover, 
    // this values go one by one.

    // therefore, I use list::unique ()

/*	list<float>::iterator iterZ0 = zRes->begin ();
	for (list<float>::iterator iterX0 = xRes->begin (); iterX0 != xRes->end (); ++iterX0, ++iterZ0) {
	const float x0 = *iterX0;
	const float z0 = *iterZ0;

	//
	list<float>::iterator iterX	= iterX0; ++iterX;
	list<float>::iterator iterZ	= iterZ0; ++iterZ;

	while ( iterX != xRes->end() ) {
	const float x = *iterX;
	const float z = *iterZ;
	if (fabs (x0 - x) < 1e-6 && fabs (z0 - z) < 1e-6) {
	iterX = xRes->erase (iterX);
	iterZ = zRes->erase (iterZ);
	} else {
	++iterX;
	++iterZ;
	}
	}
	}
*/

    // FINISH
    allPoints.remove_if (deleteAll);
    escPoints.remove_if (deleteAll);

    free (xPanel);
    free (tPanel);	

    return;

}

bool ITracer2D::checkTriangle (float curX, float curT, int ix, int iz, const int mode, float* xPanel, float* tPanel, float* xres, float* zres, float& dist) {

    const int   ind0 = (ix + mode) * zNum_ + iz;
    const float dpx0 = xPanel [ind0];
    const float dpt0 = tPanel [ind0];
    if (fabs (curX - dpx0) < 1e-6 && fabs (curT - dpt0) < 1e-6) { 
	*xres = xStart_ + (ix + mode) * xStep_;
	*zres = zStart_ + iz * zStep_;
	return true;
    }

    const int   ind1 = ix * zNum_ + iz + mode;
    const float dpx1 = xPanel [ind1];
    const float dpt1 = tPanel [ind1];
    if (fabs (curX - dpx1) < 1e-6 && fabs (curT - dpt1) < 1e-6) { 
	*xres = xStart_ + ix * xStep_;
	*zres = zStart_ + (iz + mode) * zStep_;
	return true;
    }

    const int   ind2 = ix * zNum_ + iz;
    const float dpx2 = xPanel [ind2];
    const float dpt2 = tPanel [ind2];

    if ( !isPointInsideTriangle (curX, curT, dpx0, dpt0, dpx1, dpt1, dpx2, dpt2) )
	return false;
		
    const float x13 = dpx0 - dpx2;
    const float x32 = dpx2 - dpx1;

    const float y23 = dpt1 - dpt2;
    const float y13 = dpt0 - dpt2;

    const float denom = y23 * x13 + x32 * y13; 
    if (fabs (denom) < 1e-9) return false;

    const float fX = curX - dpx2;
    const float fT = curT - dpt2;

    dist = sqrt (fX*fX + 9*fT*fT);  // average velocity - 3 km/s
    // x in m
    // t in ms -> 1ms ~ 3m

    const float w0 = (  y23 * fX + x32 * fT) / denom;
    const float w1 = ( -y13 * fX + x13 * fT) / denom;
    const float w2 = 1.f - w0 - w1;

    const float ipx0 = xStart_ + (ix + mode) * xStep_;
    const float ipy0 = zStart_ + iz * zStep_;
    const float ipx1 = xStart_ + ix * xStep_;
    const float ipy1 = zStart_ + (iz + mode) * zStep_;
    const float ipx2 = xStart_ + ix * xStep_;
    const float ipy2 = zStart_ + iz * zStep_;

    *xres = w0 * ipx0 + w1 * ipx1 + w2 * ipx2;
    *zres = w0 * ipy0 + w1 * ipy1 + w2 * ipy2;

    return true;
}
