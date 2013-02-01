#include "iTracer2D.hh"
#include <rsf.hh>
#include "support.hh"
#include <algorithm>
#include <list>
#include <iostream>
using namespace std;

static bool deleteAll (ImagePoint2D* p) { delete p; return true; }
static bool pred      (const ImagePoint2D* lhs, const ImagePoint2D* rhs) { return lhs->x_ < rhs->x_; }

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
					  int xNum, float xStart, float xStep) {
	
	zNum_   = zNum;
	zStep_  = zStep;
	zStart_ = zStart;

	pNum_   = pNum;
	pStep_  = pStep;
	pStart_ = pStart;

	xNum_   = xNum;
	xStep_  = xStep;
	xStart_ = xStart;	

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

void ITracer2D::traceImage (float* xVol, float* tVol, float x0, float z0, float p0, float* xRes, float* zRes) {

	const int zInd = (z0 - zStart_) / zStep_;
	const int pInd = (p0 - pStart_) / pStep_;
	const int xInd = (x0 - xStart_) / xStep_;

	// escape line for the diffraction point
	float* et = sf_floatalloc (pNum_);
	memset ( et, 0, pNum_ * sizeof (float) );
	float* ex = sf_floatalloc (pNum_);
	memset ( ex, 0, pNum_ * sizeof (float) );

	for (int ip = 0; ip < pNum_; ++ip) {
		const int pind = xInd * pNum_ * zNum_ + ip * zNum_ + zInd;
		ex[ip] = xVol [pind];
		et[ip] = 2 * tVol [pind];
	}

	// constant-dip panel extraction
	const int panelSize = zNum_ * xNum_;	
	float* xPanel = sf_floatalloc (panelSize);
	float* tPanel = sf_floatalloc (panelSize);
	memset ( xPanel, 0, panelSize * sizeof (float) );
	memset ( tPanel, 0, panelSize * sizeof (float) );

	list<ImagePoint2D*> allPoints;

	float* pXPanel = xPanel;
	float* pTPanel = tPanel;
	for (int ix = 0; ix < xNum_; ++ix) {
		for (int iz = 0; iz < zNum_; ++iz, ++pXPanel, ++pTPanel) {
			const int pind = ix * pNum_ * zNum_ + pInd * zNum_ + iz;
			*pXPanel = xVol [pind];
			*pTPanel = 2 * tVol [pind];

		    ImagePoint2D* p = new ImagePoint2D (*pXPanel, *pTPanel, ix, iz);
			allPoints.push_back (p);
		}
	}

	// CHANGE ME
	// sorting by x_
	allPoints.sort (pred);

	const float dx = 20;
	const float dt = 5;

	// loop over escape points
	for (int ip = 0; ip < pNum_; ++ip) {

//		sf_warning ("%d", ip);

		float curT = et[ip];
		float curX = ex[ip];

		float x1 = curX - dx;
		float x2 = curX + dx;

		list<ImagePoint2D*>::iterator iter = allPoints.begin ();

		bool found1 (false);
		list<ImagePoint2D*>::iterator iterMin;
		while (!found1 && iter != allPoints.end()) {
			ImagePoint2D* p1 = *iter;
			if (p1->x_ > x1) { iterMin = iter; found1 = true; }
			++iter;
		}
		if (!found1) continue;

		bool found2 = false;
		list<ImagePoint2D*>::iterator iterMax;
		while ( !found2 && iter != allPoints.end () ) {
			ImagePoint2D* p2 = *iter;
			if (p2->x_ > x2) {iterMax = iter; found2 = true; }
			++iter;
		}
		if (!found2) continue;

		bool isOk = false;
		for (iter = iterMin; iter != iterMax && !isOk; ++iter) {
			ImagePoint2D* iPoint = *iter;			
			int ix = iPoint->ix_;			
			int iz = iPoint->iz_;			

			if (ix && iz) {	
			    const int mode = -1; // upper triangle
			    isOk = this->checkTriangle (curX, curT, ix, iz, mode, xPanel, tPanel, xRes + ip, zRes + ip);
			}
			if (ix < xNum_ - 1 && iz < zNum_ - 1) {
				const int mode = 1; // lower triangle
				isOk = this->checkTriangle (curX, curT, ix, iz, mode, xPanel, tPanel, xRes + ip, zRes + ip);
			}
		}
	}		

	// FINISH
	allPoints.remove_if (deleteAll);

	free (ex);
	free (et);
	free (xPanel);
	free (tPanel);	

	return;

}

bool ITracer2D::checkTriangle (float curX, float curT, int ix, int iz, const int mode, float* xPanel, float* tPanel, float* xres, float* zres) {

	const int   ind0 = (ix + mode) * zNum_ + iz;
	const float dpx0 = xPanel [ind0];
	const float dpt0 = tPanel [ind0];
	const int   ind1 = ix * zNum_ + iz + mode;
	const float dpx1 = xPanel [ind1];
	const float dpt1 = tPanel [ind1];
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
