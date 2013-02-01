#include "iTracer2D.hh"
#include <rsf.hh>
#include "support.hh"
#include <algorithm>
#include <list>
#include <iostream>
using namespace std;

bool pred (const ImagePoint2D* lhs, const ImagePoint2D* rhs)
{
  return lhs->x_ < rhs->x_;
}

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

void ITracer2D::traceImage2 (float* xVol, float* tVol, float x0, float z0, float p0, float* xRes, float* zRes) {

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

		list<ImagePoint2D*>::iterator iter1 = allPoints.begin ();
		list<ImagePoint2D*>::iterator iter2 = allPoints.end ();
		--iter2;

		bool found1 (false);
		list<ImagePoint2D*>::iterator iterMin;
		while (!found1 && iter1 != allPoints.end()) {
			ImagePoint2D* p1 = *iter1;
			if (p1->x_ > x1) { iterMin = iter1; found1 = true; }
			++iter1;
		}
		if (!found1) continue;

		bool found2 = false;
		list<ImagePoint2D*>::iterator iterMax;
		while ( !found2 && iter2 != allPoints.begin () ) {
			ImagePoint2D* p2 = *iter2;
			if (p2->x_ < x2) {iterMax = iter2; found2 = true; }
			--iter2;
		}
		if (!found2) continue;

		list<ImagePoint2D*>::iterator iter;
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

	free (ex);
	free (et);
	free (xPanel);
	free (tPanel);	

	return;

}

bool ITracer2D::checkTriangle (float curX, float curT, int ix, int iz, const int mode, float* xPanel, float* tPanel, float* xRes, float* zRes) {

	Point2D dataPoint [3];	
	int ind0 = (ix + mode) * zNum_ + iz;
	dataPoint[0].setX (xPanel [ind0]);
	dataPoint[0].setY (tPanel [ind0]);
	int ind1 = ix * zNum_ + iz + mode;
	dataPoint[1].setX (xPanel [ind1]);
	dataPoint[1].setY (tPanel [ind1]);
	int ind2 = ix * zNum_ + iz;
	dataPoint[2].setX (xPanel [ind2]);
	dataPoint[2].setY (tPanel [ind2]);

	if ( isPointInsideTriangle (curX, curT,
								dataPoint[0].getX (), dataPoint[0].getY (), 
								dataPoint[1].getX (), dataPoint[1].getY (),
								dataPoint[2].getX (), dataPoint[2].getY ()) ) {
		
		const float x13 = dataPoint[0].getX () - dataPoint[2].getX ();
		const float x32 = dataPoint[2].getX () - dataPoint[1].getX ();

		const float y23 = dataPoint[1].getY () - dataPoint[2].getY ();
		const float y13 = dataPoint[0].getY () - dataPoint[2].getY ();

		const float denom = y23 * x13 + x32 * y13; 

		float w0 = ( y23 * (curX - dataPoint[2].getX())  + x32 * (curT - dataPoint[2].getY()) ) / denom;
		float w1 = ( -y13 * (curX - dataPoint[2].getX()) + x13 * (curT - dataPoint[2].getY()) ) / denom;
		float w2 = 1 - w0 - w1;

		Point2D imagePoint [3];	
		imagePoint[0].setX (xStart_ + (ix + mode) * xStep_);
		imagePoint[0].setY (zStart_ + iz * zStep_);
		imagePoint[1].setX (xStart_ + ix * xStep_);
		imagePoint[1].setY (zStart_ + (iz + mode) * zStep_);
		imagePoint[2].setX (xStart_ + ix * xStep_);
		imagePoint[2].setY (zStart_ + iz * zStep_);

		float x = w0 * imagePoint[0].getX () + w1 * imagePoint[1].getX () + w2 * imagePoint[2].getX ();
		float z = w0 * imagePoint[0].getY () + w1 * imagePoint[1].getY () + w2 * imagePoint[2].getY ();						

		*xRes  = x;
		*zRes  = z;

		return true;
	}

	return false;
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

	float* pXPanel = xPanel;
	float* pTPanel = tPanel;
	for (int ix = 0; ix < xNum_; ++ix) {
		for (int iz = 0; iz < zNum_; ++iz, ++pXPanel, ++pTPanel) {
			const int pind = ix * pNum_ * zNum_ + pInd * zNum_ + iz;
			*pXPanel = xVol [pind];
			*pTPanel = 2 * tVol [pind];
		}
	}

	const int xNumRed = xNum_ - 1;
	const int zNumRed = zNum_ - 1;

	for (int ix = 0; ix < xNumRed; ++ix) {
		for (int iz = 0; iz < zNumRed; ++iz) {

			// upper triangle

			Point2D dataPoint [3];	
			int ind0 = (ix + 1) * zNum_ + iz;
			dataPoint[0].setX (xPanel [ind0]);
			dataPoint[0].setY (tPanel [ind0]);
			int ind1 = ix * zNum_ + iz + 1;
			dataPoint[1].setX (xPanel [ind1]);
			dataPoint[1].setY (tPanel [ind1]);
			int ind2 = ix * zNum_ + iz;
			dataPoint[2].setX (xPanel [ind2]);
			dataPoint[2].setY (tPanel [ind2]);

			// go along data line

			for (int ip = 0; ip < pNum_; ++ip) {

				float curT = et[ip];
				float curX = ex[ip];

				Point2D p0 (curX, curT);

				if ( isPointInsideTriangle (curX, curT,
											dataPoint[0].getX (), dataPoint[0].getY (), 
											dataPoint[1].getX (), dataPoint[1].getY (),
											dataPoint[2].getX (), dataPoint[2].getY ()) ) {
		
					const float x13 = dataPoint[0].getX () - dataPoint[2].getX ();
					const float x32 = dataPoint[2].getX () - dataPoint[1].getX ();
	
					const float y23 = dataPoint[1].getY () - dataPoint[2].getY ();
					const float y13 = dataPoint[0].getY () - dataPoint[2].getY ();

					const float denom = y23 * x13 + x32 * y13; 

					float w0 = ( y23 * (curX - dataPoint[2].getX())  + x32 * (curT - dataPoint[2].getY()) ) / denom;
					float w1 = ( -y13 * (curX - dataPoint[2].getX()) + x13 * (curT - dataPoint[2].getY()) ) / denom;
					float w2 = 1 - w0 - w1;

					Point2D imagePoint [3];	
					imagePoint[0].setX (xStart_ + (ix + 1) * xStep_);
					imagePoint[0].setY (zStart_ + iz * zStep_);
					imagePoint[1].setX (xStart_ + ix * xStep_);
					imagePoint[1].setY (zStart_ + (iz + 1) * zStep_);
					imagePoint[2].setX (xStart_ + ix * xStep_);
					imagePoint[2].setY (zStart_ + iz * zStep_);


					float x = w0 * imagePoint[0].getX () + w1 * imagePoint[1].getX () + w2 * imagePoint[2].getX ();
					float z = w0 * imagePoint[0].getY () + w1 * imagePoint[1].getY () + w2 * imagePoint[2].getY ();						

					xRes [ip] = x;
					zRes [ip] = z;
				}
			}
	
			// lower triangle

			ind0 = (ix + 1) * zNum_ + iz;
			dataPoint[0].setX (xPanel [ind0]);
			dataPoint[0].setY (tPanel [ind0]);
			ind1 = ix * zNum_ + iz + 1;
			dataPoint[1].setX (xPanel [ind1]);
			dataPoint[1].setY (tPanel [ind1]);
			ind2 = (ix + 1) * zNum_ + iz + 1;
			dataPoint[2].setX (xPanel [ind2]);
			dataPoint[2].setY (tPanel [ind2]);

			// go along data line

			for (int ip = 0; ip < pNum_; ++ip) {
				float curT = et[ip];
				float curX = ex[ip];

				Point2D p0 (curX, curT);

				if ( isPointInsideTriangle (curX, curT,
											dataPoint[0].getX (), dataPoint[0].getY (), 
											dataPoint[1].getX (), dataPoint[1].getY (),
											dataPoint[2].getX (), dataPoint[2].getY ()) ) {
					
					const float x13 = dataPoint[0].getX () - dataPoint[2].getX ();
					const float x32 = dataPoint[2].getX () - dataPoint[1].getX ();
	
					const float y23 = dataPoint[1].getY () - dataPoint[2].getY ();
					const float y13 = dataPoint[0].getY () - dataPoint[2].getY ();

					const float denom = y23 * x13 + x32 * y13; 

					float w0 = ( y23 * (curX - dataPoint[2].getX())  + x32 * (curT - dataPoint[2].getY()) ) / denom;
					float w1 = ( -y13 * (curX - dataPoint[2].getX()) + x13 * (curT - dataPoint[2].getY()) ) / denom;
					float w2 = 1 - w0 - w1;

					Point2D imagePoint [3];	
					imagePoint[2].setX (xStart_ + ix * xStep_);
					imagePoint[2].setY (zStart_ + iz * zStep_);
					imagePoint[0].setX (xStart_ + (ix + 1) * xStep_);
					imagePoint[0].setY (zStart_ + iz * zStep_);
					imagePoint[1].setX (xStart_ + ix * xStep_);
					imagePoint[1].setY (zStart_ + (iz + 1) * zStep_);

					float x = w0 * imagePoint[0].getX () + w1 * imagePoint[1].getX() + w2 * imagePoint[2].getX();
					float z = w0 * imagePoint[0].getY () + w1 * imagePoint[1].getY() + w2 * imagePoint[2].getY();						

					xRes [ip] = x;
					zRes [ip] = z;
				}
			}
		}
	}

	free (ex);
	free (et);
	free (xPanel);
	free (tPanel);	

	return;
}
