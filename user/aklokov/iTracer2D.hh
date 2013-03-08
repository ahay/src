#ifndef ITRACER2D_H
#define ITRACER2D_H

#include <list>
using namespace std;

class ImagePoint2D {
public:

	 ImagePoint2D ();
	 ImagePoint2D (float x, float z, int ix, int iz);
	~ImagePoint2D ();

	 ImagePoint2D& operator= (const ImagePoint2D& point);
  //   bool operator < (const ImagePoint2D& p) const;

	float x_;
	float z_;		
	int   ix_;
	int   iz_;
};

class ITracer2D {

public:

	 ITracer2D ();
	~ITracer2D ();

	void init (int zNum, float zStart, float zStep, 
  			   int pNum, float pStart, float pStep,
			   int xNum, float xStart, float xStep,
			   float dx, float dt);

	void  traceImage (float* xVol, float* tVol, float x0, float z0, float p0, float sa, list<float>* xRes, list<float>* zRes);

private: 

	bool  checkTriangle (float curX, float curT, int ix, int iz, const int mode, float* xPanel, float* tPanel, float* xRes, float* zRes, float& dist);

	bool  isPointInsideTriangle (float x0, float y0, float x1, float y1, 
								 float x2, float y2, float x3, float y3);

	int   zNum_;
	float zStep_;
	float zStart_;

	int   pNum_;
	float pStep_;
	float pStart_;

	int   xNum_;
	float xStep_;
	float xStart_;

	float dx_;
	float dt_;
};
#endif
