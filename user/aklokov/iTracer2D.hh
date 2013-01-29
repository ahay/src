#ifndef ITRACER2D_H
#define ITRACER2D_H

class ITracer2D {

public:

	 ITracer2D ();
	~ITracer2D ();

	void init (int zNum, float zStart, float zStep, 
  			   int pNum, float pStart, float pStep,
			   int xNum, float xStart, float xStep);

	void  traceImage (float* xVol, float* tVol, float x0, float z0, float p0, float* xRes, float* zRes);

private: 

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
};
#endif
