#ifndef DEPTH_BACKFIT_MIGRATOR_2D_H
#define DEPTH_BACKFIT_MIGRATOR_2D_H

#include "iTracer2D.hh"

class DepthBackfitMigrator2D {

public:

	 DepthBackfitMigrator2D ();
	~DepthBackfitMigrator2D ();

	void init (int zNum, float zStart, float zStep, 
  			   int pNum, float pStart, float pStep,
			   int xNum, float xStart, float xStep,
			   int rNum, float rStart, float rStep,
			   float* xVol, float* tVol, bool isAA);

	void processPartialImage (float* piData, float curP, float* piImage);

private:

	bool getSample (float* data, const float curX, const float curZ, const float curP, float &sample);
	bool getSampleFromImage (float* data, const float curX, const float curZ, const float curP, float &sample);
	void getImageSample (float* piData, float curX, float curZ, float curP, float* sample);
	void processData    (float* piData); 

	void applyCasualIntegration (float *trace, int n);
	void applyAnticasualIntegration (float *trace, int n);

	ITracer2D iTracer_;

	float* xVol_;
	float* tVol_;

	int   zNum_;
	float zStep_;
	float zStart_;

	int   pNum_;
	float pStep_;
	float pStart_;

	int   xNum_;
	float xStep_;
	float xStart_;

	int   rNum_;
	float rStep_;
	float rStart_;
	
	bool  isAA_;

};
#endif
