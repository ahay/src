#ifndef DEPTH_MIGRATOR_2D_H
#define DEPTH_MIGRATOR_2D_H

#include "dmigratorBase.hh"

class DepthMigrator2D : public DepthMigratorBase {

public:

	       DepthMigrator2D ();
	      ~DepthMigrator2D ();

	void   processGather             (Point2D& curGatherCoords, const float* const data, float* image, float* gather, float* aCig);
	void   calcTravelTimes           (float curZ, float curX, EscapePoint* escPoints);
	// ray tracing functions
	void   getEscPointByDirection    (EscapePoint* escPoints, int size, float pRec, EscapePoint& resEscPoint);
	// calculate ray touching the current receiver
	int    getRayToPoint             (float curRecPos, float dir1, float dir2, float& timeToRec, float& recAbsP, bool& full);	
	// get sample by two-rays beam
	int    getSampleByBeam  	     (float curScatAngle, float curDipAngle, float& sample); 
	// get sample by only ray trace; implemented for zero-offset section only
	void   getSampleByRay            (float dipAngle, float& sample);
  
 	float  getSampleFromData         (const float h, const float geoY, const float geoX, const float t, const float trf = 0.f);
	// transfer parameters to wavefrontTracer
	void   setWavefrontTracerParams  (int ttRayNum, float ttRayStep, float ttRayStart, int ttNum, float ttStep, float ttStart);

private:
	
	EscapePoint* travelTimes_;

	// velocity model
	float**  velField_;
	// travel-time-tables parameters
	int      ttRayNum_;
	float    ttRayStart_;
	float    ttRayStep_;
	// used by ray tracing; define an angle corridor
	float    startDirMin_;
	float    startDirMax_;
}; 

#endif
