#ifndef DEPTH_MIGRATOR_2D_H
#define DEPTH_MIGRATOR_2D_H

#include "dmigratorBase.hh"

class DepthMigrator2D : public DepthMigratorBase {

public:

	       DepthMigrator2D ();
	      ~DepthMigrator2D ();

	void   processGather             (Point2D& curGatherCoords, const float* const data, float* image, float* gather, float* aCig, float* mCig,
									  float* xEsc, float* tEsc);
	void   processDepthSample        (const float curX, const float curZ, const float* const data, double* image, double* dag, double* aCig,
									  float* curMCig, float* xEsc, float* tEsc);
	void   calcTravelTimes           (float curZ, float curX, EscapePoint* escPoints);
	// ray tracing functions
	void   getEscPointByDirection    (EscapePoint* travelTimes, const float pRec, EscapePoint& resEscPoint);
	// calculate ray touching the current receiver
	bool   getRayToPoint             (EscapePoint* travelTimes, float curRecPos, float dir1, float dir2, float& timeToRec, float& recAbsP, bool& full);	
	// get sample by two-rays beam
	bool   getSampleByBeam  	     (EscapePoint* travelTimes, float curScatAngle, float curDipAngle, float& sample); 
	// get sample by only ray trace; implemented for zero-offset section only
	bool   getSampleByRay            (EscapePoint* travelTimes, float dipAngle, float& sample);
  
 	bool   getSampleFromData         (const float h, const float geoY, const float geoX, const float t, const float p, float& sample);
	// transfer parameters to wavefrontTracer
	void   setWavefrontTracerParams  (int ttRayNum, float ttRayStep, float ttRayStart, int ttNum, float ttStep, float ttStart);
	// setup velocity model
	void   setVelModel               (float** velField) { wavefrontTracer_.setVelModel (velField); velField_ = velField; } 	
	void   setWavefrontTracerAxes    ();

private:

	float  getVel                    (float curZ, float xCIG);

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
