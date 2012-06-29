#ifndef DEPTH_MIGRATOR_BASE_H
#define DEPTH_MIGRATOR_BASE_H

#include "support.hh"

class CurveDefinerBase;

struct WavefrontParams {
  
    int   rNum;
    float rStart;
    float rStep; 

    int   tNum;
    float tStart;
    float tStep; 
};

struct VelocityParams {
  
    int   zNum;
    float zStart;
    float zStep; 

    int   xNum;
    float xStart;
    float xStep; 
};

class WavefrontTracer {
public:
		WavefrontTracer ();
		~WavefrontTracer ();
		
		void getEscapePoints (float xSource, float zSource, EscapePoint* ep);

	
	void setParams (int raysNum, float raysStep, float raysStart) {wp_.rNum = raysNum; wp_.rStart = raysStart; wp_.rStep = raysStep;}
	void setVelModelParams (int zNum, float zStep, float zStart, 
							int xNum, float xStep, float xStart) {vp_.zNum = zNum; vp_.zStep = zStep; vp_.zStart = zStart;
																  vp_.xNum = xNum; vp_.xStep = xStep; vp_.xStart = xStart;}

	void setVelModel (float** velField) {velField_ = velField;}
	
    WavefrontParams    wp_;
    VelocityParams     vp_;

private:

    float** velField_;

};

class DepthMigratorBase {

public:

                 DepthMigratorBase ();
        virtual ~DepthMigratorBase ();

		virtual void processGather  (Point2D& curGatherCoords, const float* const data, float* image, float* gather, float* aCig);


		void setImagingParams (VolumeParams* dp, float* ptrToData, bool isAA, bool isCMP,
							   VolumeParams* vp, VolumeParams* ip, GatherParams* gp) {dp_ = dp; ptrToData_ = ptrToData; isAA_ = isAA; isCMP_ = isCMP;
																				  vp_ = vp; ip_ = ip; gp_ = gp;}

		float getMigVel (const float* const velTrace, const float curZ);
		void  initCurveDefiner (bool is3D);
		void setDataLimits ();

		void setVelModel (float** velField) {wavefrontTracer_.setVelModel (velField);} 	

		CurveDefinerBase* curveDefiner_;

		int             curOffset_;

		WavefrontTracer wavefrontTracer_;

protected:

 
 //	virtual float  getSampleFromData            (const float h, const float geoY, const float geoX, const float t, const float trf = 0.f) = 0;


		void getAzDipFromXY (float curDipY, float curDipX, float& curDip, float&curAz);

		bool isPointInsidePoly (Point2D* poly, int nPoly, Point2D& p0);


		VolumeParams*    dp_;
		GatherParams*    gp_;
		VolumeParams*      ip_;
		VolumeParams*      vp_;

		float*           ptrToData_;
	
		bool             isAA_;
		bool             isCMP_;

		// data limits
        float            dataZMin_;
        float            dataXMin_;
    	float            dataYMin_;

        float            dataZMax_;
        float            dataXMax_;
	    float            dataYMax_;

        int              dataZMaxSamp_;
        int              dataXMaxSamp_;
        int              dataYMaxSamp_;
};

#endif
