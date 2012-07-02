#ifndef DEPTH_MIGRATOR_BASE_H
#define DEPTH_MIGRATOR_BASE_H

#include "support.hh"

struct WavefrontParams {
    int rNum; float rStart; float rStep; 
    int tNum; float tStart; float tStep; 
};
// velocity model parameters for wavefront-tracer
struct VelocityParams {
    int  zNum; float zStart; float zStep; 
    int  xNum; float xStart; float xStep; 
};

class WavefrontTracer {

public:
	  	       WavefrontTracer ();
  	          ~WavefrontTracer ();
		
		void   getEscapePoints    (float xSource, float zSource, EscapePoint* ep);
	
 		void   setParams          (int raysNum, float raysStep, float raysStart) { wp_.rNum = raysNum; wp_.rStart = raysStart; 
																				   wp_.rStep = raysStep; }
		void   setVelModelParams  (int zNum, float zStep, float zStart, 
							  	   int xNum, float xStep, float xStart) {vp_.zNum = zNum; vp_.zStep = zStep; vp_.zStart = zStart;
									 							         vp_.xNum = xNum; vp_.xStep = xStep; vp_.xStart = xStart;}
		void   setVelModel        (float** pVelField) { pVelField_ = pVelField; }
	
		int    getRaysNum         () { return wp_.rNum;   }
		float  getRaysStart       () { return wp_.rStart; }
		float  getRaysStep        () { return wp_.rStep;  }

private:

    	WavefrontParams    wp_;
	    VelocityParams     vp_;
		// pointer to velocity model
	    float**            pVelField_;
};

class DepthMigratorBase {

public:

                 DepthMigratorBase ();
        virtual ~DepthMigratorBase ();

		virtual void processGather  (Point2D& curGatherCoords, const float* const data, float* image, float* gather, float* aCig);

		virtual void setWavefrontTracerParams (int ttNum, float ttStep, float ttStart);

		void   setImagingParams   (VolumeParams* dp, float* ptrToData, bool isAA, bool isCMP,
							       VolumeParams* vp, VolumeParams* ip, GatherParams* gp) { dp_ = dp; ptrToData_ = ptrToData; isAA_ = isAA; 
									  													   isCMP_ = isCMP; vp_ = vp; ip_ = ip; gp_ = gp; }
		void   setVelModelParams  (int zNum, float zStep, float zStart, 
								   int xNum, float xStep, float xStart) { wavefrontTracer_.setVelModelParams ( zNum, zStep, zStart, 
																											   xNum, xStep, xStart); }
		void   setVelModel        (float** velField) { wavefrontTracer_.setVelModel (velField); } 	
		void   setDataLimits      ();

protected:

		WavefrontTracer  wavefrontTracer_; 

		VolumeParams*    dp_;
		GatherParams*    gp_;
		VolumeParams*    ip_;
		VolumeParams*    vp_;

		float*           ptrToData_;
	
		bool             isAA_;
		bool             isCMP_;

		// data limits
        float            dataTMin_;
        float            dataXMin_;
    	float            dataYMin_;

        float            dataTMax_;
        float            dataXMax_;
	    float            dataYMax_;

        int              dataTMaxSamp_;
        int              dataXMaxSamp_;
        int              dataYMaxSamp_;
};

#endif
