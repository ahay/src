#ifndef TIME_MIGRATOR_BASE_H
#define TIME_MIGRATOR_BASE_H

#include "support.h"

class CurveDefinerBase;

class TimeMigratorBase {

public:

                  TimeMigratorBase ();
                 ~TimeMigratorBase ();

		virtual void processGather  (Point2D& curGatherCoords, float curOffset, const float* const velTrace, const bool isAzDip,
  								     float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq);

		void setImagingParams (VolumeParams* dp, float* ptrToData, bool isAA, bool isCMP,
							   VolumeParams* vp, VolumeParams* ip, GatherParams* gp) {dp_ = dp; ptrToData_ = ptrToData; isAA_ = isAA; isCMP_ = isCMP;
																				  vp_ = vp; ip_ = ip; gp_ = gp;}

		float getMigVel (const float* const velTrace, const float curZ);
		void  initCurveDefiner (bool is3D);
		void setDataLimits ();

		CurveDefinerBase* curveDefiner_;

		int             curOffset_;

protected:

		virtual int   getSampleByBeam    (const float yCIG, const float xCIG, const float curZeroTime, const float curDip, 
						 				  const float curAz, const float migVel, const bool isAzDip, float& sample) = 0;

		virtual float getSampleByRay     (const float yCIG, const float xCIG, const float curZeroTime, 
					    				  const float curDip, const float curAz, const float migVel, const bool isAzDip,
									      const float yEmAngle, const float xEmAngle) = 0;

		virtual float getSampleFromData  (const float geoY, const float geoX, const float ti, const float p) = 0;

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
