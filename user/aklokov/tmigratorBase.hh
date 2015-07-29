#ifndef TIME_MIGRATOR_BASE_H
#define TIME_MIGRATOR_BASE_H

#include "support.hh"

class CurveDefinerBase;

class TimeMigratorBase {

public:

                  TimeMigratorBase ();
                 virtual ~TimeMigratorBase ();

		virtual void processGather  (Point2D& curGatherCoords, float curOffset, const float* const velTrace, const bool isAzDip,
  								     float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq);

		void setImagingParams (VolumeParams* dp, float* ptrToData, bool isAA, int axis2label,
							   VolumeParams* vp, VolumeParams* ip, GatherParams* gp) { dp_ = dp; ptrToData_ = ptrToData; isAA_ = isAA;
																					   axis2label_ = axis2label; vp_ = vp; ip_ = ip; gp_ = gp; }

		float getMigVel (const float* const velTrace, const float curZ);
		void  initCurveDefiner (bool is3D);
		void setDataLimits ();

		virtual void getStackTaper (const float edgeTaper, const bool isDipAz) = 0;


		CurveDefinerBase* curveDefiner_;

		float             curOffset_;

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
		VolumeParams*    ip_;
		VolumeParams*    vp_;

		float*           ptrToData_;
		float*           stackTaper_;	

		bool             isAA_;
		bool             useRay_;
		int              axis2label_;

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
