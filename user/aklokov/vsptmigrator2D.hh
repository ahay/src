#ifndef VSP_TIME_MIGRATOR_2D_H
#define VSP_TIME_MIGRATOR_2D_H

#include "dmigratorBase.hh"

class VSPTimeMigrator2D : public DepthMigratorBase {

public:

	    VSPTimeMigrator2D ();
	   ~VSPTimeMigrator2D ();

    void  processGather (Point2D& curGatherCoords, const float* const data, float* image, float* dag, float* aCig, float* mCig, 
			 float* xEsc, float* tEsc);

	void   setVelModel (float** velField) { velField_ = velField; } 	

private:

	bool   getSampleFromData (const float geoY, const float geoX, const float t, const float p, float& sample);



	bool   getSampleByBeam   (float curX, float curZ, float vel, float curScatAngle, float curDipAngle, float& sample);

	float getSampleByRay (const float yCIG, const float xCIG, const float curZeroTime, 
					      const float curDip, const float curAz, const float migVel, const bool isAzDip,

					      const float yEmAngle, const float xEmAngle);

	bool getReceiverByDirection (float curX, float curZ, float dir2, float& rec);
	bool getSourceByDirection (float curX, float curZ, float dir2, float& src);

	void processTimeSample (const float curX, const float curZ, const float* const data, 
											float* curImage, float* curDag, float* curCig);

	float getVel (float curZ, float xCIG);



	// velocity model
	float**  velField_;

};
#endif
