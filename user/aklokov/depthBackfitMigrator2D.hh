#ifndef DEPTH_BACKFIT_MIGRATOR_2D_H
#define DEPTH_BACKFIT_MIGRATOR_2D_H

class DepthBackfitMigrator2D {

public:

	 DepthBackfitMigrator2D ();
	~DepthBackfitMigrator2D ();

	void processParialImage (float* piData, float curP, float* xVol, float* tVol, float* piImage);

};
#endif
