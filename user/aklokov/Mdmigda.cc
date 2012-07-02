/* 2D depth scattering-angle Kirchhoff migration  */
/*
  Copyright (C) 2012 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.hh>

#include "support.hh"
#include "dmigratorBase.hh"
#include "dmigrator2D.hh"

//#ifdef _OPENMP
//#include <omp.h>
//#endif

RunParamsTmigda  rp; // migration (running) params
VolumeParams     dp; // data params
VolumeParams     vp; // velocity params
VolumeParams     ip; // image params
GatherParams     gp; // gather params

// files
sf_file dataFile;
sf_file velFile;

sf_file imageFile;
sf_file dagFile;
sf_file acigFile;

//  Causal integration of a trace[n]
void applyCasualIntegration (float *trace, int n) {

    for (int i = 1; i < n; ++i)
        trace[i] += trace[i - 1];
        
    return;      
}

// Anticausal integrations of a trace[n]
void applyAnticasualIntegration (float *trace, int n) {

    for (int i = n - 2; i >= 0; --i)
        trace[i] += trace[i + 1];

    return;
}

void readData (float* data) {

    const int tracesNum = dp.xNum * dp.yNum * dp.hNum;
    const int dataSize = dp.zNum * tracesNum;
    memset (data, 0, dataSize * sizeof (float));

    const size_t startPos = 0;

    sf_seek (dataFile, startPos, SEEK_SET);
    sf_floatread (data, dataSize, dataFile);

    if (rp.isAA) {
		float* ptrData = data;
		for (int it = 0; it < tracesNum; ++it, ptrData += dp.zNum) {
		    applyCasualIntegration (ptrData, dp.zNum);
		    applyAnticasualIntegration (ptrData, dp.zNum);
		}
    }
	
    return;
}

void readVelocity (float** velModel) {

	for (int ix = 0; ix < vp.xNum; ++ix) {
		const size_t startPos = ix * vp.zNum * sizeof(float);

	    sf_seek (velFile, startPos, SEEK_SET);
	    sf_floatread (velModel[ix], vp.zNum, velFile);

	    if (!rp.isVelMS) {	// velocity is in km/s - we reduce it to m/s
		float* ptrVel = velModel[ix];
		for (int iv = 0; iv < vp.zNum; ++iv, ++ptrVel)
		    *ptrVel *= 1000;
		}
	}

	return;
}

void checkImageParams () {

    const float zVelMax = vp.zStart + vp.zStep * (vp.zNum - 1);
    const float xVelMax = vp.xStart + vp.xStep * (vp.xNum - 1);
    const float yVelMax = vp.yStart + vp.yStep * (vp.yNum - 1);

    // time checking

    if (vp.zStart - ip.zStart > 1E-6) {
	float diff = vp.zStart - ip.zStart;
	ip.zNum -= (diff / ip.zStep + 1);
	ip.zStart += (diff / ip.zStep + 1) * ip.zStep;
	sf_warning ("first migrated time was changed to %g", ip.zStart);			
    }
    if (ip.zNum <= 0) sf_error ("wrong image time sample number");	
    float zImageMax = ip.zStart + ip.zStep * (ip.zNum - 1);
    if (zImageMax - zVelMax > 1E-6) {
	float diff = zImageMax - zVelMax; 
	ip.zNum -= (diff / ip.zStep + 1);
	sf_warning ("time sample number was changed to %d", ip.zNum);	
    }
    if (ip.zNum <= 0) sf_error ("wrong image time sample number");	

    // inline checking

    if (vp.xStart - ip.xStart > 1E-6) {
	float diff = vp.xStart - ip.xStart;
	ip.xNum -= (diff / ip.xStep + 1);
	ip.xStart += (diff / ip.xStep + 1) * ip.xStep;
	sf_warning ("first inline coord was changed to %g", ip.xStart);			
    }
    if (ip.xNum <= 0) sf_error ("wrong image inline number");	
    float xImageMax = ip.xStart + ip.xStep * (ip.xNum - 1);
    if (xImageMax - xVelMax > 1E-6) {
	float diff = xImageMax - xVelMax; 
	ip.xNum -= (diff / ip.xStep + 1);
	sf_warning ("inline number was changed to %d", ip.xNum);	
    }
    if (ip.xNum <= 0) sf_error ("wrong image inline number");	

    // crossline checking

	if (!rp.is3D)
		return; // 2D migration is runnig

    if (vp.yStart - ip.yStart > 1E-6) {
	float diff = vp.yStart - ip.yStart;
	ip.yNum -= (diff / ip.yStep + 1);
	ip.yStart += (diff / ip.yStep + 1) * ip.yStep;
	sf_warning ("first crossline coord was changed to %g", ip.yStart);			
    }
    if (ip.yNum <= 0) sf_error ("wrong image crossline number");	
    float yImageMax = ip.yStart + ip.yStep * (ip.yNum - 1);
    if (yImageMax - yVelMax > 1E-6) {
	float diff = yImageMax - yVelMax; 
	ip.yNum -= (diff / ip.yStep + 1);
	sf_warning ("crossline number was changed to %d", ip.yNum);	
    }
    if (ip.yNum <= 0) sf_error ("wrong image crossline number");	
	
    return;
}

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc,argv);

// INPUT FILES
    dataFile = sf_input ("in");
    /* common-offset sections */

    if ( NULL != sf_getstring("vel") ) {
	/* velocity model file (velocity in m/s) */ 
	velFile  = sf_input ("vel");
	float firstvel;
	sf_floatread (&firstvel, 1, velFile);
	rp.isVelMS = true;		
	if (firstvel < 15) {
	    sf_warning ("it seems that velocity is in km/s - will be multiplied on 1000");	
	    rp.isVelMS = false;					
	}			

    } else { sf_error ("Need input: velocity model"); }

// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile) ) sf_error ("Need float input: common-offset sections");
    if ( SF_FLOAT != sf_gettype (velFile)  ) sf_error ("Need float input: velocity model");

// OUTPUT FILES

    imageFile  = sf_output("out");
    /* image */

    if ( NULL != sf_getstring("dag") ) {
	/* output file containing CIGs in the dip-angle domain */ 
	dagFile  = sf_output ("dag"); rp.isDag = true;
    } else { rp.isDag = false; }

    if ( NULL != sf_getstring("cig") ) {
	/* output file containing CIGs in the surface-offset domain */ 
	acigFile  = sf_output ("cig"); rp.isCig = true;
    } else { rp.isCig = false; }

    // data params

    if ( !sf_histint   (dataFile, "n1", &dp.zNum) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile, "d1", &dp.zStep) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile, "o1", &dp.zStart) ) sf_error ("Need o1= in input");

    if ( !sf_histint   (dataFile, "n2", &dp.xNum) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile, "d2", &dp.xStep) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile, "o2", &dp.xStart) ) sf_error ("Need o2= in input");

    if ( !sf_histint   (dataFile, "n3", &dp.yNum) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile, "d3", &dp.yStep) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile, "o3", &dp.yStart) ) sf_error ("Need o3= in input");

    if ( !sf_histint   (dataFile, "n4", &dp.hNum) )   sf_error ("Need n4= in input");
    if ( !sf_histfloat (dataFile, "d4", &dp.hStep) )  sf_error ("Need d4= in input");
    if ( !sf_histfloat (dataFile, "o4", &dp.hStart) ) sf_error ("Need o4= in input");

    //
    char* corUnit;
    char* unit;
	
    // data
    // time - in ms
    corUnit = (char*) "ms"; unit = sf_histstring (dataFile, "unit1"); if (!unit) sf_error ("unit1 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.zStep *= 1000; dp.zStart *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit2"); if (!unit) sf_error ("unit2 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.xStep *= 1000; dp.xStart *= 1000; }
    // crossline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit3"); if (!unit) sf_error ("unit3 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.yStep *= 1000; dp.yStart *= 1000; }
    // offset - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit4"); if (!unit) sf_error ("unit4 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.hStep *= 1000; dp.hStart *= 1000; }
	
    // VELOCITY MODEL PARAMS

    if ( !sf_histint   (velFile, "n1", &vp.zNum) )   sf_error ("Need n1= in velocity file");
    if ( !sf_histfloat (velFile, "d1", &vp.zStep) )  sf_error ("Need d1= in velocity file");
    if ( !sf_histfloat (velFile, "o1", &vp.zStart) ) sf_error ("Need o1= in velocity file");

    if ( !sf_histint   (velFile, "n2", &vp.xNum) )   sf_error ("Need n2= in velocity file");
    if ( !sf_histfloat (velFile, "d2", &vp.xStep) )  sf_error ("Need d2= in velocity file");
    if ( !sf_histfloat (velFile, "o2", &vp.xStart) ) sf_error ("Need o2= in velocity file");

    if ( !sf_histint   (velFile, "n3", &vp.yNum) )   sf_error ("Need n3= in velocity file");
    if ( !sf_histfloat (velFile, "d3", &vp.yStep) )  sf_error ("Need d3= in velocity file");
    if ( !sf_histfloat (velFile, "o3", &vp.yStart) ) sf_error ("Need o3= in velocity file");

    // time - in ms
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit1"); if (!unit) sf_error ("unit1 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.zStep *= 1000; vp.zStart *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit2"); if (!unit) sf_error ("unit2 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.xStep *= 1000; vp.xStart *= 1000; }
    // crossline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit3"); if (!unit) sf_error ("unit3 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.yStep *= 1000; vp.yStart *= 1000; }

// Migration parameters
    if (!sf_getbool ("is3d",       &rp.is3D))       rp.is3D = false;
    /* if y, apply 3D migration */
    if (!sf_getbool ("isCMP",      &rp.isCMP))      rp.isCMP = false;
    /* if y, data traces have coordinates of CMP position */
    if (!sf_getbool ("isAA",       &rp.isAA))       rp.isAA = true;
    /* if y, apply anti-aliasing */
    if (!sf_getbool ("isDipAz",    &rp.isDipAz))    rp.isDipAz = true;
    /* if y, apply dip/azimuth mode; if n, apply inline/crossline angle mode */

    // IMAGE PARAMS
    if (!sf_getint ("izn", &ip.zNum))        ip.zNum = dp.zNum;	
    /* number of imaged depth samples */
    if (!sf_getint ("ixn", &ip.xNum))        ip.xNum = dp.xNum;	
    /* number of imaged inlines */
    if (!sf_getint ("iyn", &ip.yNum))        ip.yNum = rp.is3D ? vp.yNum : 1;	
    /* number of imaged crosslines */
    if (!sf_getint ("iscatn", &gp.scatNum))   gp.scatNum = 1;	
    /* number of scattering-angles */
    if (!sf_getfloat ("izo", &ip.zStart))    ip.zStart = dp.zStart;
    /* first imaged depth (in meters) */
    if (!sf_getfloat ("ixo", &ip.xStart))    ip.xStart = dp.xStart;
    /* first imaged inline (in meters) */
    if (!sf_getfloat ("iyo", &ip.yStart))    ip.yStart = dp.yStart;	
    /* first imaged crossline (in meters) */
    if (!sf_getfloat ("iscato", &gp.scatStart)) gp.scatStart = 0;	
    /* first scattering-angle (in degree) */
    if (!sf_getfloat ("izd", &ip.zStep))     ip.zStep = dp.zStep;
    /* step in depth (in meters) */
    if (!sf_getfloat ("ixd", &ip.xStep))     ip.xStep = dp.xStep;	
    /* step in inlines (in meters) */
    if (!sf_getfloat ("iyd", &ip.yStep))     ip.yStep = dp.yStep;
    /* step in crosslines (in meters) */
	if (!sf_getfloat ("iscatd", &gp.scatStep)) gp.scatStep = 10.f;	
    /* scattering-angle increment (in degree) */

	// CHECK IMAGE PARAMETERS

    checkImageParams ();

    // GATHER PARAMS

    gp.zNum = ip.zNum;
    if (!sf_getint ("dipn" , &gp.dipNum))      gp.dipNum = 161;	
    /* number of dip-angles */
    if (!sf_getint ("sdipn", &gp.sdipNum))     gp.sdipNum = 1;	
    /* number of secondary (azimuth or crossline) angles */
    gp.zStart = ip.zStart;
    if (!sf_getfloat ("dipo",  &gp.dipStart))   gp.dipStart = -80.f;	
    /* first dip-angle */
    if (!sf_getfloat ("sdipo", &gp.sdipStart))  gp.sdipStart = 90.f;
    /* first secondary (azimuth or crossline) angle */
    gp.zStep = ip.zStep;
    if (!sf_getfloat ("dipd",  &gp.dipStep))   gp.dipStep = 1.f;	
    /* step in dip-angle */
    if (!sf_getfloat ("sdipd", &gp.sdipStep))  gp.sdipStep = 1.f;	
    /* step in secondary (azimuth or crossline) angle */

	// TRAVEL TIMES TABLES
	int ttRayNum (0); float ttRayStep (0.f); float ttRayStart (0.f);
	int ttNum (0); float ttStep (0.f); float ttStart (0.f);

	const float maxTime = dp.zStart + (dp.zNum - 1) * dp.zStep;
    if ( !sf_getfloat ("ttd",  &ttStep) ) ttStep = 0.002f; // not sure that this value is the optimal
    /* travel-times increment */
    if ( !sf_getint  ("ttn", &ttNum) ) ttNum = 0.001 * 0.5 * maxTime / ttStep + 1;
    /* travel-times number */
    if ( !sf_getint  ("ttrayn", &ttRayNum) ) ttRayNum =  2 * gp.dipNum + 1;
    /* travel-times rays number */
    if ( !sf_getfloat ("ttrayd",  &ttRayStep) ) ttRayStep = gp.dipStep / 2.f;
    /* travel-times rays increment */
    if ( !sf_getfloat ("ttrayo",  &ttRayStart) ) ttRayStart = gp.dipStart - ttRayStep + 180.f; 	
    /* travel-times rays start */
	// calculating ttStart we do "+180" because the direction to the top corresponds to 180 degree

    // INITIATE OUTPUT

    // image file
    sf_putint (imageFile, "n1", ip.zNum); sf_putint (imageFile, "n2", ip.xNum); sf_putint (imageFile, "n3", ip.yNum); sf_putint (imageFile, "n4", 1);
    sf_putfloat (imageFile, "d1", ip.zStep); sf_putfloat (imageFile, "d2", ip.xStep); sf_putfloat (imageFile, "d3", ip.yStep); 
    sf_putfloat (imageFile, "d4", 1);   
    sf_putfloat (imageFile, "o1", ip.zStart); sf_putfloat (imageFile, "o2", ip.xStart); sf_putfloat (imageFile, "o3", ip.yStart);    
    sf_putfloat (imageFile, "o4", 0);   
    sf_putstring(imageFile, "label1", "depth"); sf_putstring(imageFile, "label2", "inline"); sf_putstring(imageFile, "label3", "crossline");
    sf_putstring(imageFile, "unit1", "m"); sf_putstring(imageFile, "unit2", "m"); sf_putstring(imageFile, "unit3", "m");

    if (rp.isCig) {
		// angle gathers file
    	sf_putint (acigFile, "n1", ip.zNum); sf_putint (acigFile, "n2", gp.scatNum); sf_putint (acigFile, "n3", ip.xNum); 
		sf_putint (acigFile, "n4", ip.yNum); 
    	sf_putfloat (acigFile, "d1", ip.zStep); sf_putfloat (acigFile, "d2", gp.scatStep); sf_putfloat (acigFile, "d3", ip.xStep);
		sf_putfloat (acigFile, "d4", ip.yStep); 
    	sf_putfloat (acigFile, "o1", ip.zStart); sf_putfloat (acigFile, "o2", gp.scatStart); sf_putfloat (acigFile, "o3", ip.xStart);
		sf_putfloat (acigFile, "o4", ip.yStart);    
		sf_putstring(acigFile, "label1", "depth"); sf_putstring(acigFile, "label2", "scattering angle"); 
		sf_putstring(acigFile, "label3", "inline"); sf_putstring(acigFile, "label4", "crossline");
		sf_putstring(acigFile, "unit1", "m"); sf_putstring(acigFile, "unit2", "deg"); 
		sf_putstring(acigFile, "unit3", "m"); sf_putstring(acigFile, "unit4", "m");
    }

    if (rp.isDag) {
		// dip-angle gathers file
		sf_putint (dagFile, "n1", ip.zNum); sf_putint (dagFile, "n2", gp.dipNum); sf_putint (dagFile, "n3", gp.sdipNum);
		sf_putint (dagFile, "n4", ip.xNum); sf_putint (dagFile, "n5", ip.yNum);
    	sf_putfloat (dagFile, "d1", ip.zStep); sf_putfloat (dagFile, "d2", gp.dipStep); sf_putfloat (dagFile, "d3", gp.sdipStep);
		sf_putfloat (dagFile, "d4", ip.xStep); sf_putfloat (dagFile, "d5", ip.yStep);    
    	sf_putfloat (dagFile, "o1", ip.zStart); sf_putfloat (dagFile, "o2", gp.dipStart); sf_putfloat (dagFile, "o3", gp.sdipStart);
		sf_putfloat (dagFile, "o4", ip.xStart); sf_putfloat (dagFile, "o5", ip.yStart);    
		sf_putstring(dagFile, "label1", "depth");
		if (rp.isDipAz) {
		    sf_putstring(dagFile, "label2", "dip angle"); sf_putstring(dagFile, "label3", "azimuth");
		} else {
		    sf_putstring(dagFile, "label2", "inline slope"); sf_putstring(dagFile, "label3", "crossline slope");
		}
		sf_putstring(dagFile, "label4", "inline"); 	sf_putstring(dagFile, "label5", "crossline");
		sf_putstring(dagFile, "unit1", "m"); sf_putstring(dagFile, "unit2", "deg"); sf_putstring(dagFile, "unit3", "deg");
		sf_putstring(dagFile, "unit4", "m"); sf_putstring(dagFile, "unit5", "m");
    }

	// SIZES

    // dip-angle gather size
    const int dataSize = dp.zNum * dp.xNum * dp.yNum * dp.hNum;
	// dip-angle gather size
    const int dagSize  = gp.zNum * gp.dipNum;  // * gp.sdipNum; - for 3D migration
	// scattering-angle-gather size
    const int acigSize = gp.zNum * gp.scatNum;

	// MEMORY ALLOCATION

    // data
    float* data = sf_floatalloc (dataSize);
    // velocity model
    float** velModel = sf_floatalloc2 (vp.zNum, vp.xNum);
    // dip-angle gather
    float* dag  = sf_floatalloc (dagSize);
	// scattering-angle gather    
	float* acig = sf_floatalloc (acigSize);
	// image
	float* image = sf_floatalloc (gp.zNum);

	// DEFINE MIGRATOR

    // set migrator
    DepthMigratorBase* migrator;
//    if (rp.is3D) {
//		migrator = new DepthMigrator3D ();
//		migrator->initCurveDefiner (true);
//    } else {
		migrator = new DepthMigrator2D ();
//    }
    migrator->setImagingParams (&dp, data, rp.isAA, rp.isCMP, &vp, &ip, &gp);
    migrator->setDataLimits ();
	migrator->setWavefrontTracerParams (ttRayNum, ttRayStep, ttRayStart, ttNum, ttStep, ttStart);
	migrator->setVelModelParams ( vp.zNum, vp.zStep, vp.zStart,
 							      vp.xNum, vp.xStep, vp.xStart );
	// read data
	readData     (data);
	readVelocity (velModel);	
    migrator->setVelModel (velModel);

	// MAIN LOOP

    const int fullGatherNum = ip.yNum * ip.xNum;

    for (int ipy = 0; ipy < ip.yNum; ++ipy) {
		for (int ipx = 0; ipx < ip.xNum; ++ipx) {
			sf_warning ("gather %d of %d;", ipx + 1, fullGatherNum);	
			Point2D curGatherPos (ipx, ipy);
		    const size_t startInd = ipx + ipy * ip.xNum;

			memset (dag,   0, dagSize  * sizeof (float));
			memset (acig,  0, acigSize * sizeof (float));
			memset (image, 0, gp.zNum  * sizeof (float));

			migrator->processGather (curGatherPos, data, image, dag, acig);

			// output the image trace
			size_t startPos = startInd * gp.zNum * sizeof(float);
			sf_seek (imageFile, startPos, SEEK_SET);
			sf_floatwrite (image, gp.zNum, imageFile);	
			// output the dip-angle gather
		    if (rp.isDag) {
				startPos = startInd * dagSize * sizeof(float);
				sf_seek (dagFile, startPos, SEEK_SET);
				sf_floatwrite (dag, dagSize, dagFile);
		    }			
			// output the scattering-angle gather
			if (rp.isCig) {
			    // write trace into angle CIG file
			    startPos = startInd * acigSize * sizeof(float);
			    sf_seek (acigFile, startPos, SEEK_SET);
			    sf_floatwrite (acig, acigSize, acigFile);
			}
		}
	}

    free (dag);
    free (acig);
    free (image);

    free (*velModel);
    free (data);

    sf_fileclose (dagFile);
    sf_fileclose (acigFile);
    sf_fileclose (imageFile);

    sf_fileclose (velFile);
	sf_fileclose (dataFile);

    return 0;
}
